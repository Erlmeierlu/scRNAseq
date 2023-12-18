library(SingleR)
library(celldex)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(monocle3)
library(Seurat)
library(foreach)
library(scRNAseq)
library(data.table)
library(readr)
library(forcats)
library(scran)
library(readr)
library(biomaRt)

#personal theme
theme_my <- function() {
  
  theme(
    panel.grid.major = element_line(colour = "lightgray"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    ),
    panel.spacing.x = unit(10,"mm"),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.key = element_blank(),
    text = element_text(size = 15),
    strip.text.x = element_text(size = 10, margin = margin(b = 2, t = 2)),
    strip.background = element_rect(fill = "#9FD7D2", colour = "black", size = 1)
  )
}

#setting up directories
gfsDir <- '/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin'
plotsDir <- file.path(gfsDir, 'plots')
tablesDir <- file.path(gfsDir, 'tables')
oldDir <- "/vscratch/scRNAseq/data/old"
shinyDir <- 'dge-app'
dataDir <-"data"
resDir <- "results"

# Load Data ---------------------------------------------------------------

monocle.obj <- read_rds(file.path(dataDir, "2_unsupervised.cds"))

reference.obj <- readRDS(file.path(oldDir, "TwoGroups_celltypes_group.rds"))
reference.obj <- UpdateSeuratObject(object = reference.obj)
reference.obj <- as.SingleCellExperiment(reference.obj)
colnames(colData(reference.obj))[colnames(colData(reference.obj)) == "celltype"] <- "label.fine"

old.mon <- readRDS(file.path(oldDir, "scRNAseq_3_monocle_more_hash_cutoff.cds"))
old_meta <- as.data.frame(colData(old.mon))
old_counts <- assay(old.mon)
old.seurat <- CreateSeuratObject(counts = old_counts,
                                    project = "ol",
                                    assay = "counts",
                                    meta.data = old_meta)
old.mon <-  as.SingleCellExperiment(old.seurat)

colnames(colData(old.mon))[colnames(colData(old.mon)) == "celltype"] <- "label.fine"

# SingleR annotation ------------------------------------------------------

## Reference data ----------------------------------------------------------

### making a list of reference data sets  ----------------------------------------------

ref_data <- list(
  #Human Data Sets: 
  ##HPCA generally / Keratinocytes
  # hpca = HumanPrimaryCellAtlasData(),
  # 
  ##DB ImmuneCells, comprehensive CD4+ subsets; only one B cell subset, no dendritic cells
  dice = DatabaseImmuneCellExpressionData(),
  # 
  # ##Monaco Immunc cells
  # monaco = MonacoImmuneData(),

  #MouseData:
  ##MouseRNAseqData
  mrsd = MouseRNAseqData(),
  
  ##ImmGen 
  immgen = ImmGenData(),
  
  ##reference Dataset
  struc = reference.obj,
  
  ##old dataset
  old = old.mon
)

convert_human_to_mouse <- function(gene_list) {
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  setDT(mouse_human_genes)
  mouse_human_genes <- mouse_human_genes[Common.Organism.Name %in% c('human', "mouse, laboratory")]
  
  vapply(gene_list, function(x) {
    class_key = mouse_human_genes[Symbol == x, DB.Class.Key]
    if(identical(class_key, integer(0))) return(x)
    mouse_human_genes[DB.Class.Key %in% class_key & Common.Organism.Name == "mouse, laboratory", Symbol][1]
  },  FUN.VALUE = character(1), USE.NAMES = F)
}

#Human Gene names as Mouse Gene names
# rownames(ref_data$hpca) <- str_to_title(rownames(ref_data$hpca))
rownames(ref_data$dice) <- convert_human_to_mouse(rownames(ref_data$dice))
# rownames(ref_data$monaco) <- str_to_title(rownames(ref_data$monaco))

metadata <- as.data.frame(colData(monocle.obj))

counts_matrix <- assay(monocle.obj)

obj.as.seurat <- CreateSeuratObject(counts = counts_matrix,
                                project = "integrated_scRNAseq",
                                assay = "integrated",
                                meta.data = metadata)


sce <-  as.SingleCellExperiment(obj.as.seurat)



### run singleR with additional saving of original results data -------------------------------
# parallel computation
doParallel::registerDoParallel(cores = 5)
foreach(ref = names(ref_data)) %dopar% {
  message(ref)
  for(labelx in c("label.main", "label.fine")){
    message(paste(".", labelx))
    
    ref.file <- paste0("cell_types_", ref, "_", labelx, ".csv")
    
    ref.file.orig <- paste0("orig_cell_types_", ref, "_", labelx, ".rds")
    
    if(!labelx %in% colnames(colData(ref_data[[ref]]))) next
    
    message(paste(".", "running SingleR"))
    
    results <- SingleR(
      test = sce,
      ref = ref_data[[ref]],
      labels = colData(ref_data[[ref]])[, labelx],   
      de.method="wilcox"
    )
    
    # save results
    write_rds(results, file.path(vDir, "data", ref.file.orig))
    
    res <- data.table(
      as_tibble(results, rownames = "cell"),
      ref = ref,
      labels = labelx
    )
    
    #
    colnames(res) <- gsub("\\.", "_", colnames(res))
    
    res <- res[,c("cell", "labels", "tuning_scores_first", "tuning_scores_second"), with=F]
    
    for(cx in colnames(res)){
      if(is.numeric(res[[cx]])) res[[cx]] <- round(res[[cx]], 2)
    }
    
    write_csv(res, file.path(tablesDir, ref.file))
  }
  TRUE
}

### read in singleR results ------------------------------------------------

# read in the results

for (ref in names(ref_data)){

  for(labelx in c("label.main", "label.fine")){
      
    if(!labelx %in% colnames(colData(ref_data[[ref]]))) next
    assign(paste0("cell_types_", ref, "_", labelx), 
           read.csv(file.path(tablesDir, paste0("cell_types_", ref, "_", labelx, ".csv"))))
    
    assign(paste0("orig_cell_types_", ref, "_", labelx), 
           readRDS(file.path(vDir, "data", paste0("orig_cell_types_", ref, "_", labelx, ".rds"))))
    
  }
  
}

### assign labels to colData ----------------------------------------------------

# ref based labels
for(variable in ls(pattern = "^cell_types")){
  colData(monocle.obj)[[variable]] <- get(variable)$labels
}

# monocle.obj@colData$celltype_raw <- as.factor(
#   case_when(
#     monocle.obj@colData$cell_types_immgen_label.main == "Epithelial cells" ~ "Keratinocytes",
#     TRUE ~ monocle.obj@colData$cell_types_immgen_label.main
#   )
# )

monocle.obj@colData$celltype_raw <- monocle.obj@colData$cell_types_old_label.fine

ct_majority <-
  monocle.obj@colData %>% as_tibble() %>% 
  group_by(Cluster, celltype_raw) %>%  
  summarize(n = n()) %>% 
  slice_max(n = 1, order_by = n, with_ties = F)

monocle.obj@colData$celltype <-
  str_replace_all(ct_majority[match(monocle.obj@colData$Cluster, ct_majority$Cluster),]$celltype_raw,
                  " |\\-",
                  "_") %>%
  case_when(grepl('CD45', .) ~ "T_CD45.1",
            grepl('GC', .) ~ 'B_GC',
            .default = .) %>% 
  factor()

monocle.obj@colData$ct_cluster <- paste(
  monocle.obj@colData$Cluster, monocle.obj@colData$celltype, sep = "_"
  )

### heatmaps of ref. based ct anno. --------------------------------------

# setup list of singleR results 

results <- list()
for(variable in ls(pattern = "^orig")){
    results[[variable]] <- get(variable)
}

# compute the percentages of each cell type prediction in each cluster
freq_list <- list()
for (ref in names(results)) {
    freq_list[[ref]] <-
        as.data.frame(table(
            cluster = colData(monocle.obj)$Cluster,
            label = results[[ref]]$labels
        ) / (as.numeric(rep(
            table(colData(monocle.obj)$Cluster),
            length(table(results[[ref]]$labels))
        ))) * 100) %>%
        filter(Freq > 5) %>%  # filter for results with more than 5% prevelance
        mutate(reference = gsub("^.*?s_", "", ref))
    
}

freq_list <- bind_rows(freq_list)
colnames(freq_list)[1] <- "Cluster"
freq_list <- merge(freq_list, ct_majority, by = intersect(names(freq_list),names(ct_majority))) %>% dplyr::select(-n)

freq_list <- freq_list %>% 
  mutate(celltype_raw = str_replace_all(celltype_raw, " |\\-", "_") %>% 
           case_when(grepl('CD45', .) ~ "T_CD45.1",
                     grepl('GC', .) ~ 'B_GC',
                     .default = .)) %>% 
  arrange(celltype_raw, Cluster) %>% 
  unite("ct_cluster", c(Cluster, celltype_raw)) %>% 
  mutate(ct_cluster = fct_inorder(ct_cluster))

write_rds(freq_list, file.path(dataDir, "all_labels_ref_based_for_HM.rds"))

# plot a heatmap
freq_list %>% 
  ggplot(aes(x = label, y = ct_cluster, fill = Freq)) +
  geom_tile(colour = "white", show.legend = FALSE) +
  scale_fill_gradient(low = "ivory2", high = "red") +
  theme_my() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      size = 8,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  facet_grid(~reference,
    scales = "free",
    space = "free"
  )

ggsave(file.path(plotsDir, "all_labels_HM.pdf"), height = 13, width = 30)


# Save Object -------------------------------------------------------------
monocle.obj@colData$ct_cluster <-
  factor(monocle.obj@colData$ct_cluster,
         levels = levels(freq_list$ct_cluster))

write_rds(monocle.obj, file.path(dataDir, "3_annotated_monocle.cds"))


# Subclustering -----------------------------------------------------------

monocle.obj <- read_rds(file.path(dataDir, "3_annotated_monocle.cds"))

t_subst <- choose_cells(monocle.obj)
t_subst <- cluster_cells(t_subst, resolution = 0.0005)

write_rds(t_subst, file.path(dataDir, '3_t_subset.cds'))

b_subst <- choose_cells(monocle.obj)
b_subst <- cluster_cells(b_subst, resolution = 0.00009)

write_rds(b_subst, file.path(dataDir, '3_b_subset.cds'))


# Ref Data for T subsets --------------------------------------------------

th_express <- fread(file.path(tablesDir, 'data-2.csv'))[,-1]
th_express <- as.data.frame(th_express)
rownames(th_express) <- th_express$gene_symbol
th_express <- as.matrix(th_express[,-1])

th_express <- as(th_express, 'sparseMatrix')


# Integrating New Clusters back into data ---------------------------------

monocle.obj <- read_rds(file.path(dataDir, "3_annotated_monocle.cds"))
t_subset <- read_rds(file.path(dataDir, '3_t_subset.cds'))

increase_fac <- function(factor, increase){
  names <- names(factor)
  values <- unname(factor) %>% as.numeric
  values <- values + increase
  factor(setNames(values, names))
}

t_subset@clusters@listData$UMAP$clusters <- increase_fac(clusters(t_subset), 69)

comb_factors <- function(factor1, factor2){
  names <- c(names(factor1), names(factor2))
  fact <- forcats::fct_c(factor1, factor2) %>% droplevels
  
  names(fact) <- names
  fact
}
monocle.obj@clusters$UMAP$clusters <- 
  comb_factors(monocle.obj@clusters$UMAP$clusters[!names(monocle.obj@clusters$UMAP$clusters) %in% names(clusters(t_subset))],
             clusters(t_subset))

monocle.obj@clusters$UMAP$clusters <- 
  monocle.obj@clusters$UMAP$clusters[match(colnames(monocle.obj), names(monocle.obj@clusters$UMAP$clusters))]

tab <- monocle.obj %>% clusters %>% unname %>% table
r <- rank(-tab, ties.method = 'first')
levels(monocle.obj@clusters$UMAP$clusters) <- unname(r) 

b_subset <- read_rds(file.path(dataDir, '3_b_subset.cds'))

b_subset@clusters@listData$UMAP$clusters <- increase_fac(clusters(b_subset), 69)

monocle.obj@clusters$UMAP$clusters <- 
  comb_factors(monocle.obj@clusters$UMAP$clusters[!names(monocle.obj@clusters$UMAP$clusters) %in% names(clusters(b_subset))],
               clusters(b_subset))

monocle.obj@clusters$UMAP$clusters <- 
  monocle.obj@clusters$UMAP$clusters[match(colnames(monocle.obj), names(monocle.obj@clusters$UMAP$clusters))]

tab <- monocle.obj %>% clusters %>% unname %>% table
r <- rank(-tab, ties.method = 'first')
levels(monocle.obj@clusters$UMAP$clusters) <- unname(r) 
