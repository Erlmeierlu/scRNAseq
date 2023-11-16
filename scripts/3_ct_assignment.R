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
vDir <- ("/vscratch/scRNAseq")
plotsDir <- file.path(vDir, "plots")
tablesDir <- file.path(vDir, "tables")
oldDir <- file.path(vDir, "data/old")
dataDir <-("data")
resDir <- ("results")

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
                                    project = "integrated_scRNAseq",
                                    assay = "integrated",
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

#Human Gene names as Mouse Gene names
# rownames(ref_data$hpca) <- str_to_title(rownames(ref_data$hpca))
rownames(ref_data$dice) <- str_to_title(rownames(ref_data$dice))
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
doParallel::registerDoParallel(cores = 7)
foreach(ref = names(ref_data)) %dopar% {
  print(ref)
  labelx <- "label.main"
  for(labelx in c("label.main", "label.fine")){
    print(paste(".", labelx))
    
    ref.file <- paste0("cell_types_", ref, "_", labelx, ".csv")
    
    ref.file.orig <- paste0("orig_cell_types_", ref, "_", labelx, ".rds")
    
    if(!labelx %in% colnames(colData(ref_data[[ref]]))) next
    
    print(paste(".", "running SingleR"))
    
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
  str_replace_all(ct_majority[match(monocle.obj@colData$Cluster, ct_majority$Cluster), ]$celltype_raw,
                  " |\\-",
                  "_")

monocle.obj@colData$ct_cluster <- paste(monocle.obj@colData$Cluster, monocle.obj@colData$celltype, sep = "_")

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
  arrange(celltype_raw, Cluster) %>% 
  mutate(celltype_raw = str_replace_all(celltype_raw, " |\\-", "_")) %>% 
  unite("ct_cluster", c(Cluster, celltype_raw)) %>% 
  mutate(ct_cluster = fct_inorder(ct_cluster))

saveRDS(freq_list, file.path(vDir, "data", "all_labels_ref_based_for_HM.rds"))

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

ggsave(file.path(plotsDir, "all_labels_HM.pdf"), height = 13, width = 25)


# Save Object -------------------------------------------------------------
monocle.obj@colData$ct_cluster <-
  factor(monocle.obj@colData$ct_cluster,
         levels = levels(freq_list$ct_cluster))

write_rds(monocle.obj, file.path(dataDir, "3_annotated_monocle.cds"))
