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
      linewidth = 1
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

reference.obj@metadata <- list(type = 'sc')

old.mon <- readRDS(file.path(oldDir, "scRNAseq_3_monocle_more_hash_cutoff.cds"))
old_meta <- as.data.frame(colData(old.mon))
old_counts <- assay(old.mon)
old.mon <- CreateSeuratObject(counts = old_counts,
                              project = "ol",
                              assay = "counts",
                              meta.data = old_meta)
old.mon <-  as.SingleCellExperiment(old.mon)

colnames(colData(old.mon))[colnames(colData(old.mon)) == "celltype"] <- "label.fine"

old.mon@metadata <- list(type = 'sc')
# SingleR annotation ------------------------------------------------------

## Reference data ----------------------------------------------------------

### making a list of reference data sets  ----------------------------------------------

ref_data <- list(
  #Human Data Sets: 
  ##HPCA generally / Keratinocytes
  # hpca = HumanPrimaryCellAtlasData(),
  # 
  ##DB ImmuneCells, comprehensive CD4+ subsets; only one B cell subset, no dendritic cells
  # dice = DatabaseImmuneCellExpressionData(),
  # 
  # ##Monaco Immunc cells
  # monaco = MonacoImmuneData(),
  
  #MouseData:
  ##MouseRNAseqData
  # mrsd = MouseRNAseqData(),
  
  ##ImmGen 
  immgen = ImmGenData(),
  
  ##reference Dataset
  struc = reference.obj,
  
  ##old dataset
  old = old.mon
)

rm(old.mon)
rm(old_counts)
rm(old_meta)
rm(reference.obj)
gc()

convert_human_to_mouse <- function(gene_list) {
  mouse_human_genes = fread("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt")
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
# rownames(ref_data$dice) <- convert_human_to_mouse(rownames(ref_data$dice))
# rownames(ref_data$monaco) <- str_to_title(rownames(ref_data$monaco))

sce <- CreateSeuratObject(counts = assay(monocle.obj),
                          project = "integrated_scRNAseq",
                          assay = "integrated",
                          meta.data = as.data.frame(colData(monocle.obj)))


sce <- as.SingleCellExperiment(sce)



### run singleR with additional saving of original results data -------------------------------
# parallel computation
res <- list()
doParallel::registerDoParallel(cores = 7)
sc_res <- foreach(ref = names(ref_data)) %dopar% {
  message(ref)
  for(labelx in c("label.main", "label.fine")){
    message(paste(".", labelx))
    
    if(!labelx %in% colnames(colData(ref_data[[ref]]))) next
    
    message(paste(".", "running SingleR"))
    
    use_method <- fifelse(is.null(ref_data[[ref]]@metadata$type),
                          'classic',
                          'wilcox')
    
    use_aggr <- fifelse(use_method == 'classic',
                        FALSE,
                        TRUE)
    
    results <- SingleR(
      test = sce,
      ref = ref_data[[ref]],
      labels = colData(ref_data[[ref]])[, labelx],   
      aggr.ref = use_aggr,
      de.method = use_method
    )
    
    res[[labelx]] <- data.table(
      as_tibble(results, rownames = "cell"),
      ref = ref,
      labels = labelx
    )
    
    colnames(res[[labelx]]) <- gsub("\\.", "_", colnames(res[[labelx]]))
    
    res[[labelx]] <- res[[labelx]][, 
                                   .(cell, 
                                     labels, 
                                     tuning_scores_first, 
                                     tuning_scores_second)]
    
    cols <- res[[labelx]][,
                          .SD,
                          .SDcols = is.numeric] %>%
      colnames()
    
    res[[labelx]][,
                  (cols) := round(.SD, 2),
                  .SDcols = cols]
    
  }
  res
}

sc_res <- setNames(sc_res, names(ref_data))

sc_res <- rbindlist(lapply(sc_res, rbindlist,
                           idcol = 'label'),
                    idcol = 'ref')

# Save Raw assignment data ------------------------------------------------

fwrite(sc_res, file.path(tablesDir, 'cell_types_all_ref.csv'))


### assign labels to colData ----------------------------------------------------
sc_res <- fread(file.path(tablesDir, 'cell_types_all_ref.csv'))

sc_res <- sc_res[, dcast(.SD,
                         cell ~ ref + label,
                         value.var = 'labels')]


colData(monocle.obj) <- merge(colData(monocle.obj), 
                              sc_res, 
                              by.x = 0, 
                              by.y = 'cell',
                              sort = F)

rownames(colData(monocle.obj)) <- colData(monocle.obj)$Row.names
colData(monocle.obj)$Row.names <- NULL


monocle.obj@colData$celltype_raw <- monocle.obj@colData$old_label.fine

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

res <- fread(file.path(tablesDir, 'cell_types_all_ref.csv'))

res[, cluster := rep(clusters(monocle.obj), 4)]
res <- res[,
           keyby = .(ref, label, cluster, labels),
           .N][,
               by = .(ref, label, cluster),
               freq := N / sum(N) * 100][freq > 5,-c('N')]

res[,
    label := do.call(paste, c(.SD, sep = '_')),
    .SDcols = 1:2]


res <- merge(res, ct_majority, by.x = 'cluster', by.y = 'Cluster')


res <- res %>%
  mutate(celltype_raw = str_replace_all(celltype_raw, " |\\-", "_") %>%
           case_when(grepl('CD45', .) ~ "T_CD45.1",
                     grepl('GC', .) ~ 'B_GC',
                     .default = .)) %>%
  arrange(celltype_raw, cluster) %>%
  unite("ct_cluster", c(cluster, celltype_raw)) %>%
  mutate(ct_cluster = fct_inorder(ct_cluster),
         n = NULL,
         ref = NULL)

write_rds(res, file.path(dataDir, "all_labels_ref_based_for_HM.rds"))

# plot a heatmap
res %>% 
  ggplot(aes(x = labels, y = ct_cluster, fill = freq)) +
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
  facet_grid(~ label,
             scales = "free",
             space = "free"
  )

ggsave(file.path(plotsDir, "all_labels_HM.pdf"), height = 13, width = 20)


# Save Object -------------------------------------------------------------
monocle.obj@colData$ct_cluster <-
  factor(monocle.obj@colData$ct_cluster,
         levels = levels(res$ct_cluster))

write_rds(monocle.obj, file.path('/vscratch/scRNAseq/data', "pre_annotation_monocle.cds"))


# Subclustering -----------------------------------------------------------

monocle.obj <- read_rds(file.path('/vscratch/scRNAseq/data', "pre_annotation_monocle.cds"))

t_subst <- choose_cells(monocle.obj)
t_subst <- cluster_cells(t_subst, resolution = 0.0005)

write_rds(t_subst, file.path('/vscratch/scRNAseq/data', '3_t_subset.cds'))

b_subst <- choose_cells(monocle.obj)
b_subst <- cluster_cells(b_subst, resolution = 0.00009)

write_rds(b_subst, file.path('/vscratch/scRNAseq/data', '3_b_subset.cds'))

# Integrating New Clusters back into data ---------------------------------

#Most of the following is done to order the factor levels in the final monocle
#object by factor size 
monocle.obj <- read_rds(file.path('/vscratch/scRNAseq/data', "pre_annotation_monocle.cds"))
t_subset <- read_rds(file.path('/vscratch/scRNAseq/data', '3_t_subset.cds'))

#First the factor levels of the subsets are increased
increase_fac <- function(factor, increase){
  names <- names(factor)
  values <- unname(factor) %>% as.numeric
  values <- values + increase
  factor(setNames(values, names))
}

by_cl <- length(levels(monocle.obj@clusters@listData$UMAP$clusters))
t_subset@clusters@listData$UMAP$clusters <- increase_fac(clusters(t_subset), by_cl)

#Then the factors get combined in the monocle object
comb_factors <- function(factor1, factor2){
  names <- c(names(factor1), names(factor2))
  fact <- forcats::fct_c(factor1, factor2) %>% droplevels
  
  names(fact) <- names
  fact
}
monocle.obj@clusters$UMAP$clusters <- 
  comb_factors(monocle.obj@clusters$UMAP$clusters[!names(monocle.obj@clusters$UMAP$clusters) %in% names(clusters(t_subset))],
               clusters(t_subset))

#Make sure the values are in the right order
monocle.obj@clusters$UMAP$clusters <- 
  monocle.obj@clusters$UMAP$clusters[match(colnames(monocle.obj), names(monocle.obj@clusters$UMAP$clusters))]

#Lastly, the factor levels are adjusted and ordered
tab <- fct_count(clusters(monocle.obj))
r <- rank(-tab$n, ties.method = 'first')

monocle.obj@clusters$UMAP$clusters <- clusters(monocle.obj) %>% 
  lvls_revalue(r %>% as.character) %>% 
  fct_inseq()

#Same for B cell subset
b_subset <- read_rds(file.path('/vscratch/scRNAseq/data', '3_b_subset.cds'))

by_cl <- length(levels(monocle.obj@clusters@listData$UMAP$clusters))
b_subset@clusters@listData$UMAP$clusters <- increase_fac(clusters(b_subset), by_cl)

monocle.obj@clusters$UMAP$clusters <- 
  comb_factors(monocle.obj@clusters$UMAP$clusters[!names(monocle.obj@clusters$UMAP$clusters) %in% names(clusters(b_subset))],
               clusters(b_subset))

monocle.obj@clusters$UMAP$clusters <- 
  monocle.obj@clusters$UMAP$clusters[match(colnames(monocle.obj), names(monocle.obj@clusters$UMAP$clusters))]

tab <- fct_count(clusters(monocle.obj))
r <- rank(-tab$n, ties.method = 'first')

monocle.obj@clusters$UMAP$clusters <- clusters(monocle.obj) %>% 
  lvls_revalue(r %>% as.character) %>% 
  fct_inseq()

#Reassiging right cluster levels in colData
monocle.obj$Cluster <- clusters(monocle.obj) %>% unname()

#also for ct_cluster
monocle.obj$ct_cluster <-
  paste(monocle.obj$Cluster, monocle.obj$celltype, sep = "_")

#sorting by character AND numeric in same string is tedious
fct_order <-
  data.table(stringi::stri_split_fixed(unique(monocle.obj$ct_cluster),
                                       "_",
                                       2,
                                       simplify = T))[
                                         order(V2, as.numeric(V1))
                                       ][,
                                         do.call(paste,
                                                 c(.SD,
                                                   sep = '_'))]

monocle.obj$ct_cluster <- factor(monocle.obj$ct_cluster, levels = fct_order)


# Ref Data for T subsets --------------------------------------------------

### Scraped data was removed now, because raw data provided similar(or even 
### better) results. Also, there is at least 1 replicate of samples.

### TH scraped data
# th_express <- fread(file.path(tablesDir, 'th_ref.csv'))
# th_express <- as.data.frame(th_express)
# rownames(th_express) <- th_express$V1
# th_express <- as.matrix(th_express[,-1])
# 
# th_express <- as(th_express, 'sparseMatrix')
# met <- data.frame(label.fine = colnames(th_express), row.names = colnames(th_express))
# 
# th_scrape <- CreateSeuratObject(counts = th_express, meta.data = met)
# th_scrape <- as.SingleCellExperiment(th_scrape)

### TH Raw data
th_counts <- fread('https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/582/E-MTAB-2582/Files/Teichmann-ThExpress_rawCounts.txt')

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host = 'http://oct2014.archive.ensembl.org')
genes <- getBM(filters = 'ensembl_gene_id', 
               attributes= c("ensembl_gene_id",
                             'external_gene_name'),
               values = th_counts$ENSEMBL_ID,
               mart = ensembl)

th_counts <- th_counts[ENSEMBL_ID %in% genes$ensembl_gene_id, 
                       gene_id := genes$external_gene_name
][, 
  na.omit(.SD)
][,
  lapply(.SD, mean),
  by = gene_id,
  .SDcols = -1
][gene_id != '']

th_counts <- as.data.frame(th_counts)
rownames(th_counts) <- th_counts$gene_id
th_counts$gene_id <- NULL

th_counts <- as(as.matrix(th_counts), 'sparseMatrix')

th_meta <- fread('https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/582/E-MTAB-2582/Files/E-MTAB-2582.sdrf.txt')
th_meta <- th_meta[, .(`Extract Name`, label.fine = str_extract(`Extract Name`, '^[:alnum:]+'))] %>% unique()
th_meta <- data.frame(label.fine = th_meta$label.fine, row.names = colnames(th_counts))

th_raw <- CreateSeuratObject(counts = th_counts, meta.data = th_meta)
th_raw <- as.SingleCellExperiment(th_raw)
th_raw <- scuttle::logNormCounts(th_raw)

### Rest of T data

t_sort_info <- readxl::read_xlsx(file.path(tablesDir, 'Sort format for Single cells HDM.xlsx'))
t_sort_info <- setDT(t_sort_info[-1, -2])
t_sort_info <- t_sort_info[, melt(.SD,
                                  id.vars = 1,
                                  value.name = 'label.fine')][,
                                                              sample := do.call(paste, c(.SD, sep = '')),
                                                              .SDcols = 1:2][, -(1:2)]
t_sort_info[is.na(label.fine), label.fine := 'unknown']
t_sort_info <- as.data.frame(t_sort_info)
rownames(t_sort_info) <- t_sort_info$sample
t_sort_info$sample <- NULL

t_counts <- fread(file.path(tablesDir, 'GSE131935_SS2_17_449_rpkms.tab'))
t_counts <- t_counts[, lapply(.SD, mean), by = gene, .SDcols = -1]
t_counts <- as.data.frame(t_counts)
rownames(t_counts) <- t_counts$gene
t_counts$gene <- NULL
t_counts <- as(as.matrix(t_counts), 'sparseMatrix')


th_rest <- CreateSeuratObject(counts = t_counts, meta.data = t_sort_info)
th_rest <- as.SingleCellExperiment(th_rest)
th_rest <- scuttle::logNormCounts(th_rest)


ref_data <- list(
  th_raw = th_raw,
  th_rest = th_rest
)

sce <- CreateSeuratObject(counts = assay(monocle.obj),
                          project = "integrated_scRNAseq",
                          assay = "integrated",
                          meta.data = as.data.frame(colData(monocle.obj)))


sce <- as.SingleCellExperiment(sce)


res <- list()
doParallel::registerDoParallel(cores = 7)
sc_res <- foreach(ref = names(ref_data)) %dopar% {
  for(labelx in c("label.main", "label.fine")){
    
    if(!labelx %in% colnames(colData(ref_data[[ref]]))) next
    
    use_method <- fifelse(is.null(ref_data[[ref]]@metadata$type),
                          'classic',
                          'wilcox')
    
    use_aggr <- fifelse(use_method == 'classic',
                        FALSE,
                        TRUE)
    
    results <- SingleR(
      test = sce,
      ref = ref_data[[ref]],
      labels = colData(ref_data[[ref]])[, labelx],   
      aggr.ref = use_aggr,
      de.method = use_method
    )
    
    res[[labelx]] <- data.table(
      as_tibble(results, rownames = "cell"),
      ref = ref,
      labels = labelx
    )
    
    colnames(res[[labelx]]) <- gsub("\\.", "_", colnames(res[[labelx]]))
    
    res[[labelx]] <- res[[labelx]][, 
                                   .(cell, 
                                     labels, 
                                     tuning_scores_first, 
                                     tuning_scores_second)]
    
    cols <- res[[labelx]][,
                          .SD,
                          .SDcols = is.numeric] %>%
      colnames()
    
    res[[labelx]][,
                  (cols) := round(.SD, 2),
                  .SDcols = cols]
    
  }
  res
}

sc_res <- setNames(sc_res, names(ref_data))

sc_res <- rbindlist(lapply(sc_res, rbindlist,
                           idcol = 'label'),
                    idcol = 'ref')

# Save Raw assignment data ------------------------------------------------

fwrite(sc_res, file.path(tablesDir, 'cell_types_t_sub_ref.csv'))

### assign labels to colData ----------------------------------------------------
sc_res <- fread(file.path(tablesDir, 'cell_types_t_sub_ref.csv'))

sc_res <- sc_res[, dcast(.SD,
                         cell ~ ref + label,
                         value.var = 'labels')]


colData(monocle.obj) <- merge(colData(monocle.obj), 
                              sc_res, 
                              by.x = 0, 
                              by.y = 'cell',
                              sort = F)

rownames(colData(monocle.obj)) <- colData(monocle.obj)$Row.names
colData(monocle.obj)$Row.names <- NULL

res <- fread(file.path(tablesDir, 'cell_types_t_sub_ref.csv'))
res2 <- fread(file.path(tablesDir, 'cell_types_all_ref.csv'))

res <- rbind(res,res2)
res[, cluster := rep(clusters(monocle.obj), 6)]
res <- res[,
           keyby = .(ref, label, cluster, labels),
           .N][,
               by = .(ref, label, cluster),
               freq := N / sum(N) * 100][freq > 5,-c('N')]

res[,
    label := do.call(paste, c(.SD, sep = '_')),
    .SDcols = 1:2]


#Maybe newer version of ct_majority with majority of all ref data bases
cols <- colnames(monocle.obj@colData) %>% grep(pattern = 'label', value = T)
ct_majority <- as.data.table(monocle.obj@colData)[, c('Cluster', ..cols)]
ct_majority <-
  ct_majority[, melt(.SD, 'Cluster', var = 'ref', val = 'label')]
ct_majority <- ct_majority[ct_majority[, .I[which.max(.N)], by= .(Cluster, ref)]$V1]
ct_majority <- ct_majority[, dcast(.SD, Cluster ~ ref)]
#combination of old data and one t subset ref data
ct_majority[, th_raw_label.fine := fcase(th_raw_label.fine == 'Naive', 'T_naive',
                                         th_raw_label.fine == 'iTreg', 'Treg_induced',
                                         rep(TRUE, .N), th_raw_label.fine)]
ct_majority[, celltype_raw := fcase(old_label.fine %in% c('Tgd', 'T cells'), th_raw_label.fine,
                                    !old_label.fine %in% c('Tgd', 'T cells'), old_label.fine)]
ct_majority <- ct_majority[, .(Cluster, celltype_raw)]



res <- merge(res, ct_majority, by.x = 'cluster', by.y = 'Cluster')

res[, celltype_raw := str_replace_all(celltype_raw, " |\\-", "_")]
res[, celltype_raw := fcase(grepl('CD45', celltype_raw), 'T_CD45.1',
                            grepl('GC', celltype_raw), 'B_GC',
                            rep(TRUE, .N), celltype_raw)]
res[,
    ct_cluster := do.call(paste, c(.SD, sep = '_')), 
    .SDcols = c('cluster', 'celltype_raw')]

fc_levels <- unique(res[order(celltype_raw, cluster), ct_cluster])
res[,
    ct_cluster := factor(ct_cluster,
                         levels = fc_levels)]

write_rds(res, file.path(dataDir, "all_labels_all_clusters_ref_based_for_HM.rds"))

# plot a heatmap
res %>% 
  ggplot(aes(x = labels, y = ct_cluster, fill = freq)) +
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
  facet_grid(~ label,
             scales = "free",
             space = "free"
  )

ggsave(file.path(plotsDir, "all_clusters_all_ct_labels_HM.pdf"), height = 15, width = 28)


# Assign to monocle object ------------------------------------------------

monocle.obj@colData$celltype <- res$celltype_raw[match(monocle.obj$Cluster, res$cluster)]

monocle.obj@colData$ct_cluster <- res$ct_cluster[match(monocle.obj$Cluster, res$cluster)]

write_rds(monocle.obj, file.path(dataDir, '3_annotated_monocle.cds'))

