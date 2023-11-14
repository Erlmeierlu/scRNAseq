library(renv)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(matrixStats)
library(hexbin)
library(patchwork)
library(monocle3)

#set seed
set.seed(42)
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

#Functions
make_col_stats <- function(counts.df) {
  x <- as.matrix(counts.df)
  
  hashtags <-
    data.frame(sum = colSums(x[grepl("-[0-9]+$", rownames(x)),]),
               max = colMaxs(x[grepl("-[0-9]+$", rownames(x)),]),
               ratio = colMaxs(x[grepl("-[0-9]+$", rownames(x)),]) / colSums(x[grepl("-[0-9]+$", rownames(x)),]))
  
  cd45.1 <- 
    data.frame(sum = colSums(x[grepl("CD45.*$", rownames(x)),]),
               cd45.1 = x[grepl("CD45.1$", rownames(x)),],
               ratio = x[grepl("CD45.1$", rownames(x)),] / colSums(x[grepl("CD45.*$", rownames(x)),]))
  
  cd4 <-
    data.frame(sum = colSums(x[grepl("CD[1-9]$", rownames(x)),]),
               cd4 = x[grepl("CD4$", rownames(x)),],
               ratio = x[grepl("CD4$", rownames(x)),] / colSums(x[grepl("CD[1-9]$", rownames(x)),]))
  
  return(list(hashtags = hashtags, cd45x = cd45.1, cd4_cd8 = cd4))
}

make_row_stats <- function(counts.df) {
  x <- as.matrix(counts.df)
  
  hashtags <-
    data.frame(sum = rowSums(x[grepl("-[0-9]+$", rownames(x)),]),
               max = rowMaxs(x[grepl("-[0-9]+$", rownames(x)),]),
               ratio = rowMaxs(x[grepl("-[0-9]+$", rownames(x)),]) / rowSums(x[grepl("-[0-9]+$", rownames(x)),]))
  
  cd45x <- 
    data.frame(sum = rowSums(x[grepl("CD45.*$", rownames(x)),]),
               max = rowMaxs(x[grepl("CD45.*$", rownames(x)),]),
               ratio = rowMaxs(x[grepl("CD45.*$", rownames(x)),]) / rowSums(x[grepl("CD45.*$", rownames(x)),]))
  
  cd4_cd8 <-
    data.frame(sum = rowSums(x[grepl("CD[1-9]$", rownames(x)),]),
               max = rowMaxs(x[grepl("CD[1-9]$", rownames(x)),]),
               ratio = rowMaxs(x[grepl("CD[1-9]$", rownames(x)),]) / rowSums(x[grepl("CD[1-9]$", rownames(x)),]))
  
  return(list(hashtags = hashtags, cd45x = cd45x, cd4_cd8 = cd4_cd8))
}

make_CITE_stats <- function(counts.df) {
  x <- as.matrix(counts.df)
  
  ratios <- 
    data.frame(cd45.1_.2_ratio = (x[grepl("CD45.1$", rownames(x)),]+1)/(x[grepl("CD45.2$", rownames(x)),]+1),
               cd4_8_ratio = (x[grepl("CD4$", rownames(x)),]+1)/(x[grepl("CD8$", rownames(x)),]+1))
  
  return(ratios)
}

DetermineDimensionality <- function(object){
  pct <- object[["pca"]]@stdev / sum(object[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  pcs <- min(co1, co2)
  return(pcs)
}

# Load Data ---------------------------------------------------------------

#setting up directories
baseDir <- getwd()
rawDir <- ("/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin/raw_data/data_from_BSF2/COUNT/")
plotsDir <- file.path(baseDir, "plots/")
tablesDir <- file.path(baseDir, "tables/")
dataDir <- file.path(baseDir, "data/")

#directories to loop over
fileDirs <- paste0(grep(list.files(rawDir), pattern = "AGG", invert = T, value = T), "/filtered_feature_bc_matrix/")


#Seurat Object
pbmc_data_list <- list()
pbmc_list <- list()
ab_list <- list()

for (file in fileDirs) {
  x <- str_extract(file, "[A-z]+[0-9]")
  
  pbmc_data_list[[x]] <- Read10X(file.path(rawDir, file))
  
  pbmc_list[[x]] <-
    CreateSeuratObject(
      counts = pbmc_data_list[[x]]$'Gene Expression',
      project = "scRNAseq",
      min.cells = 0,
      min.features = 200
    )
  
  ab_list[[x]] <- CreateSeuratObject(
    counts = pbmc_data_list[[x]]$'Antibody Capture',
    project = "scRNAseq",
    min.cells = 0,
    min.features = 0
  )
}

#creating list of counts
CountsAB <- lapply(ab_list, GetAssayData, slot = "counts")

# Quality Control ---------------------------------------------------------

ab_col_statlist <- lapply(CountsAB, make_col_stats)
ab_row_statlist <- lapply(CountsAB, make_row_stats)


ab_col_stats <- bind_rows(lapply(ab_col_statlist, bind_rows, .id = "abtype"), .id = "sample")
ab_row_stats <- bind_rows(lapply(ab_row_statlist, bind_rows, .id = "abtype"), .id = "sample")


CITE_statlist <- lapply(CountsAB, make_CITE_stats)
CITE_stats <- bind_rows(CITE_statlist, .id = "sample")

# *All ABs ---------------------------------------------------------------

ab_row_stats %>% group_by(abtype, sample) %>% summarise(sum = sum(sum)) %>% 
  ggplot() +
  geom_col(aes(abtype, sum)) + 
  facet_wrap(~ sample, scales = "free") + 
  theme_my() +
  labs(x = "Antibody Type", y = "Reads", title = "Antibody reads per sample")

#ggsave(file.path(plotsDir, "AB_reads_per_sample.pdf"))


# *Hashtag ABs ------------------------------------------------------------

ab_col_stats %>% filter(abtype == "hashtags") %>%  group_by(sum, sample) %>% dplyr::count() %>% 
  ggplot() +
  geom_col(aes(sum, n), col = "gray25", fill = "ivory2") +
  facet_wrap(~ sample, scales = "free") +
  theme_my() +
  labs(title = "Hashtag reads per sample", x = "Reads", y = "Number of Cells")

#ggsave(file.path(plotsDir, "Hashtag_reads_per_sample.pdf"))


ab_row_stats %>% filter(abtype == "hashtags") %>% 
  ggplot() +
  geom_col(aes(rownames(ab_row_stats %>% filter(abtype == "hashtags")) %>% str_remove("-[A-z]+-+[A-z]+[1-9]"), sum),
           col = "gray25",
           fill = "ivory2") +
  facet_wrap(~ sample, scales = "free_y") +
  theme_my() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1,
    size = 10
  )) +
  labs(title = "Hashtag subtype distribution", x = "Hashtag subtype", y = "Reads")

#ggsave(file.path(plotsDir, "Hashtag_subtype_distribution.pdf"))


ab_col_stats %>% filter(abtype == "hashtags", sum != 0) %>% 
  ggplot() +
  stat_binhex(aes(sum, ratio, fill = log(..count..)), bins = 30, col = "black") +
  facet_wrap( ~ sample, scales = "free_x") +
  scale_fill_gradient(low = "ivory2", high = "red") +
  theme_my() +
  labs(x = "Reads per Cell", y = "Ratio [MaxReads/SumReads]", title = "Hashtag ratio per sample")

#ggsave(file.path(plotsDir, "Hashtag_ratio_per_sample.pdf"))


ab_row_stats %>% filter(abtype == "hashtags", sum != 0) %>% 
  ggplot() +
  geom_point(aes(sum, ratio), shape = 21, alpha = 0.5, col = "black", size = 2, fill = "ivory2") +
  facet_wrap( ~ sample, scales = "free_x") +
  theme_my() +
  labs(title = "Distribution of reads per hashtag", x = "Reads per hashtag", y = "Ratio [MaxReads/SumReads]")

#ggsave(file.path(plotsDir, "Distribution_of_reads_per_hashtag.pdf"))


# *CD45 ABs ---------------------------------------------------------------

ab_col_stats %>%
  filter(abtype == "cd45x") %>%
  group_by(sample, sum == 0, ratio %in% c(0, 1), sum != 0 & !(ratio %in% c(0, 1))) %>%
  summarise(
    sum = sum(`sum == 0`),
    ratio = sum(`ratio %in% c(0, 1)`),
    both = sum(`sum != 0 & !(ratio %in% c(0, 1))`)
  ) %>%
  mutate(
    count = sum(c(sum, ratio, both), na.rm = T),
    ratio = replace_na(ratio, 0),
    Distinct_Tags = case_when(
      sum == 0 &
        ratio == 0 ~ "2",
      sum == 0 & both == 0 ~ "1",
      ratio == 0 & both == 0 ~ "0"
    )
  ) %>%
  ggplot() +
  geom_col(aes(Distinct_Tags, count)) +
  facet_wrap(~ sample, scales = "free_y") +
  labs(x = "Number of distinct tags", title = "Distinct CD45.1 / CD45.2 tags per sample") +
  theme_my()

#ggsave(file.path(plotsDir, "Distinct_CD45_tags_per_sample.pdf"))


ab_col_stats %>% filter(abtype == "cd45x", !(ratio %in% c(NA, NaN))) %>% 
  ggplot() + 
  geom_violin(aes(sample, ratio, fill = sample)) +
  theme_my() +
  labs(title = "Ratio of CD45.1 / CD45.x reads", y = "Ratio [CD45.1 Reads/ Sum Reads)]")

#ggsave(file.path(plotsDir, "Ratio_of_CD45_reads.pdf"))


ab_col_stats %>% filter(abtype == "cd45x") %>%  group_by(sum, sample) %>% dplyr::count() %>% 
  ggplot() +
  geom_col(aes(sum, n), col = "gray25", fill = "ivory2") +
  facet_wrap(~ sample, scales = "free") +
  theme_my() +
  labs(title = "CD45.x reads per sample", x = "Reads", y = "Number of Cells")

#ggsave(file.path(plotsDir, "CD45_reads_per_sample.pdf"))


ab_row_stats %>% filter(abtype == "cd45x") %>% 
  ggplot() +
  geom_col(aes(rownames(ab_row_stats %>% filter(abtype == "cd45x")) %>% str_extract("CD45.[1-2]"), sum),
           col = "gray25",
           fill = "ivory2") +
  facet_wrap(~ sample, scales = "free_y") +
  theme_my() +
  theme(axis.text.x = element_text(
    size = 10
  )) +
  labs(title = "CD45.1 / CD45.2 subtype distribution", x = "CITE subtype", y = "Reads")

#ggsave(file.path(plotsDir, "CD45_subtype_distribution.pdf"))


ab_col_stats %>% filter(abtype == "cd45x", sum != 0) %>% 
  ggplot() +
  stat_binhex(aes(sum, ratio, fill = log(..count..)), bins = 30, col = "black") +
  facet_wrap( ~ sample, scales = "free_x") +
  scale_fill_gradient(low = "ivory2", high = "red") +
  theme_my() +
  scale_x_log10() +
  labs(x = "Reads", title = "CD45.1 / CD45.x ratio per sample", y = "Ratio [CD45.1 Reads/ Sum Reads)]")

#ggsave(file.path(plotsDir, "CD45_ratio_per_sample.pdf"))


# *CD8/CD4 ABs -------------------------------------------------------------

ab_col_stats %>%
  filter(abtype == "cd4_cd8") %>%
  group_by(sample, sum == 0, ratio %in% c(0, 1), sum != 0 & !(ratio %in% c(0, 1))) %>%
  summarise(
    sum = sum(`sum == 0`),
    ratio = sum(`ratio %in% c(0, 1)`),
    both = sum(`sum != 0 & !(ratio %in% c(0, 1))`)
  ) %>%
  mutate(
    count = sum(c(sum, ratio, both), na.rm = T),
    ratio = replace_na(ratio, 0),
    Distinct_Tags = case_when(
      sum == 0 &
        ratio == 0 ~ "2",
      sum == 0 & both == 0 ~ "1",
      ratio == 0 & both == 0 ~ "0"
    )
  ) %>%
  ggplot() +
  geom_col(aes(Distinct_Tags, count)) +
  facet_wrap(~ sample, scales = "free_y") +
  labs(x = "Number of distinct tags", title = "Distinct CD4 / CD8 tags per sample") +
  theme_my()

#ggsave(file.path(plotsDir, "Distinct_CD4_tags_per_sample.pdf"))


ab_col_stats %>% filter(abtype == "cd4_cd8", !(ratio %in% c(NA, NaN))) %>% 
  ggplot() + 
  geom_violin(aes(sample, ratio, fill = sample)) +
  theme_my() +
  labs(title = "Ratio of CD4 / CDx reads", y = "Ratio [CD4 Reads / Sum Reads)]")

#ggsave(file.path(plotsDir, "Ratio_of_CD4_reads.pdf"))


ab_col_stats %>% filter(abtype == "cd4_cd8") %>%  group_by(sum, sample) %>% dplyr::count() %>% 
  ggplot() +
  geom_col(aes(sum, n), col = "gray25", fill = "ivory2") +
  facet_wrap(~ sample, scales = "free") +
  theme_my() +
  labs(title = "CD4 / CD8 reads per sample", x = "Reads", y = "Number of Cells")

#ggsave(file.path(plotsDir, "CD4_reads_per_sample.pdf"))


ab_row_stats %>% filter(abtype == "cd4_cd8") %>% 
  ggplot() +
  geom_col(aes(rownames(ab_row_stats %>% filter(abtype == "cd4_cd8")) %>% str_extract("CD[0-9]"), sum),
           col = "gray25",
           fill = "ivory2") +
  facet_wrap(~ sample, scales = "free_y") +
  theme_my() +
  theme(axis.text.x = element_text(
    size = 10
  )) +
  labs(title = "CD4 / CD8 subtype distribution", x = "CITE subtype", y = "Reads")

#ggsave(file.path(plotsDir, "CD4_subtype_distribution.pdf"))


ab_col_stats %>% filter(abtype == "cd4_cd8", sum != 0) %>% 
  ggplot() +
  stat_binhex(aes(sum, ratio,fill = log(..count..)), bins = 30, col = "black") +
  facet_wrap( ~ sample, scales = "free_x") +
  scale_fill_gradient(low = "ivory2", high = "red") +
  theme_my() +
  scale_x_log10() +
  labs(x = "Reads", title = "CD4 / CDx ratio per sample", y = "Ratio [CD4 Reads / Sum Reads)]")

#ggsave(file.path(plotsDir, "CD4_ratio_per_sample.pdf"))


# *CD45.1/.2 - CD4/8 combinations ------------------------------------------

CITE_stats %>% 
  ggplot() + 
  stat_bin_hex(aes(cd45.1_.2_ratio, cd4_8_ratio, fill = log(..count..)), bins = 30, col = "black") + 
  facet_wrap(~sample, scales = "free") + 
  scale_fill_gradient(low = "ivory2", high = "red") +
  scale_x_log10() +
  scale_y_log10() +
  theme_my() +
  labs(x = "Ratio [CD45.1:CD45.2]", y = "Ratio [CD4:CD8]", title = "Combination of CITE-tags")

#ggsave(file.path(plotsDir, "CITE_tag_combinations.pdf"))


# Assign metadata --------------------------------------------------------

AssignMetadata <- function(object, abcounts){
  sx <- object
  barcodes <- row.names(sx@meta.data)
  ab <- abcounts[,barcodes]
  ab.hash <- ab[!grepl('CD', row.names(ab)),]
  
  idx <- apply(t(t(as.matrix(ab.hash))/colSums(ab.hash)), 2, function(col) which(col > 0.6))
  
  stopifnot(all(colnames(ab.hash) == row.names(sx@meta.data)))
  
  hash <- lapply(idx, function(i) row.names(ab.hash)[i])
  hash[sapply(hash, length) == 0] <- NA
  
  sx@meta.data$hash <- unlist(hash)
  sx@meta.data$hash.sum <- colSums(ab.hash)
  sx@meta.data$hash.ratio <- (colMaxs(as.matrix(ab.hash))/colSums(ab.hash))
  sx@meta.data$CD45 <- ab['HTO-CD45.1',]/(ab['HTO-CD45.1',] + ab['HTO-CD45.2',])
  sx@meta.data$CD45.sum <- (ab['HTO-CD45.1',] + ab['HTO-CD45.2',])
  sx@meta.data$CD4 <- ab['HTO-CD4',]/(ab['HTO-CD4',] + ab['HTO-CD8',])
  sx@meta.data$CD4.sum <- (ab['HTO-CD4',] + ab['HTO-CD8',])
  object <- sx
  return(object)
}

pbmc_list <- mapply(AssignMetadata, pbmc_list, CountsAB)


for(sample in names(pbmc_list)){
  pbmc_list[[sample]]$percent.mt <- PercentageFeatureSet(pbmc_list[[sample]], pattern = "^mt-")
}


meta <- bind_rows(lapply(pbmc_list, function(sx) sx@meta.data), .id = 'sample')

ggplot(meta, aes(x=sample, y=percent.mt)) + geom_violin() + geom_hline(yintercept = 10)



# QC plots ----------------------------------------------------------------

violin_plot_list <-
  lapply(
    pbmc_list,
    VlnPlot,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  )

scatter_plot_list <- mapply(
  wrap_plots,
  lapply(
    pbmc_list,
    FeatureScatter,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  ),
  lapply(
    pbmc_list,
    FeatureScatter,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  ),
  SIMPLIFY = F
)

# CELL CYCLE
for(sx in names(pbmc_list)){
  # add sample name to metadata 
  pbmc_list[[sx]]$sample <- sx 
  # Cell cycle scoring (will be added to metadata)
  pbmc_list[[sx]] <- NormalizeData(pbmc_list[[sx]], verbose = FALSE)
  pbmc_list[[sx]] <- CellCycleScoring(pbmc_list[[sx]], 
                                   s.features = cc.genes$s.genes, 
                                   g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
}

#Filtering
pbmc_subset <-
  lapply(pbmc_list,
         subset,
         subset = nFeature_RNA > 200 &
           percent.mt < 10 & hash.sum > 10 & hash.ratio > 0.6)

monocle.obj.list <- list()
for(sx in names(pbmc_subset)) {
  seurat.obj <- pbmc_subset[[sx]]
  mat.use <- seurat.obj@assays$RNA@counts
  monocle.obj.list[[sx]] <-
    new_cell_data_set(expression_data = mat.use,
                      cell_metadata = seurat.obj@meta.data)
}

monocle.obj <- combine_cds(cds_list = monocle.obj.list, cell_names_unique = FALSE)

# Process dataset
monocle.obj <-
  preprocess_cds(monocle.obj, verbose = TRUE) %>%
  reduce_dimension(preprocess_method = "PCA", verbose = TRUE)


monocle.obj <- align_cds(monocle.obj,
                         alignment_group = "sample", 
                         residual_model_formula_str = "~G2M.Score + percent.mt + S.Score",
                         verbose = TRUE
)
monocle.obj <- reduce_dimension(monocle.obj,
                                reduction_method = "UMAP",
                                preprocess_method = "Aligned",
                                verbose = TRUE)

# Clustering

monocle.obj = cluster_cells(monocle.obj)

plot_cells(monocle.obj, color_cells_by = "sample")
plot_cells(monocle.obj, color_cells_by = "Phase")
plot_cells(monocle.obj, color_cells_by = "hash")

monocle.obj@colData$hash <- as.factor(monocle.obj@colData$hash)
monocle.obj@colData$sample <- as.factor(monocle.obj@colData$sample)
monocle.obj@colData$Phase <- as.factor(monocle.obj@colData$Phase)
monocle.obj@colData$Cluster <- unname(clusters(monocle.obj[,rownames(colData(monocle.obj))]))

rowData(monocle.obj)$gene_short_name <- row.names(rowData(monocle.obj))


# Store full dataset
saveRDS(monocle.obj, file = file.path(dataDir, "scRNAseq_1_monocle.cds"))

saveRDS(pbmc_subset, file = file.path(dataDir,"scRNAseq_1.rds"))



# Seurat workflow ---------------------------------------------------------

#Normalize Data
# pbmc_subset <- lapply(pbmc_subset, NormalizeData)
# 
# #highly variable features
# pbmc_subset <- lapply(pbmc_subset, FindVariableFeatures, selection.method = "vst", nfeatures = 2000)
# 
# top10 <- bind_rows(lapply(lapply(pbmc_subset,VariableFeatures), head, 10), .id = "sample")
# 
# variable_feature_plot_list <- mapply(
#   wrap_plots,
#   lapply(pbmc_subset, VariableFeaturePlot),
#   mapply(
#     LabelPoints,
#     plot = lapply(pbmc_subset, VariableFeaturePlot),
#     points = top10,
#     repel = T,
#     SIMPLIFY = F
#   ),
#   SIMPLIFY = F
# )
# 
# #Scale Data
# pbmc_subset <-
#   mapply(ScaleData, pbmc_subset, features = lapply(pbmc_subset, rownames))
# 
# #Linear Dim Reduction
# pbmc_subset <-
#   mapply(
#     RunPCA,
#     pbmc_subset,
#     features = lapply(pbmc_subset,
#                       VariableFeatures, ),
#     npcs = lapply(sapply(pbmc_subset, ncol) - 1, min, 50),
#     SIMPLIFY = F
#   )
# 
# pca_plot_list <- lapply(pbmc_subset, DimPlot, reduction = "pca")
# 
# dim_heatmap_list <-
#   lapply(
#     pbmc_subset,
#     DimHeatmap,
#     dims = 1:4,
#     cells = 500,
#     balanced = T,
#     fast = F,
#     combine = T,
#     ncol = 2
#   )
# 
# #Determine Dimensionality
# pbmc_subset <- lapply(pbmc_subset, JackStraw, num.replicate = 100)
# pbmc_subset <- lapply(pbmc_subset, ScoreJackStraw, dims = 1:20)
# 
# jack_straw_plot_list <- lapply(pbmc_subset, JackStrawPlot, dims = 1:20)
# elbow_plot_list <- lapply(pbmc_subset, ElbowPlot)
# 
# det_dim_seq_list <- mapply(seq, 1, lapply(pbmc_subset, DetermineDimensionality))
# 
# pbmc_subset <-
#   mapply(FindNeighbors, pbmc_subset, dims = det_dim_seq_list)
# 
# pbmc_subset <- lapply(pbmc_subset, FindClusters, resolution = 0.5)
# 
# #UMAP
# pbmc_subset <-
#   mapply(RunUMAP,
#          pbmc_subset,
#          dims = det_dim_seq_list,
#          n.neighbors = lapply(sapply(pbmc_subset, ncol) - 2, min, 30),
#          verbose = F)
# 
# umap_plot_list <- lapply(pbmc_subset, DimPlot, reduction = "umap")


