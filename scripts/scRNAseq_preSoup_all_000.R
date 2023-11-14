library(dplyr)
library(stringr)
library(Seurat)
library(monocle3)
library(gt)
library(Matrix)
library(scds)
library(DoubletFinder)
library(scCustomize)


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
vDir <- ("/vscratch/scRNAseq")
rawDir <- ("/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin/raw_data/scRNA_from_BSF/COUNT")
plotsDir <- file.path(vDir, "plots")
tablesDir <- file.path(vDir, "tables")
dataDir <- file.path(vDir, "data")
resDir <- file.path(baseDir, "results")

#directories to loop over
fileDirs <- list()
fileDirs[["filtered"]] <- paste0(grep(list.files(rawDir), pattern = "AGG", invert = T, value = T), "/filtered_feature_bc_matrix/")
fileDirs[["raw"]] <- paste0(grep(list.files(rawDir), pattern = "AGG", invert = T, value = T), "/raw_feature_bc_matrix/")

fileDirs <- lapply(fileDirs, function(x) {
  x[!grepl("fLN_40B3", x)]
})

#Seurat Object
pbmc_data_list <- list()
pbmc_list <- list()
# ab_list <- list()

for (file in fileDirs[["filtered"]]) {
    
    x <- file %>% str_remove("_trans.*")
    
    experiment <- case_when(!grepl("[0-9]+B[0-9]", x) ~ "HDAC1",
                            grepl("40", x) ~ "HDAC1",
                            TRUE ~ "HDAC2")
    
    pbmc_data_list[[experiment]][[x]] <- Read10X(file.path(rawDir, file))
    
    pbmc_list[[experiment]][[x]] <-
        CreateSeuratObject(
            counts = pbmc_data_list[[experiment]][[x]]$'Gene Expression',
            project = "scRNAseq",
            min.cells = 0,
            min.features = 0
        )
    
    # ab_list[[experiment]][[x]] <- CreateSeuratObject(
    #   counts = pbmc_data_list[[experiment]][[x]]$'Antibody Capture',
    #   project = "scRNAseq",
    #   min.cells = 0,
    #   min.features = 0
    # )
}

pbmc_data_list_raw <- list()
pbmc_list_raw <- list()
for (file in fileDirs[["raw"]]) {
  
  x <- file %>% str_remove("_trans.*")
  
  experiment <- case_when(!grepl("[0-9]+B[0-9]", x) ~ "HDAC1",
                          grepl("40", x) ~ "HDAC1",
                          TRUE ~ "HDAC2")
  
  pbmc_data_list_raw[[experiment]][[x]] <- Read10X(file.path(rawDir, file))
  
  pbmc_list_raw[[experiment]][[x]] <-
    CreateSeuratObject(
      counts = pbmc_data_list_raw[[experiment]][[x]]$'Gene Expression',
      project = "scRNAseq",
      min.cells = 0,
      min.features = 0
    )
  
  # ab_list[[experiment]][[x]] <- CreateSeuratObject(
  #   counts = pbmc_data_list[[experiment]][[x]]$'Antibody Capture',
  #   project = "scRNAseq",
  #   min.cells = 0,
  #   min.features = 0
  # )
}


# CountsAB <- list()
# #creating list of counts
# for (experiment in names(ab_list)){
# CountsAB[[experiment]] <- lapply(ab_list[[experiment]], GetAssayData, slot = "counts")
# }

# saveRDS(CountsAB, file.path(dataDir, "ab_counts.rds"))
CountsAB <- readRDS(file.path(dataDir, "ab_counts.rds"))

  #removal of bad sample
CountsAB <- lapply(CountsAB, function(x) {
    x[!grepl("fLN_40B3", names(x))]
  })

AssignMetadata <- function(object, abcounts){
    sx <- object
    barcodes <- row.names(sx@meta.data)
    ab <- abcounts[,barcodes]
    ab.hash <- ab[grepl("-[0-9]+$", row.names(ab)),]
    
    ab.aggr <-
        t(matrix(
            data = c(colSums(ab.hash[1:3,]), colSums(ab.hash[4:6,]), colSums(ab.hash[7:9,])) ,
            ncol = 3,
            dimnames = list(
                colnames(ab.hash),
                ifelse(
                    !grepl("\\-[0-9]\\-", rownames(ab.hash)),
                    case_when(
                        grepl("30[1-3]", rownames(ab.hash)) ~ paste0(rownames(ab.hash), "-1"),
                        grepl("30[4-6]", rownames(ab.hash)) ~ paste0(rownames(ab.hash), "-2"),
                        grepl("30[7-9]", rownames(ab.hash)) ~ paste0(rownames(ab.hash), "-3")
                    ) %>% str_remove("\\-30[0-9]"),
                    rownames(ab.hash) %>% str_remove("\\-30[0-9]")
                ) %>% unique()
            )
        ))
    ab.aggr <- as(ab.aggr,"sparseMatrix")
    
    stopifnot(colSums(ab.aggr) == colSums(ab.hash))
    
    ab.cd4 <- ab[grepl("CD[1-9]$", row.names(ab)),]
    ab.cd45 <- ab[grepl("CD45.*$", row.names(ab)),]
    
    if(any(grepl("GD", rownames(ab)))){
        ab.gd <- ab[grepl("GD", rownames(ab)),]
        gd.counts <- ab.gd
        ab.gdcd4 <- ab[grepl("GD|CD4$", rownames(ab)),]
    }
    
    idx.aggr <- apply(t(t(as.matrix(ab.aggr))/colSums(ab.aggr)), 2, function(col) which(col > 0.6))
    idx.hash <- apply(t(t(as.matrix(ab.hash))/colSums(ab.hash)), 2, function(col) which(col > 0.6))
    # idx.hash.55 <- apply(t(t(as.matrix(ab.hash))/colSums(ab.hash)), 2, function(col) which(col > 0.55))
    # idx.hash.07 <- apply(t(t(as.matrix(ab.hash))/colSums(ab.hash)), 2, function(col) which(col > 0.7))
    # idx.hash.08 <- apply(t(t(as.matrix(ab.hash))/colSums(ab.hash)), 2, function(col) which(col > 0.8))
    idx.cd4 <- apply(t(t(as.matrix(ab.cd4))/colSums(ab.cd4)), 2, function(col) which(col > 0.6))
    idx.cd45 <- apply(t(t(as.matrix(ab.cd45))/colSums(ab.cd45)), 2, function(col) which(col > 0.6))
    if(exists("ab.gd")){
        # ab.gd[which(ab.gd >= 10)] <- ">10reads"
        # ab.gd[which(ab.gd == "0")] <- "0reads"
        # ab.gd[which(!ab.gd %in% c(">10reads", "0reads"))] <- "1-9reads"
        idx.gdcd4 <- apply(t(t(as.matrix(ab.gdcd4))/colSums(ab.gdcd4)), 2, function(col) which(col > 0.6))
    }
    
    
    stopifnot(all(colnames(ab.aggr) == row.names(sx@meta.data)))
    stopifnot(all(colnames(ab.hash) == row.names(sx@meta.data)))
    
    aggr <- lapply(idx.aggr, function(i) row.names(ab.aggr)[i])
    aggr[sapply(aggr, length) == 0] <- "Undefined"
    
    hash <- lapply(idx.hash, function(i) row.names(ab.hash)[i])
    hash[sapply(hash, length) == 0] <- "Undefined"
    # 
    # hash55 <- lapply(idx.hash.55, function(i) row.names(ab.hash)[i])
    # hash55[sapply(hash55, length) == 0] <- "Undefined"
    # 
    # hash07 <- lapply(idx.hash.07, function(i) row.names(ab.hash)[i])
    # hash07[sapply(hash07, length) == 0] <- "Undefined"
    # 
    # hash08 <- lapply(idx.hash.08, function(i) row.names(ab.hash)[i])
    # hash08[sapply(hash08, length) == 0] <- "Undefined"
    
    cd4cd8 <- lapply(idx.cd4, function(i) row.names(ab.cd4)[i])
    cd4cd8[sapply(idx.cd4, length) == 0] <- "Undefined"
    
    cd45x <- lapply(idx.cd45, function(i) row.names(ab.cd45)[i])
    cd45x[sapply(idx.cd45, length) == 0] <- "Undefined"
    
    if(exists("ab.gd")){
        gdcd4 <- lapply(idx.gdcd4, function(i) row.names(ab.gdcd4)[i])
        gdcd4[sapply(idx.gdcd4, length) == 0] <- "Undefined"
    }
    
    sx@meta.data$hash <- unlist(hash)
    # sx@meta.data$hash55 <- unlist(hash55)
    # sx@meta.data$hash07 <- unlist(hash07)
    # sx@meta.data$hash08 <- unlist(hash08)
    if(exists("ab.gd")) sx@meta.data$gd <- ab.gd
    if(exists("ab.gd")) sx@meta.data$gd.counts <- gd.counts
    if(exists("ab.gd")) sx@meta.data$gdcd4 <- unlist(gdcd4)
    if(exists("ab.gd")) sx@meta.data$gd.ratio <- ab.gdcd4["HTO-GD.TCR",]/colSums(ab.gdcd4)
    sx@meta.data$hash.agg <- unlist(aggr)
    sx@meta.data$cd4cd8 <- unlist(cd4cd8)
    sx@meta.data$cd45x <- unlist(cd45x)
    sx@meta.data$hash.sum <- colSums(ab.hash)
    sx@meta.data$hash.ratio <- colMaxs(as.matrix(ab.hash))/colSums(ab.hash)
    sx@meta.data$hash.agg.ratio <- colMaxs(as.matrix(ab.aggr))/colSums(ab.aggr)
    sx@meta.data$CD45 <- ab['HTO-CD45.1',]/(ab['HTO-CD45.1',] + ab['HTO-CD45.2',])
    sx@meta.data$CD45.sum <- ab['HTO-CD45.1',] + ab['HTO-CD45.2',]
    sx@meta.data$CD4 <- ab['HTO-CD4',]/(ab['HTO-CD4',] + ab['HTO-CD8',])
    sx@meta.data$CD4.sum <- ab['HTO-CD4',] + ab['HTO-CD8',]
    object <- sx
    return(object)
}

for (experimentx in unique(names(pbmc_list))){
    pbmc_list[[experimentx]] <- mapply(AssignMetadata, pbmc_list[[experimentx]], CountsAB[[experimentx]])
}
for (experimentx in names(pbmc_list)){
    for(sample in names(pbmc_list[[experimentx]])){
        pbmc_list[[experimentx]][[sample]]$percent.mt <- PercentageFeatureSet(pbmc_list[[experimentx]][[sample]], pattern = "^mt-")
    }
}

# CELL CYCLE
for (experimentx in unique(names(pbmc_list))){
    for(sx in names(pbmc_list[[experimentx]])){
        # add sample name to metadata 
        pbmc_list[[experimentx]][[sx]]$sample <- sx 
        # Cell cycle scoring (will be added to metadata)
        pbmc_list[[experimentx]][[sx]] <- NormalizeData(pbmc_list[[experimentx]][[sx]], verbose = TRUE)
        pbmc_list[[experimentx]][[sx]] <- CellCycleScoring(pbmc_list[[experimentx]][[sx]], 
                                                           s.features = str_to_title(cc.genes$s.genes), 
                                                           g2m.features = str_to_title(cc.genes$g2m.genes), set.ident = TRUE)
    }
}

pbmc_list <- pbmc_list[c("HDAC1","HDAC2")]
pbmc_subset <- list()
#Filtering
for (experimentx in names(pbmc_list)){
    pbmc_subset[[experimentx]] <-
        lapply(pbmc_list[[experimentx]],
               subset,
               subset = nFeature_RNA > 50
               # &
               #     percent.mt < 10 & hash.sum > 10
               )
}
# 
# 
# pbmc_subset_raw <- list()
# for (experimentx in names(pbmc_list)){
#     pbmc_subset_raw[[experimentx]] <-
#         lapply(pbmc_list[[experimentx]],
#                subset,
#                subset = nFeature_RNA >= 0 &
#                    percent.mt <= 100)
# }
#table of filtered data
# 
# 
# filtered_data <- tibble(
#     .rows = length(unlist(pbmc_list)),
#     experiment = c(rep(names(pbmc_list)[1], length(names(
#         pbmc_list[[1]]
#     ))), rep(names(pbmc_list)[2], length(names(
#         pbmc_list[[2]]
#     )))),
#     sample = lapply(pbmc_list, names) %>% unlist() %>% unname(),
#     start = lapply(pbmc_list %>% unlist(), ncol) %>% unlist() %>% unname(),
#     "nFeature_RNA" = lapply(
#         lapply(pbmc_list %>% unlist(), subset, subset = nFeature_RNA > 200),
#         ncol
#     ) %>% unlist() %>% unname(),
#     "percent.mt" = lapply(
#         lapply(pbmc_list %>% unlist(), subset, subset = percent.mt < 10),
#         ncol
#     ) %>% unlist() %>% unname(),
#     "hash.sum" = lapply(lapply(pbmc_list %>% unlist(), subset, subset = hash.sum > 10),
#                         ncol) %>% unlist() %>% unname(),
#     "hash.ratio" = lapply(
#         lapply(pbmc_list %>% unlist(), subset, subset = hash.ratio > 0.6),
#         ncol
#     ) %>% unlist() %>% unname(),
#     "treat.ratio" = lapply(
#         lapply(pbmc_list %>% unlist(), subset, subset = hash.agg.ratio > 0.6),
#         ncol
#     ) %>% unlist() %>% unname(),
#     filtered = lapply(pbmc_subset %>% unlist(), ncol) %>% unlist() %>% unname(),
#     assigned.hash = lapply(
#         lapply(pbmc_subset %>% unlist(), subset, subset = hash.ratio > 0.6),
#         ncol
#     ) %>% unlist() %>% unname(),
#     assigned.treat = lapply(
#         lapply(pbmc_subset %>% unlist(), subset, subset = hash.agg.ratio > 0.6),
#         ncol
#     ) %>% unlist() %>% unname()
# )
# 
# gt_table <- gt(filtered_data, rowname_col = "sample", groupname_col = "experiment")
# 
# gt_table <-
#     gt_table %>%  tab_header(title = "Filtered Cells") %>%
#     tab_spanner(
#         label = "Cell Counts",
#         columns = c(start, filtered, assigned.hash, assigned.treat),
#         gather = T
#     ) %>%  tab_spanner(
#         label = "Filter Criteria",
#         columns = c(
#             "nFeature_RNA",
#             "percent.mt",
#             "hash.sum"),
#         gather = T) %>% 
#     tab_spanner(
#         label = "Hashtag Ratios",
#         columns = c(
#             "hash.ratio",
#             "treat.ratio"),
#         gather = T
#     ) %>%
#     fmt_number(use_seps = T,
#                decimals = 0) %>%
#     tab_style(style = list(cell_text(weight = "bold")),
#               locations = cells_body(columns = filtered)) %>%
#     tab_style(style = list(cell_text(weight = "bold", style = "italic")),
#               locations = cells_stub()) %>%
#     tab_style(style = cell_fill(color = "gray90"),
#               locations = list(cells_body(rows = contains(c(
#                   "fSkin", "mSkin"
#               ))),
#               cells_stub(rows = contains(c(
#                   "fSkin", "mSkin"
#               ))))) %>%
#     tab_footnote(footnote = md("Filtered with 'Filter Criteria'"),
#                  location = cells_column_labels(columns = filtered)) %>% 
#     tab_footnote(footnote = md("Filtered with 'Filter Criteria' & 'hash.ratio'"),
#                  location = cells_column_labels(columns = assigned.hash)) %>% 
#     tab_footnote(footnote = md("Filtered with 'Filter Criteria' & 'treat.ratio'"),
#                  location = cells_column_labels(columns = assigned.treat)) %>% 
#     tab_footnote(footnote = md("Cells expressing more than **200 unique genes**"),
#                  location = cells_column_labels(columns = contains("nFeat"))) %>%
#     tab_footnote(footnote = md("Cells with less than **10% mitochondrial reads**"),
#                  location = cells_column_labels(columns = contains("percent.mt"))) %>%
#     tab_footnote(footnote = md("Cells with at least **10 hashtag reads**"),
#                  location = cells_column_labels(columns = contains("hash.sum"))) %>%
#     tab_footnote(footnote = md("Cells with more than **60% reads** of **one hashtag**"),
#                  location = cells_column_labels(columns = contains("hash.ratio"))) %>%
#     tab_footnote(footnote = md("Cells with more than **60% reads** of hashtags corresponding to a **treatment**"),
#                  location = cells_column_labels(columns = contains("treat.ratio"))) %>% 
#     tab_options(row_group.as_column = T)
# 
#gt_table %>% gtsave(file.path(tablesDir,"filtered_cells.png"), expand = 10)

# monocle.obj.list <- list()


# Seurat ------------------------------------------------------------------

seurat.obj <- list()
# pbmc_subset <- pbmc_list
# monocle.obj.list <- list()
for (experimentx in names(pbmc_subset)) {
    for (sx in names(pbmc_subset[[experimentx]])) {
        seurat.obj[[sx]] <- pbmc_subset[[experimentx]][[sx]]
        # mat.use <- seurat.obj[[sx]]@assays$RNA@counts
        # monocle.obj.list[[sx]] <-
        # new_cell_data_set(expression_data = mat.use,
        # cell_metadata = seurat.obj[[sx]]@meta.data)
    }
}
# monocle.obj <- combine_cds(cds_list = monocle.obj.list, cell_names_unique = FALSE)
# rowData(monocle.obj)$gene_short_name <- row.names(rowData(monocle.obj))


seurat_analysis <- function(object){
  #Standard Workflow
  seu.obj <- object
  seu.obj <- NormalizeData(seu.obj)
  seu.obj <- FindVariableFeatures(seu.obj)
  seu.obj <- ScaleData(seu.obj, vars.to.regress = c("Phase"))
  seu.obj <- RunPCA(seu.obj)
  
  #Determine Dimensionality + Unsup Ana
  dims <- seq(1, DetermineDimensionality(seu.obj))
  seu.obj <- FindNeighbors(seu.obj, dims = dims)
  seu.obj <- FindClusters(seu.obj)
  seu.obj <- RunUMAP(seu.obj, dims = dims)

  #Doublet Identification
  sweep.list <- paramSweep_v3(seu.obj, PCs = dims, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk <- bcmvn[which.max(bcmvn$BCmetric), "pK"] %>% droplevels() %>% levels %>% as.numeric()
  nExp <- round(0.075*nrow(seu.obj@meta.data))
  seu.obj <- doubletFinder_v3(seu.obj, PCs = dims, pK = pk, nExp = nExp)
  return(seu.obj)
}

seu <- lapply(seurat.obj, seurat_analysis)
# 
# seu.obj <- Merge_Seurat_List(seurat.obj, add.cell.ids = names(seurat.obj))
# 
# new_names <-
#   colnames(seu.obj) %>% 
#   str_c(str_extract(., "^[:alpha:]+_[:alnum:]+") , sep = "_") %>% 
#   str_remove("^[:alpha:]+_[:alnum:]+_")
# 
# seu.obj <- RenameCells(seu.obj, new.names = new_names)
# 
# seu.obj <- NormalizeData(seu.obj)


monocle.raw.list <- list()
seurat.obj <- list()
pbmc_subset_raw <- pbmc_list_raw

for (experimentx in names(pbmc_subset_raw)) {
    for (sx in names(pbmc_subset_raw[[experimentx]])) {
        seurat.obj[[sx]] <- pbmc_subset_raw[[experimentx]][[sx]]
        mat.use <- seurat.obj[[sx]]@assays$RNA@counts
        monocle.raw.list[[sx]] <-
            new_cell_data_set(expression_data = mat.use,
                              cell_metadata = seurat.obj[[sx]]@meta.data)
    }
}

monocle.raw <- combine_cds(cds_list = monocle.raw.list, cell_names_unique = FALSE)
rowData(monocle.raw)$gene_short_name <- row.names(rowData(monocle.raw))

saveRDS(monocle.raw, file.path(dataDir, "scRNAseq_1_monocle_raw.cds"))

# monocle.obj@colData$hash <- as.factor(monocle.obj@colData$hash)
# monocle.obj@colData$hash55 <- as.factor(monocle.obj@colData$hash55)
# monocle.obj@colData$hash07 <- as.factor(monocle.obj@colData$hash07)
# monocle.obj@colData$hash08 <- as.factor(monocle.obj@colData$hash08)
# monocle.obj@colData$gd <- as.factor(monocle.obj@colData$gd)
# monocle.obj@colData$hash.agg <- as.factor(monocle.obj@colData$hash.agg)
# monocle.obj@colData$gdcd4 <- as.factor(monocle.obj@colData$gdcd4)
# monocle.obj@colData$cd4cd8 <- as.factor(monocle.obj@colData$cd4cd8)
# monocle.obj@colData$cd45x <- as.factor(monocle.obj@colData$cd45x)
# monocle.obj@colData$sample <- as.factor(monocle.obj@colData$sample)
# monocle.obj@colData$organ <- as.factor(str_extract(monocle.obj@colData$sample, "[A-Z][:alpha:]+"))
# levels1 <- monocle.obj@colData %>% as_tibble %>% filter(organ == "Skin") %>% pull(sample) %>% droplevels() %>%  levels
# levels2 <- monocle.obj@colData %>% as_tibble %>% filter(organ == "LN") %>% pull(sample) %>% droplevels() %>%  levels
# levelsx <- paste(c(levels1, levels2))
# monocle.obj@colData$sample <- monocle.obj@colData$sample %>% forcats::fct_relevel(levelsx)
# monocle.obj@colData$Phase <- as.factor(monocle.obj@colData$Phase)
# monocle.obj@colData$batch <-
#     as.factor(str_extract(monocle.obj@colData$sample, "[A-Z][0-9]$"))
# monocle.obj@colData$treatment <-
#     as.factor(case_when(
#         grepl("30[1-3]$", monocle.obj@colData$hash) ~ "HDAC_WT",
#         grepl("30[4-6]$", monocle.obj@colData$hash) ~ "HDAC_cKO",
#         grepl("30[7-9]$", monocle.obj@colData$hash) ~ "NoT",
#         TRUE ~ "Undefined"
#     ))
# monocle.obj@colData$treatment55 <-
#     as.factor(case_when(
#         grepl("30[1-3]$", monocle.obj@colData$hash55) ~ "HDAC_WT",
#         grepl("30[4-6]$", monocle.obj@colData$hash55) ~ "HDAC_cKO",
#         grepl("30[7-9]$", monocle.obj@colData$hash55) ~ "NoT",
#         TRUE ~ "Undefined"
#     ))
# monocle.obj@colData$treatment07 <-
#     as.factor(case_when(
#         grepl("30[1-3]$", monocle.obj@colData$hash07) ~ "HDAC_WT",
#         grepl("30[4-6]$", monocle.obj@colData$hash07) ~ "HDAC_cKO",
#         grepl("30[7-9]$", monocle.obj@colData$hash07) ~ "NoT",
#         TRUE ~ "Undefined"
#     ))
# monocle.obj@colData$treatment08 <-
#     as.factor(case_when(
#         grepl("30[1-3]$", monocle.obj@colData$hash08) ~ "HDAC_WT",
#         grepl("30[4-6]$", monocle.obj@colData$hash08) ~ "HDAC_cKO",
#         grepl("30[7-9]$", monocle.obj@colData$hash08) ~ "NoT",
#         TRUE ~ "Undefined"
#     ))
# monocle.obj@colData$treatment.agg <-
#     as.factor(case_when(
#         grepl("\\-3$", monocle.obj@colData$hash.agg) ~ "NoT",
#         grepl("\\-1$", monocle.obj@colData$hash.agg) ~ "HDAC_WT",
#         grepl("\\-2$", monocle.obj@colData$hash.agg) ~ "HDAC_cKO",
#         TRUE ~ "Undefined"
#     ))
# monocle.obj@colData$hashid <-
#     as.factor(ifelse(
#         monocle.obj@colData$hash == "Undefined",
#         "Undefined",
#         str_extract(monocle.obj@colData$hash, "[0-9]+$")
#     ))
# monocle.obj@colData$hashid.agg <-
#     as.factor(ifelse(
#         monocle.obj@colData$hash.agg == "Undefined",
#         "Undefined",
#         str_extract(monocle.obj@colData$hash.agg, "[0-9]+$")
#     ))
# monocle.obj@colData$sex <- as.factor(str_to_upper(str_extract(monocle.obj@colData$sample,"^[A-z]")))
# monocle.obj@colData$treatment <-
#     monocle.obj@colData$treatment %>% forcats::fct_relevel("NoT", "HDAC_WT", "HDAC_cKO", "Undefined")
# monocle.obj@colData$treatment07 <-
#     monocle.obj@colData$treatment07 %>% forcats::fct_relevel("NoT", "HDAC_WT", "HDAC_cKO", "Undefined")
# monocle.obj@colData$treatment08 <-
#     monocle.obj@colData$treatment08 %>% forcats::fct_relevel("NoT", "HDAC_WT", "HDAC_cKO", "Undefined")
# monocle.obj@colData$treatment55 <-
#     monocle.obj@colData$treatment55 %>% forcats::fct_relevel("NoT", "HDAC_WT", "HDAC_cKO", "Undefined")
# monocle.obj@colData$treatment.agg <-
#     monocle.obj@colData$treatment.agg %>% forcats::fct_relevel("NoT", "HDAC_WT", "HDAC_cKO", "Undefined")
# monocle.obj@colData$sample_treat <- as.factor(paste(monocle.obj@colData$sample, monocle.obj@colData$treatment, sep = "___"))
# monocle.obj@colData$sample_treat.agg <-
#     as.factor(paste(
#         monocle.obj@colData$sample,
#         monocle.obj@colData$treatment.agg,
#         sep = "___"
#     ))
# monocle.obj@colData$experiment <-
#     factor(case_when(
#         !grepl("[0-9]+B[0-9]", monocle.obj@colData$sample) ~ "HDAC1",
#         grepl("40", monocle.obj@colData$sample) ~ "HDAC1",
#         TRUE ~ "HDAC2"
#     ))
# monocle.obj@colData$mouse <-
#     paste0(
#         "M",
#         monocle.obj@colData$sample %>%  str_remove_all("Skin|LN") %>% forcats::fct_inorder() %>% dense_rank()
#     )
# 
# monocle.obj@colData %>% as_tibble() %>% filter(!treatment.agg %in% c("Undefined")) %>% dplyr::select(batch, experiment, sex, treatment.agg) %>% unique()
# 
# 
# samples_run1 <- str_remove(grep(
#     list.files(file.path(
#         dirname(dirname(rawDir)), "scRNA_from_BSF1_2022_11/COUNT"
#     )),
#     pattern = "AGG",
#     invert = T,
#     value = T
# ), "_tra.*")
# 
# samples_run2 <- str_remove(grep(
#     list.files(file.path(
#         dirname(dirname(rawDir)), "scRNA_from_BSF_2023_03_01/COUNT"
#     )),
#     pattern = "AGG",
#     invert = T,
#     value = T
# ), "_tra.*")
# 
# samples_run2 <- samples_run2[!samples_run2 %in% samples_run1]
# 
# samples_run3 <- str_remove(grep(
#     list.files(rawDir),
#     pattern = "AGG",
#     invert = T,
#     value = T
# ), "_tra.*")
# 
# samples_run3 <- samples_run3[!samples_run3 %in% c(samples_run1, samples_run2)]
# 
# monocle.obj@colData$run <-
#     case_when(monocle.obj@colData$sample %in% samples_run1 ~ "run1",
#               monocle.obj@colData$sample %in% samples_run2 ~ "run2",
#               monocle.obj@colData$sample %in% samples_run3 ~ "run3")

# additional.meta <- monocle.obj@colData %>% subset(sample != "fLN_40B3", c(treatment55, hash55))
# saveRDS(additional.meta, file.path(dataDir,"additional_meta.RDS"))


# Process dataset
monocle.obj <-
  preprocess_cds(monocle.obj, verbose = TRUE) %>%
  reduce_dimension(preprocess_method = "PCA", verbose = TRUE)

monocle.obj <- align_cds(monocle.obj,
                         alignment_group = "sample", 
                         residual_model_formula_str = "~Phase",
                         verbose = TRUE
)

monocle.obj <- reduce_dimension(monocle.obj,
                                reduction_method = "UMAP",
                                preprocess_method = "Aligned",
                                verbose = TRUE)

# Clustering
monocle.obj <-  cluster_cells(monocle.obj)

monocle.obj@colData$Cluster <- unname(clusters(monocle.obj[,rownames(colData(monocle.obj))]))

# Store old dataset
# saveRDS(monocle.obj, file = file.path(dataDir, "scRNAseq_1_monocle_old.cds"))

as.Seurat(monocle.obj, counts = "counts", logcounts = NULL)

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


# Old CT Assignment -------------------------------------------------------

library(renv)
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
library(patchwork)
library(scRNAseq)
library(data.table)
library(readr)
library(forcats)

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

# Load Data ---------------------------------------------------------------

# monocle.obj <- readRDS(file.path(dataDir, "/scRNAseq_1_monocle.cds"))

reference.obj <- readRDS(file.path(dataDir, "TwoGroups_celltypes_group.rds"))
reference.obj <- UpdateSeuratObject(object = reference.obj)
reference.obj <- as.SingleCellExperiment(reference.obj)

colnames(colData(reference.obj))[colnames(colData(reference.obj)) == "celltype"] <- "label.fine"


# SingleR annotation ------------------------------------------------------

## Reference data ----------------------------------------------------------

### making a list of reference data sets  ----------------------------------------------

ref_data <- list(
    
    #Human Data Sets: 
    ##HPCA generally / Keratinocytes
    hpca = HumanPrimaryCellAtlasData(),
    
    ##DB ImmuneCells, comprehensive CD4+ subsets; only one B cell subset, no dendritic cells
    dice = DatabaseImmuneCellExpressionData(),
    
    ##Monaco Immunc cells
    monaco = MonacoImmuneData(),
    
    #MouseData:
    ##MouseRNAseqData
    mrsd = MouseRNAseqData(),
    
    ##ImmGen 
    immgen = ImmGenData(),
    
    ##reference Dataset
    struc = reference.obj
)


#Human Gene names as Mouse Gene names

rownames(ref_data$hpca) <- str_to_title(rownames(ref_data$hpca))
rownames(ref_data$dice) <- str_to_title(rownames(ref_data$dice))
rownames(ref_data$monaco) <- str_to_title(rownames(ref_data$monaco))


#
metadata <- as.data.frame(colData(monocle.obj))

counts_matrix <- assay(monocle.obj)

obj.as.seurat <- CreateSeuratObject(counts = counts_matrix,
                                    project = "integrated_scRNAseq",
                                    assay = "integrated",
                                    meta.data = metadata)


sce <-  as.SingleCellExperiment(obj.as.seurat)



### run singleR with additional saving of original results data -------------------------------

# parallel computation
doParallel::registerDoParallel(cores=7)
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
        # saveRDS(results, file.path(dataDir, ref.file.orig))
        
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
               readRDS(file.path(dataDir, paste0("orig_cell_types_", ref, "_", labelx, ".rds"))))
        
    }
    
}

### assign labels to colData ----------------------------------------------------

# ref based labels

for(variable in ls(pattern = "^cell_types")){
    colData(monocle.obj)[[variable]] <- get(variable)$labels
}

monocle.obj@colData$celltype_raw <- as.factor(
    case_when(
        monocle.obj@colData$cell_types_immgen_label.main == "Epithelial cells" ~ "Keratinocytes",
        TRUE ~ monocle.obj@colData$cell_types_immgen_label.main
    )
)

ct_majority <-
    monocle.obj@colData %>% as_tibble() %>% 
    group_by(Cluster, celltype_raw) %>%  
    summarize(n = n()) %>% 
    slice_max(n = 1, order_by = n, with_ties = F)

for(r in seq(nrow(monocle.obj@colData))) {
    monocle.obj@colData$ct_cluster[r] <-
        paste0(
            as.character(monocle.obj@colData$Cluster[r]),
            " (",
            as.character(ct_majority[as.numeric(ct_majority$Cluster) == as.numeric(monocle.obj@colData$Cluster)[r], ]$celltype_raw),
            ")"
        )
}

monocle.obj@colData$celltype <-
    monocle.obj@colData$ct_cluster %>% str_extract_all("\\w+( \\w+)|[:alpha:]+") %>% unlist()

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
freq_list <- freq_list %>% arrange(celltype, Cluster)
freq_list$ct_cluster <- paste0(freq_list$Cluster," (",freq_list$celltype,")")
freq_list$ct_cluster <- as.factor(freq_list$ct_cluster)
freq_list$ct_cluster <- fct_inorder(freq_list$ct_cluster)
# saveRDS(freq_list, file.path(dataDir,"all_labels_ref_based_for_HM.rds"))

# plot a heatmap
freq_list %>% filter(
    reference %in% c(
        "hpca_label.fine",
        "immgen_label.fine",
        "immgen_label.main",
        "mrsd_label.fine",
        "mrsd_label.main",
        "struc_label.fine"
    )
) %>%
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

# ggsave(file.path(plotsDir, "all_labels_HM.pdf"), height = 15, width = 30)


# Save Object -------------------------------------------------------------
monocle.obj@colData$ct_cluster <-
    factor(monocle.obj@colData$ct_cluster,
           levels = levels(freq_list$ct_cluster))

# saveRDS(monocle.obj, file.path(dataDir, "scRNAseq_2_monocle.cds"))


# Old CLuster refinement --------------------------------------------------

library(tidyr)
library(dplyr)
library(stringr)

new.clusters <- read.csv(file.path(tablesDir, "NewClusterML.csv"))
laia.clusters <- read.csv(file.path(tablesDir, "Ludwig4.csv"))

new.clusters <- new.clusters %>% filter(NewClusterML != "")

new.clusters <- new.clusters %>% mutate(cluster = as.numeric(str_extract(new.clusters$NewClusterML, "\\d+")))

new.clusters$Barcode <- case_when(grepl("\\-1$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-1", "-1_fLN_40B2"),
                                  grepl("\\-2$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-2", "-1_fLN_40B3"),
                                  grepl("\\-3$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-3", "-1_fLN_41B1"),
                                  grepl("\\-4$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-4", "-1_fLN_41B2"),
                                  grepl("\\-5$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-5", "-1_fLN_41B3"),
                                  grepl("\\-6$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-6", "-1_fLN_B1"),
                                  grepl("\\-7$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-7", "-1_fSkin_40B2"),
                                  grepl("\\-8$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-8", "-1_fSkin_40B3"),
                                  grepl("\\-9$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-9", "-1_fSkin_41B1"),
                                  grepl("\\-10$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-10", "-1_fSkin_41B2"),
                                  grepl("\\-11$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-11", "-1_fSkin_41B3"),
                                  grepl("\\-12$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-12", "-1_fSkin_B1"),
                                  grepl("\\-13$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-13", "-1_mLN_40B2"),
                                  grepl("\\-14$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-14", "-1_mLN_40B3"),
                                  grepl("\\-15$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-15", "-1_mLN_41B1"),
                                  grepl("\\-16$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-16", "-1_mLN_41B3"),
                                  grepl("\\-17$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-17", "-1_mLN_41B4"),
                                  grepl("\\-18$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-18", "-1_mLN_B1"),
                                  grepl("\\-19$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-19", "-1_mSkin_40B2"),
                                  grepl("\\-20$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-20", "-1_mSkin_40B3"),
                                  grepl("\\-21$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-21", "-1_mSkin_41B1"),
                                  grepl("\\-22$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-22", "-1_mSkin_41B3"),
                                  grepl("\\-23$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-23", "-1_mSkin_B1"))




clusters <- setNames(new.clusters$cluster, new.clusters$Barcode)
clusters <- clusters[colnames(monocle.obj)]
clusters <- as.factor(clusters)

clusters <- clusters %>% recode("57" = "34",
                                "42" = "14",
                                "51" = "39",
                                "9" = "7",
                                "62" = "47",
                                "40" = "12",
                                "43" = "19",
                                "27" = "26",
                                "18" = "13",
                                "30" = "2",
                                "3" = "1",
                                "5" = "1", 
                                "6" = "1",
                                "16" = "1",
                                "22" = "1"
)

#Just to check which clusters got removed and need replacement in order to have the clusters in a sequence
ol_c <- clusters[which(!(clusters %in% seq(length(levels(clusters)))))] %>% droplevels %>% levels
new_c <- seq(length(levels(clusters)))[which(!(seq(length(levels(clusters))) %in% clusters))]

ol_c
new_c

clusters <- clusters %>% recode_factor(
    "53" = "3",
    "54" = "5",
    "55" = "6",
    "56" = "9",
    "58" = "16",
    "59" = "18",
    "60" = "22",
    "61" = "27",
    "63" = "30",
    "64" = "40",
    "65" = "42",
    "66" = "43",
    "67" = "51"
)

#having levels in right order
for (i in as.character(length(levels(clusters)):1)){
    clusters <- relevel(clusters, i)
}

laia.clusters <- laia.clusters %>% filter(Ludwig3 != "")
laia.clusters <- laia.clusters %>% mutate(cluster = as.numeric(str_extract(laia.clusters$Ludwig3, "\\d+")))

laia.clusters$Barcode <- case_when(grepl("\\-1$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-1", "-1_fLN_40B2"),
                                   grepl("\\-2$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-2", "-1_fLN_40B3"),
                                   grepl("\\-3$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-3", "-1_fLN_41B1"),
                                   grepl("\\-4$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-4", "-1_fLN_41B2"),
                                   grepl("\\-5$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-5", "-1_fLN_41B3"),
                                   grepl("\\-6$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-6", "-1_fLN_B1"),
                                   grepl("\\-7$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-7", "-1_fSkin_40B2"),
                                   grepl("\\-8$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-8", "-1_fSkin_40B3"),
                                   grepl("\\-9$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-9", "-1_fSkin_41B1"),
                                   grepl("\\-10$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-10", "-1_fSkin_41B2"),
                                   grepl("\\-11$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-11", "-1_fSkin_41B3"),
                                   grepl("\\-12$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-12", "-1_fSkin_B1"),
                                   grepl("\\-13$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-13", "-1_mLN_40B2"),
                                   grepl("\\-14$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-14", "-1_mLN_40B3"),
                                   grepl("\\-15$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-15", "-1_mLN_41B1"),
                                   grepl("\\-16$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-16", "-1_mLN_41B3"),
                                   grepl("\\-17$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-17", "-1_mLN_41B4"),
                                   grepl("\\-18$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-18", "-1_mLN_B1"),
                                   grepl("\\-19$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-19", "-1_mSkin_40B2"),
                                   grepl("\\-20$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-20", "-1_mSkin_40B3"),
                                   grepl("\\-21$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-21", "-1_mSkin_41B1"),
                                   grepl("\\-22$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-22", "-1_mSkin_41B3"),
                                   grepl("\\-23$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-23", "-1_mSkin_B1"))



# setdiff(
#     clusters %>% unname %>% unique %>% sort,
#     laia.clusters %>% filter(!is.na(cluster)) %>% pull(cluster) %>% unique %>% sort
# )

#laia.clusters[which(laia.clusters$Barcode %in% names((clusters[clusters == 3]))),] %>% pull(Ludwig3) %>% unique

laia.clusters$cluster <- case_when(laia.clusters$Ludwig3 == "GC B cells" ~ 50,
                                   laia.clusters$Ludwig3 == "Plasma cells" ~ 3,
                                   laia.clusters$Ludwig3 == "Pre-plasmablasts" ~ 53,
                                   TRUE ~ laia.clusters$cluster)


clusters <- setNames(laia.clusters$cluster, laia.clusters$Barcode)
clusters <- clusters[colnames(monocle.obj)]
clusters <- as.factor(clusters)

# Assigning to monocle.obj ------------------------------------------------
monocle.obj@clusters$UMAP$clusters <- clusters

monocle.obj@colData$Cluster <- unname(clusters(monocle.obj[,rownames(colData(monocle.obj))]))



monocle.obj@colData$celltype <- case_when(monocle.obj@colData$Cluster == 34 ~ "CD45.1+ T cells",
                                          monocle.obj@colData$Cluster == 44 ~ "Unassigned1",
                                          monocle.obj@colData$Cluster == 11 ~ "Unassigned2",
                                          monocle.obj@colData$Cluster == 5  ~ "Unassigned3",
                                          monocle.obj@colData$Cluster == 18 ~ "Myofibroblasts",
                                          monocle.obj@colData$Cluster == 50 ~ "GC B cells",
                                          monocle.obj@colData$Cluster == 3 ~ "Plasma cells",
                                          monocle.obj@colData$Cluster == 53 ~ "Pre-plasmablasts",
                                          TRUE ~ monocle.obj@colData$celltype)



monocle.obj@colData <- transform(
    monocle.obj@colData,
    ct_cluster = paste0(
        as.character(monocle.obj@colData$Cluster),
        " (",
        monocle.obj@colData$celltype
        ,
        ")"
    )
)

ordered_factor <-
    monocle.obj@colData %>% as_tibble() %>% arrange(celltype, Cluster) %>% select(celltype, Cluster)  %>%
    mutate(ct_cl = forcats::fct_inorder(paste0(as.character(Cluster),
                                               " (",
                                               celltype
                                               ,
                                               ")")))

monocle.obj@colData$ct_cluster <- factor(monocle.obj@colData$ct_cluster, levels = levels(ordered_factor$ct_cl))

# Save Object -------------------------------------------------------------

saveRDS(monocle.obj, file.path(dataDir, "scRNAseq_2a_monocle_with_fLN_40B3.cds"))

monocle.obj <- monocle.obj[,colnames(monocle.obj)[!grepl("fLN_40B3", colnames(monocle.obj))]]
monocle.obj@colData$sample <- monocle.obj@colData$sample %>% droplevels()

# additional.meta <- readRDS(file.path(dataDir, "additional_meta.RDS"))
# monocle.obj@colData <- cbind(monocle.obj@colData,additional.meta)

# saveRDS(monocle.obj, file.path(dataDir, "scRNAseq_3_monocle_more_hash_cutoff.cds"))
