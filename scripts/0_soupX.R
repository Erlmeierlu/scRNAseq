library(dplyr)
library(stringr)
library(Seurat)
library(DoubletFinder)
library(scCustomize)
library(data.table)
library(monocle3)
library(SoupX)
library(foreach)
library(readr)

#setting up directories
gfsDir <- 'media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin'
rawDir <- file.path(gfsDir, 'raw_data/scRNA_from_BSF/COUNT')
plotsDir <- file.path(gfsDir, 'plots')
tablesDir <- file.path(gfsDir, 'tables')
oldDir <- "/vscratch/scRNAseq/data/old"
shinyDir <- 'dge-app'
dataDir <-"data"
resDir <- "results"

#directories to loop over
fileDirs <- list()
fileDirs[["filtered"]] <- paste0(
    grep(list.files(rawDir), pattern = "AGG", invert = T, value = T),
    "/filtered_feature_bc_matrix/"
)
fileDirs[["raw"]] <- paste0(
    grep(list.files(rawDir), pattern = "AGG", invert = T, value = T),
    "/raw_feature_bc_matrix/"
)

fileDirs <- lapply(fileDirs, function(x) {
    x[!grepl("fLN_40B3", x)]
})


# Load Raw Data -----------------------------------------------------------
load_data <- function(dat, 
                      path, 
                      slot = c("none", "gene expression", "antibody capture")){
    data <- Read10X(file.path(path, dat))
    
    slot <- match.arg(slot)
    if (slot == "gene expression")
        data <- data[[1]]
    else if (slot == "antibody capture")
        data <- data[[2]]
    
    obj <- CreateSeuratObject(
        counts = data,
        project = "scRNAseq",
        min.cells = 0,
        min.features = 0
    )
    x <- dat %>% str_extract("^[:alpha:]+_[:alnum:]+")
    
    setNames(list(obj), x)
}

data_list <- sapply(fileDirs$filtered, 
                    load_data,
                    path = rawDir, 
                    slot = "gene expression", 
                    USE.NAMES = FALSE)

ab_list <- sapply(fileDirs$filtered, 
                  load_data,
                  path = rawDir, 
                  slot = "antibody capture", 
                  USE.NAMES = FALSE)


data_list_raw <- sapply(fileDirs$raw, 
                        load_data,
                        path = rawDir, 
                        slot = "gene expression", 
                        USE.NAMES = FALSE)


#Create AB counts list
CountsAB <- lapply(ab_list, GetAssayData, slot = "counts")

write_rds(CountsAB, file.path(dataDir, "ab_counts.rds"))

# Pre-Processing ---------------------------------------------------------

# AssignMetadata <- function(sx, abcounts) {
#     stopifnot(all(colnames(sx) == colnames(abcounts)))
#     
#     ab.hash <- abcounts[grepl("-[0-9]+$", row.names(abcounts)), ]
#     ab.aggr <-
#         list(
#             "HTO-fLN-40B2-1" = 1:3,
#             "HTO-fLN-40B2-2" = 4:6,
#             "HTO-fLN-40B2-3" = 7:9
#         ) %>%
#         sapply(function(x)
#             colSums(ab.hash[x, ])) %>%
#         t() %>%
#         as("sparseMatrix")
#     
#     stopifnot(colSums(ab.aggr) == colSums(ab.hash))
#     
#     ab.cd4 <- abcounts[grepl("CD[1-9]$", row.names(abcounts)), ]
#     ab.cd45 <- abcounts[grepl("CD45.*$", row.names(abcounts)), ]
#     
#     if (any(grepl("GD", rownames(abcounts)))) {
#         ab.gd <- abcounts[grepl("GD", rownames(abcounts)), ]
#         ab.gdcd4 <- abcounts[grepl("GD|CD4$", rownames(abcounts)), ]
#     } else {
#         ab.gd <- NULL
#         ab.gdcd4 <- NULL
#     }
#     
#     assign_ab <- function(ab, threshold) {
#         ratio <- t(t(as.matrix(ab)) / colSums(ab))
#         #which ab fulfills the threshold criteria in each cell?
#         apply(ratio, 2, function(col) {
#             id <- which(col > threshold)
#             if (length(id) == 0)
#                 return("Undefined")
#             names(id)
#         })
#     }
#     
#     aggr <- assign_ab(ab.aggr, 0.6)
#     hash <- assign_ab(ab.hash, 0.6)
#     cd4cd8 <- assign_ab(ab.cd4, 0.75)
#     cd45x <- assign_ab(ab.cd45, 0.6)
#     if (!is.null(ab.gd)) {
#         gdcd4 <- assign_ab(ab.gdcd4, 0.75)
#     } else {
#         gdcd4 <- NULL
#     }
#     
#     max_ratio <- function(x) {
#         colMaxs(as.matrix(x)) / colSums(x)
#     }
#     
#     sx@meta.data <- sx@meta.data %>%
#         mutate(
#             hash,
#             aggr,
#             cd4cd8,
#             cd45x,
#             gdcd4 = if (!is.null(ab.gd)) {
#                 gdcd4
#             } else
#                 NA,
#             gd.ratio = if (!is.null(ab.gd)) {
#                 ab.gdcd4["HTO-GD.TCR",] / colSums(ab.gdcd4)
#             } else
#                 NA,
#             hash.sum = colSums(ab.hash),
#             hash.ratio = max_ratio(ab.hash),
#             hash.agg.ratio = max_ratio(ab.aggr),
#             CD45 = ab.cd45["HTO-CD45.1",] / colSums(ab.cd45),
#             CD45.sum = colSums(ab.cd45),
#             CD4 = ab.cd4['HTO-CD4',] / colSums(ab.cd4),
#             CD4.sum = colSums(ab.cd4)
#         )
#     sx
# }
# 
# #Assign meta
# data_list <- mapply(AssignMetadata, data_list, CountsAB)
# 

#exclude barcodes with less than 50 UMIs for soupx analysis
data_subset <-
    lapply(data_list,
           subset,
           subset = nFeature_RNA > 50)


seurat_analysis <- function(object){
    
    DetermineDimensionality <- function(object){
        pct <- object[["pca"]]@stdev / sum(object[["pca"]]@stdev) * 100
        cumu <- cumsum(pct)
        co1 <- which(cumu > 90 & pct < 5)[1]
        co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
        pcs <- min(co1, co2)
        return(pcs)
    }
    
    #Standard Workflow
    seu.obj <- object
    seu.obj <- NormalizeData(seu.obj)
    seu.obj <- FindVariableFeatures(seu.obj)
    seu.obj <- ScaleData(seu.obj)
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
    pk <- bcmvn[which.max(bcmvn$BCmetric), "pK"] %>% 
        droplevels() %>% 
        levels %>% 
        as.numeric()
    nExp <- round(0.075*nrow(seu.obj@meta.data))
    seu.obj <- doubletFinder_v3(seu.obj, PCs = dims, pK = pk, nExp = nExp)
    
    colnames(seu.obj@meta.data)[grepl(
        "DF.class", 
        colnames(seu.obj@meta.data))] <- "doublet_classification"
    
    colnames(seu.obj@meta.data)[grepl(
        "pANN", 
        colnames(seu.obj@meta.data))] <- "pANN"
    
    seu.obj
}

seu <- lapply(data_subset, seurat_analysis)

seu <- Merge_Seurat_List(seu, add.cell.ids = names(seu))

new_names <-
    colnames(seu) %>%
    str_c(str_extract(., "^[:alpha:]+_[:alnum:]+") , sep = "_") %>%
    str_remove("^[:alpha:]+_[:alnum:]+_")

seu <- RenameCells(seu, new.names = new_names)

# Create Objects --------------------------------------------------
#Save Seurat object
write_rds(seu, file.path(dataDir, '0_souped.rds'))

#Create Lookup-Table for later use in QC
doublet_LUT <- seu@meta.data[, "doublet_classification", drop = F]
doublet_LUT <- as.data.table(doublet_LUT, keep.rownames = T) 

write_rds(doublet_LUT, file.path(dataDir, "doublet_LUT.rds"))


#Raw reads as monocle object
monocle.raw.list <- lapply(data_list_raw, function(x){
    mat <- x@assays$RNA@counts
    cds <- new_cell_data_set(expression_data = mat,
                             cell_metadata = x@meta.data)
    return(cds)
    
})

monocle.raw <- combine_cds(cds_list = monocle.raw.list, cell_names_unique = FALSE)
rowData(monocle.raw)$gene_short_name <- row.names(rowData(monocle.raw))

write_rds(monocle.raw, file.path(dataDir, "0_monocle_raw.cds"))

# Load Data ---------------------------------------------------------------
#remove comment when needed
seu <- read_rds(file.path(dataDir, "0_souped.rds"))
monocle.raw <- read_rds(file.path(dataDir, "0_monocle_raw.cds"))

# SoupX -------------------------------------------------------------------

sample <- str_extract(
    colnames(seu), "[:alpha:]+_[:alnum:]+$"
    ) %>% unique()

doParallel::registerDoParallel(cores=6)
foreach(sx = sample) %dopar% {
    
sc <- SoupChannel(counts(monocle.raw[, grepl(sx, colnames(monocle.raw))]), 
                  GetAssayData(seu[, grepl(sx, colnames(seu))], slot = 'counts'))

sc <- setClusters(sc, Idents(seu[, grepl(sx, colnames(seu))]))

#for samples with few good markers; Just applies to one skin sample with few cells
mar <- quickMarkers(sc$toc,clusters = sc$metaData$clusters)
#to have at least 10 markers with additional 0.03 tfidf buffer
tfidf_min <- mar$tfidf[which(rank(-mar$tfidf) == 10)] - 0.03

#in most cases tfidfMin will be 1 (the default value)
sc <- autoEstCont(sc, tfidfMin = min(tfidf_min, 1))

out <- adjustCounts(sc)

DropletUtils:::write10xCounts(file.path(vDir, paste0("data/countsSoupX/", sx)), 
                              x = out, 
                              type = "sparse", 
                              version = "3")
}
