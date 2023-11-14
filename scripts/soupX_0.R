library(dplyr)
library(stringr)
library(Seurat)
library(DoubletFinder)
library(scCustomize)
library(data.table)
library(monocle3)
library(SoupX)
library(foreach)

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


# Load Raw Data -----------------------------------------------------------

#Seurat Object
pbmc_data_list <- list()
pbmc_list <- list()
ab_list <- list()

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
    
    ab_list[[experiment]][[x]] <- CreateSeuratObject(
        counts = pbmc_data_list[[experiment]][[x]]$'Antibody Capture',
        project = "scRNAseq",
        min.cells = 0,
        min.features = 0
    )
}

pbmc_data_list_raw <- list()
pbmc_list_raw <- list()
for (file in fileDirs[["raw"]]) {
    x <- file %>% str_remove("_trans.*")
    
    experiment <- case_when(!grepl("[0-9]+B[0-9]", x) ~ "HDAC1",
                            grepl("40", x) ~ "HDAC1",
                            TRUE ~ "HDAC2")
    
    pbmc_data_list_raw[[experiment]][[x]] <-
        Read10X(file.path(rawDir, file))
    
    pbmc_list_raw[[experiment]][[x]] <-
        CreateSeuratObject(
            counts = pbmc_data_list_raw[[experiment]][[x]]$'Gene Expression',
            project = "scRNAseq",
            min.cells = 0,
            min.features = 0
        )
}

CountsAB <- list()
#creating list of AB counts
for (experiment in names(ab_list)){
    CountsAB[[experiment]] <- lapply(ab_list[[experiment]], GetAssayData, slot = "counts")
}

#removal of bad sample
CountsAB <- lapply(CountsAB, function(x) {
    x[!grepl("fLN_40B3", names(x))]
})

saveRDS(CountsAB, file.path(dataDir, "ab_counts.rds"))


# Pre-Processing ---------------------------------------------------------

AssignMetadata <- function(object, abcounts){
    sx <- object
    barcodes <- row.names(sx@meta.data) %>% str_extract("[:alpha:]+\\-\\d")
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
    ab.aggr <- as(ab.aggr, "sparseMatrix")
    
    stopifnot(colSums(ab.aggr) == colSums(ab.hash))
    
    ab.cd4 <- ab[grepl("CD[1-9]$", row.names(ab)),]
    ab.cd45 <- ab[grepl("CD45.*$", row.names(ab)),]
    
    if(any(grepl("GD", rownames(ab)))){
        ab.gd <- ab[grepl("GD", rownames(ab)),]
        gd.counts <- ab.gd
        ab.gdcd4 <- ab[grepl("GD|CD4$", rownames(ab)),]
    } else {
        ab.gd <- NULL
        gd.counts <- NULL
        ab.gdcd4 <- NULL
    }
    
    idx.aggr <- apply(t(t(as.matrix(ab.aggr))/colSums(ab.aggr)), 2, function(col) which(col > 0.6))
    idx.hash <- apply(t(t(as.matrix(ab.hash))/colSums(ab.hash)), 2, function(col) which(col > 0.6))
    idx.cd4 <- apply(t(t(as.matrix(ab.cd4))/colSums(ab.cd4)), 2, function(col) which(col > 0.75))
    idx.cd45 <- apply(t(t(as.matrix(ab.cd45))/colSums(ab.cd45)), 2, function(col) which(col > 0.6))
    if(!is.null(ab.gd)){
        idx.gdcd4 <- apply(t(t(as.matrix(ab.gdcd4))/colSums(ab.gdcd4)), 2, function(col) which(col > 0.75))
    } else {
        idx.gdcd4 <- NULL
    }
    
    
    stopifnot(all(colnames(ab.aggr) == row.names(sx@meta.data) %>% str_extract("[:alpha:]+\\-\\d")))
    stopifnot(all(colnames(ab.hash) == row.names(sx@meta.data) %>% str_extract("[:alpha:]+\\-\\d")))
    
    aggr <- lapply(idx.aggr, function(i) row.names(ab.aggr)[i])
    aggr[sapply(aggr, length) == 0] <- "Undefined"
    
    hash <- lapply(idx.hash, function(i) row.names(ab.hash)[i])
    hash[sapply(hash, length) == 0] <- "Undefined"
    
    cd4cd8 <- lapply(idx.cd4, function(i) row.names(ab.cd4)[i])
    cd4cd8[sapply(idx.cd4, length) == 0] <- "Undefined"
    
    cd45x <- lapply(idx.cd45, function(i) row.names(ab.cd45)[i])
    cd45x[sapply(idx.cd45, length) == 0] <- "Undefined"
    
    if(!is.null(ab.gd)){
        gdcd4 <- lapply(idx.gdcd4, function(i) row.names(ab.gdcd4)[i])
        gdcd4[sapply(idx.gdcd4, length) == 0] <- "Undefined"
    }else {
        gdcd4 <- NULL
    }
    
    sx@meta.data$hash <- unlist(hash)
    if(!is.null(ab.gd)) sx@meta.data$gdcd4 <- unlist(gdcd4)
    if(!is.null(ab.gd)) sx@meta.data$gd.ratio <- ab.gdcd4["HTO-GD.TCR",]/colSums(ab.gdcd4)
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

#Assign meta
for (experimentx in unique(names(pbmc_list))){
    pbmc_list[[experimentx]] <- mapply(AssignMetadata, pbmc_list[[experimentx]], CountsAB[[experimentx]])
}
#Determine %mito genes
for (experimentx in names(pbmc_list)){
    for(sample in names(pbmc_list[[experimentx]])){
        pbmc_list[[experimentx]][[sample]]$percent.mt <- PercentageFeatureSet(pbmc_list[[experimentx]][[sample]], pattern = "^mt-")
    }
}
#cell cycle scoring
for (experimentx in unique(names(pbmc_list))){
    for(sx in names(pbmc_list[[experimentx]])){
        
        pbmc_list[[experimentx]][[sx]]$sample <- sx 
        pbmc_list[[experimentx]][[sx]] <- NormalizeData(pbmc_list[[experimentx]][[sx]], verbose = TRUE)
        pbmc_list[[experimentx]][[sx]] <- CellCycleScoring(pbmc_list[[experimentx]][[sx]], 
                                                           s.features = str_to_title(cc.genes$s.genes), 
                                                           g2m.features = str_to_title(cc.genes$g2m.genes), set.ident = TRUE)
    }
}

#exclude barcodes with less than 50 UMIs for soupx analysis
pbmc_list <- pbmc_list[c("HDAC1","HDAC2")]
pbmc_subset <- list()
for (experimentx in names(pbmc_list)){
    pbmc_subset[[experimentx]] <-
        lapply(pbmc_list[[experimentx]],
               subset,
               subset = nFeature_RNA > 50
        )
}

#unlist
seu <- list()
for (experimentx in names(pbmc_subset)) {
    for (sx in names(pbmc_subset[[experimentx]])) {
        seu[[sx]] <- pbmc_subset[[experimentx]][[sx]]
    }
}

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

seu <- lapply(seu, seurat_analysis)
seu <- lapply(seu, function(x){
    colnames(x@meta.data)[grepl("DF.class", colnames(x@meta.data))] <- "doublet_classification"
    colnames(x@meta.data)[grepl("pANN", colnames(x@meta.data))] <- "pANN"
    return(x)
})

seu <- Merge_Seurat_List(seu, add.cell.ids = names(seu))

new_names <-
    colnames(seu) %>%
    str_c(str_extract(., "^[:alpha:]+_[:alnum:]+") , sep = "_") %>%
    str_remove("^[:alpha:]+_[:alnum:]+_")

seu <- RenameCells(seu, new.names = new_names)


# Create Monocle Objects --------------------------------------------------
#Create monocle object from seurat object
monocle.obj <- new_cell_data_set(expression_data = seu@assays$RNA@counts, 
                                 cell_metadata = seu@meta.data)
rowData(monocle.obj)$gene_short_name <- row.names(rowData(monocle.obj))


#Create Lookup-Table for later use in QC
doublet_LUT <- monocle.obj@colData[, "doublet_classification", drop = F]
doublet_LUT <- as.data.table(doublet_LUT, keep.rownames = T) 

saveRDS(doublet_LUT, file.path(dataDir, "doublet_LUT.rds"))

# Process dataset again (I like monocle more here)
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
#assign metadata column
monocle.obj@colData$Cluster <- unname(clusters(monocle.obj[,rownames(colData(monocle.obj))]))

#Raw reads as monocle object
monocle.raw.list <- list()
for (experimentx in names(pbmc_list_raw)) {
    for (sx in names(pbmc_list_raw[[experimentx]])) {
        mat.use <- pbmc_list_raw[[experimentx]][[sx]]@assays$RNA@counts
        monocle.raw.list[[sx]] <-
            new_cell_data_set(expression_data = mat.use,
                              cell_metadata = pbmc_list_raw[[experimentx]][[sx]]@meta.data)
    }
}

monocle.raw <- combine_cds(cds_list = monocle.raw.list, cell_names_unique = FALSE)
rowData(monocle.raw)$gene_short_name <- row.names(rowData(monocle.raw))

saveRDS(monocle.raw, file.path(dataDir, "scRNAseq_0_monocle_raw.cds"))

# Load Data ---------------------------------------------------------------
#remove comment when needed
monocle.obj <- readRDS(file.path(dataDir, "scRNAseq_0_souped.cds"))
monocle.raw <- readRDS(file.path(dataDir, "scRNAseq_0_monocle_raw.cds"))

# SoupX -------------------------------------------------------------------
sc <- SoupChannel(counts(monocle.raw), 
                  counts(monocle.obj))

sc <- setClusters(sc, setNames(monocle.obj@colData$Cluster, rownames(monocle.obj@colData)))
sc <- autoEstCont(sc)
out <- adjustCounts(sc)

doParallel::registerDoParallel(cores=7)
foreach(sx = unique(monocle.obj@colData$sample)) %dopar% {
x <- out[,grepl(sx, colnames(out))]
DropletUtils:::write10xCounts(file.path(dataDir, paste0("/countsSoupX/", sx)), x = x, type = "sparse", version = "3")
}

