library(dplyr)
library(stringr)
library(monocle3)
library(Seurat)

#setting up directories
baseDir <- getwd()
plotsDir <- file.path(baseDir, "plots/")
tablesDir <- file.path(baseDir, "tables/")
dataDir <- file.path(baseDir, "data/")

# Create Data ---------------------------------------------------------------

#directories to loop over
fileDirs <- list.files(file.path(dataDir, "countsSoupX"))

#Seurat Object
pbmc_data_list <- list()
pbmc_list <- list()

for(file in fileDirs){
    experiment <- case_when(!grepl("[0-9]+B[0-9]", file) ~ "HDAC1",
                            grepl("40", file) ~ "HDAC1",
                            TRUE ~ "HDAC2")
    
    pbmc_data_list[[experiment]][[file]] <- Read10X(file.path(dataDir, "countsSoupX", file))
    
    pbmc_list[[experiment]][[file]] <-
        CreateSeuratObject(
            counts = pbmc_data_list[[experiment]][[file]],
            project = "scRNAseq",
            min.cells = 0,
            min.features = 0
        )
}

CountsAB <- readRDS(file.path(dataDir, "ab_counts.rds"))


# Assign Metadata ---------------------------------------------------------

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
               subset = nFeature_RNA > 200 &
                   percent.mt < 10 & hash.sum > 10)
}

#Create Monocle Ob
monocle.obj.list <- list()
seurat.obj <- list()
for (experimentx in names(pbmc_subset)) {
    for (sx in names(pbmc_subset[[experimentx]])) {
        seurat.obj[[sx]] <- pbmc_subset[[experimentx]][[sx]]
        mat.use <- seurat.obj[[sx]]@assays$RNA@counts
        monocle.obj.list[[sx]] <-
            new_cell_data_set(expression_data = mat.use,
                              cell_metadata = seurat.obj[[sx]]@meta.data)
    }
}
monocle.obj <- combine_cds(cds_list = monocle.obj.list, cell_names_unique = TRUE)


# Metadata refinement
monocle.obj@colData$hash <- as.factor(monocle.obj@colData$hash)
monocle.obj@colData$hash.agg <- as.factor(monocle.obj@colData$hash.agg)
monocle.obj@colData$gdcd4 <- as.factor(monocle.obj@colData$gdcd4)
monocle.obj@colData$cd4cd8 <- as.factor(monocle.obj@colData$cd4cd8)
monocle.obj@colData$cd45x <- as.factor(monocle.obj@colData$cd45x)
monocle.obj@colData$sample <- as.factor(monocle.obj@colData$sample)
monocle.obj@colData$organ <- as.factor(str_extract(monocle.obj@colData$sample, "[A-Z][:alpha:]+"))
levels1 <- monocle.obj@colData %>% as_tibble %>% filter(organ == "Skin") %>% pull(sample) %>% droplevels() %>%  levels
levels2 <- monocle.obj@colData %>% as_tibble %>% filter(organ == "LN") %>% pull(sample) %>% droplevels() %>%  levels
levelsx <- paste(c(levels1, levels2))
monocle.obj@colData$sample <- monocle.obj@colData$sample %>% forcats::fct_relevel(levelsx)
monocle.obj@colData$Phase <- as.factor(monocle.obj@colData$Phase)
monocle.obj@colData$batch <-
    as.factor(str_extract(monocle.obj@colData$sample, "[A-Z][0-9]$"))
monocle.obj@colData$treatment <-
    as.factor(case_when(
        grepl("30[1-3]$", monocle.obj@colData$hash) ~ "HDAC_WT",
        grepl("30[4-6]$", monocle.obj@colData$hash) ~ "HDAC_cKO",
        grepl("30[7-9]$", monocle.obj@colData$hash) ~ "NoT",
        TRUE ~ "Undefined"
    ))
monocle.obj@colData$treatment.agg <-
    as.factor(case_when(
        grepl("\\-3$", monocle.obj@colData$hash.agg) ~ "NoT",
        grepl("\\-1$", monocle.obj@colData$hash.agg) ~ "HDAC_WT",
        grepl("\\-2$", monocle.obj@colData$hash.agg) ~ "HDAC_cKO",
        TRUE ~ "Undefined"
    ))
monocle.obj@colData$hashid <-
    as.factor(ifelse(
        monocle.obj@colData$hash == "Undefined",
        "Undefined",
        str_extract(monocle.obj@colData$hash, "[0-9]+$")
    ))
monocle.obj@colData$hashid.agg <-
    as.factor(ifelse(
        monocle.obj@colData$hash.agg == "Undefined",
        "Undefined",
        str_extract(monocle.obj@colData$hash.agg, "[0-9]+$")
    ))
monocle.obj@colData$sex <- as.factor(str_to_upper(str_extract(monocle.obj@colData$sample,"^[A-z]")))
monocle.obj@colData$treatment <-
    monocle.obj@colData$treatment %>% forcats::fct_relevel("NoT", "HDAC_WT", "HDAC_cKO", "Undefined")
monocle.obj@colData$treatment.agg <-
    monocle.obj@colData$treatment.agg %>% forcats::fct_relevel("NoT", "HDAC_WT", "HDAC_cKO", "Undefined")
monocle.obj@colData$sample_treat <- as.factor(paste(monocle.obj@colData$sample, monocle.obj@colData$treatment, sep = "___"))
monocle.obj@colData$sample_treat.agg <-
    as.factor(paste(
        monocle.obj@colData$sample,
        monocle.obj@colData$treatment.agg,
        sep = "___"
    ))
monocle.obj@colData$experiment <-
    factor(case_when(
        !grepl("[0-9]+B[0-9]", monocle.obj@colData$sample) ~ "HDAC1",
        grepl("40", monocle.obj@colData$sample) ~ "HDAC1",
        TRUE ~ "HDAC2"
    ))

rowData(monocle.obj)$gene_short_name <- row.names(rowData(monocle.obj))


#remove bad sample
monocle.obj <- monocle.obj[, !grepl("fLN_40B3", colnames(monocle.obj))]
monocle.obj@colData <- droplevels(monocle.obj@colData)

nrow(colData(monocle.obj)) == ncol(monocle.obj)

#save data
saveRDS(monocle.obj, file.path(dataDir, "scRNAseq_0_souped.cds"))
