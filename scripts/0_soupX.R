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
library(DropletUtils)

# Directories and Used Functions -----------------------------------------------------------
source('functions/directories_and_theme.R')
source('functions/load_raw_data.R')

#Unique Functions
#This function provides the first preprocessing steps within Seurat
#This is necessary, because the DoubletFinder package works with Seurat
#only.
seurat_analysis <- function(object){
    #This function can be used to determine the dimensionality of the data.
    #Meaning, it calculates how many principle components (PCs) are needed
    #to capture the majority of the variation in the data.
    DetermineDimensionality <- function(object){
        #First Part:
        #Determine percent of variation associated with each PC
        pct <- object[["pca"]]@stdev / sum(object[["pca"]]@stdev) * 100
        #Calculate cumulative percents for each PC
        cumu <- cumsum(pct)
        # Determine which PC has cumulative percent greater than 90% and 
        # % variation associated with the PC is less than 5
        co1 <- which(cumu > 90 & pct < 5)[1]
        
        #Second Part:
        #Determine the difference between variation of PC and subsequent PC
        #and find last point where change of % of variation is more than 0.1%
        co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
                    decreasing = T)[1] + 1
        
        #Minimum of the two calculation tells us how many PCs are needed to 
        #capture majority of the variation
        pcs <- min(co1, co2)
        pcs
    }
    
    #Standard Workflow
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object)
    
    #Determine Dimensionality + Unsup Ana
    dims <- seq(1, DetermineDimensionality(object))
    object <- FindNeighbors(object, dims = dims)
    object <- FindClusters(object)
    object <- RunUMAP(object, dims = dims)
    
    #Doublet Identification
    #This is part of the standard workflow from the DoubletFinder package
    sweep.list <- paramSweep_v3(object, PCs = dims, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pk <- bcmvn[which.max(bcmvn$BCmetric), "pK"] %>% 
        droplevels() %>% 
        levels %>% 
        as.numeric()
    
    #Next we set the number of expected doublets in each sample..
    #Unfortunately this is arbitrary, and can lead to faulty results.
    #Here, 7.5% doublets are expected.
    
    #The dev talks about how to best determine this number. I chose the 
    #standard value here. However, there might be a better value. Check
    #this link, and rerun the analysis. However, it might lead to some 
    #changes down the line...
    #https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76#issuecomment-624283980
    
    nExp <- round(0.075*nrow(object@meta.data))
    object <- doubletFinder_v3(object, PCs = dims, pK = pk, nExp = nExp)
    
    #just some reorganizing of the seurat object meta data..
    colnames(object@meta.data)[grepl(
        "DF.class", 
        colnames(object@meta.data))] <- "doublet_classification"
    
    colnames(object@meta.data)[grepl(
        "pANN", 
        colnames(object@meta.data))] <- "pANN"
    
    #and return it
    object
}


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


#Excluded the sample because we think somethign was mixed up here(wrong
#hashtag?)

fileDirs <- lapply(fileDirs, function(x) {
    x[!grepl("fLN_40B3", x)]
})


# Load Raw Data -----------------------------------------------------------

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

#Save it, cause we will need it later aswell.
write_rds(CountsAB, file.path(dataDir, "0_ab_counts.rds"))

# Pre-Processing ---------------------------------------------------------

#exclude barcodes with less than 50 UMIs for SoupX analysis (later)
data_subset <-
    lapply(data_list,
           subset,
           subset = nFeature_RNA > 50)

#here we use the function defined at the top to make the seurat+doubletfinder
#steps
seu <- lapply(data_subset, seurat_analysis)

#Integrate into one large dataset
seu <- Merge_Seurat_List(seu, add.cell.ids = names(seu))

#Renaming the barcodes so it fits with the later analysis with monocle
new_names <-
    colnames(seu) %>%
    str_c(str_extract(., "^[:alpha:]+_[:alnum:]+") , sep = "_") %>%
    str_remove("^[:alpha:]+_[:alnum:]+_")
seu <- RenameCells(seu, new.names = new_names)

# Create Objects --------------------------------------------------
#Save Seurat object if you want to. Not needed in future scripts. Just
#Here for a potential checkpoint.
# write_rds(seu, file.path(vDir, '0_souped.rds'))
# seu <- read_rds(file.path(vDir, '0_souped.rds'))

#Create Lookup-Table for later use in QC
doublet_LUT <- seu@meta.data[, "doublet_classification", drop = F]
doublet_LUT <- as.data.table(doublet_LUT, keep.rownames = T) 

write_rds(doublet_LUT, file.path(dataDir, "0_doublet_LUT.rds"))

#Now we load the raw lists used for SoupX analysis. 10X genomics provides
#and raw version of the data, and an filtered version. The filtered version
#is usually used for the analysis. SoupX utilizes the raw version to
#determine ambient RNA present in each droplet..

#Raw reads as monocle object
monocle.raw.list <- lapply(data_list_raw, function(x){
    mat <- x@assays$RNA@counts
    cds <- new_cell_data_set(expression_data = mat,
                             cell_metadata = x@meta.data)
    cds
})


#integrate in one data set
monocle.raw <- combine_cds(cds_list = monocle.raw.list, cell_names_unique = FALSE)
rowData(monocle.raw)$gene_short_name <- row.names(rowData(monocle.raw))

# write_rds(monocle.raw, file.path(vDir, "0_monocle_raw.cds"))

# Load Data ---------------------------------------------------------------
#remove comment when needed
# seu <- read_rds(file.path(vDir, '0_souped.rds'))
# monocle.raw <- read_rds(file.path(vDir, "0_monocle_raw.cds"))

#time to free some memory..
data_list_raw <- data_list <- data_subset <- monocle.raw.list <- NULL
gc() #garbage collector frees memory

# SoupX -------------------------------------------------------------------

sample <- str_extract(
    colnames(seu), "[:alpha:]+_[:alnum:]+$"
    ) %>% unique()

# parallel computing for faster speed
doParallel::registerDoParallel(cores=6)
foreach(sx = sample) %dopar% {
    
    #standard soupx stuff. Removal of ambient RNA
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
    
    #here we write the corrected counts to the disk. We will use these in
    #the analysis
    DropletUtils:::write10xCounts(file.path(vDir, paste0("countsSoupX/", sx)), 
                                  x = out, 
                                  type = "sparse", 
                                  version = "3")
}
