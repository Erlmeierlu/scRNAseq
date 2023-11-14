library(monocle3)
library(SoupX)
library(dplyr)
library(foreach)

#setting up directories
baseDir <- getwd()
rawDir <- ("/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin/raw_data/scRNA_from_BSF/COUNT")
plotsDir <- file.path(baseDir, "plots/")
tablesDir <- file.path(baseDir, "tables/")
dataDir <- file.path(baseDir, "data/")


# Load Data ---------------------------------------------------------------
monocle.obj <- readRDS(file.path(dataDir, "old/scRNAseq_2a_monocle_with_fLN_40B3.cds"))
monocle.raw <- readRDS(file.path(dataDir, "old/scRNAseq_1_monocle_raw.cds"))


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

