library(monocle3)
library(magrittr)
library(readr)

#setting up directories
vDir <- ("/vscratch/scRNAseq")
plotsDir <- ("/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin/plots")
tablesDir <- file.path(vDir, "tables")
oldDir <- file.path(vDir, "data/old")
dataDir <-("data")
resDir <- ("results")

# Load Data ---------------------------------------------------------------
monocle.obj <- read_rds(file.path(dataDir, "1_qc.cds"))

# Unsupervised Analysis --------------------------------------------------------
# Process dataset
monocle.obj <-
    preprocess_cds(monocle.obj, verbose = TRUE) %>%
    reduce_dimension(preprocess_method = "PCA", verbose = TRUE)

# batch correction
monocle.obj <- align_cds(monocle.obj,
                         alignment_group = "sample", 
                         residual_model_formula_str = "~Phase",
                         verbose = TRUE
)

#UMAP
monocle.obj <- reduce_dimension(monocle.obj,
                                reduction_method = "UMAP",
                                preprocess_method = "Aligned",
                                verbose = TRUE)


# Clustering
monocle.obj <- cluster_cells(monocle.obj)

monocle.obj@colData$Cluster <- unname(clusters(monocle.obj[,rownames(colData(monocle.obj))]))

# Store dataset
write_rds(monocle.obj, file.path(dataDir, "2_unsupervised.cds"))
