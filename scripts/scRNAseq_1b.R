library(monocle3)
library(magrittr)

#setting up directories
baseDir <- getwd()
plotsDir <- file.path(baseDir, "plots/")
tablesDir <- file.path(baseDir, "tables/")
dataDir <- file.path(baseDir, "data/")

#seeeed
set.seed(42)

# Load Data ---------------------------------------------------------------
monocle.obj <- readRDS(file.path(dataDir, "scRNAseq_0_souped.cds"))

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
saveRDS(monocle.obj, file = file.path(dataDir, "scRNAseq_1_monocle.cds"))
