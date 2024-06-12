library(monocle3)
library(magrittr)
library(readr)

# Directories and Used Functions -----------------------------------------------------------
source('functions/directories_and_theme.R')

# Load Data ---------------------------------------------------------------
monocle.obj <- read_rds(file.path(vDir, "1_qc.cds"))

# Unsupervised Analysis --------------------------------------------------------

##This section basically follows the monocle3 tutorial..

# Process dataset + PCA
monocle.obj <-
    preprocess_cds(monocle.obj, verbose = TRUE) %>%
    reduce_dimension(preprocess_method = "PCA", verbose = TRUE)

# batch correction. We correct for sample (i.e. information about sex, HDAC-
#group, organ, batch) and Phase. Note that this only has an effect on the
#UMAP-coordinates, not on the raw data (this will be corrected differently 
#in the DGE)!
monocle.obj <- align_cds(monocle.obj,
                         alignment_group = "sample", 
                         residual_model_formula_str = "~Phase",
                         verbose = TRUE
)

#Calculate UMAP coordinates
monocle.obj <- reduce_dimension(monocle.obj,
                                reduction_method = "UMAP",
                                preprocess_method = "Aligned",
                                verbose = TRUE)

# Clustering
monocle.obj <- cluster_cells(monocle.obj)

#Assigning cluster information to meta data for ease of use
monocle.obj@colData$Cluster <- unname(clusters(monocle.obj[,rownames(colData(monocle.obj))]))

# Export Data for later use -----------------------------------------------
# Store dataset
write_rds(monocle.obj, file.path(vDir, "2_unsupervised.cds"))

