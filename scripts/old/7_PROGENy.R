library(monocle3)
library(decoupleR)
library(readr)
library(Seurat)
library(SeuratObject)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#setting up directories
gfsDir <- '/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin'
plotsDir <- file.path(gfsDir, 'plots')
tablesDir <- file.path(gfsDir, 'tables')
oldDir <- "/vscratch/scRNAseq/data/old"
shinyDir <- 'dge-app'
dataDir <-"data"
resDir <- "results"
# Load Data ---------------------------------------------------------------

monocle.obj <- read_rds(file.path(dataDir, '3_annotated_monocle.cds'))
t_subset <- read_rds(file.path(dataDir, '3_t_subset.cds'))

# Analysis ----------------------------------------------------------------

seurat.obj <- SeuratObject::as.Seurat(t_subset, data = NULL)
seurat.obj <- NormalizeData(seurat.obj)

net <- get_progeny(organism = 'mouse', top = 500)
mat <- as.matrix(seurat.obj@assays$originalexp@data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

seurat.obj[['pathwaysmlm']] <- acts %>%
    pivot_wider(id_cols = 'source', names_from = 'condition',
                values_from = 'score') %>%
    column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = seurat.obj) <- "pathwaysmlm"

# Scale the data
seurat.obj <- ScaleData(seurat.obj)
seurat.obj@assays$pathwaysmlm@data <- seurat.obj@assays$pathwaysmlm@scale.data


p1 <- DimPlot(seurat.obj, group.by = 'organ', reduction = "UMAP", label = TRUE, pt.size = 0.5) + 
    NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(seurat.obj, features = c("JAK-STAT")) & 
           scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
    ggtitle('JAK-STAT activity')

p1+p2
