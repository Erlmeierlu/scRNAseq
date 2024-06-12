library(tidyverse)
library(ROCit)
library(data.table)
library(monocle3)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)

# Directories and Used Functions -----------------------------------------------------------
source('functions/directories_and_theme.R')
source('functions/convert_genes.R')
source('functions/assign_ct_cluster.R')

#This function exports our object as anndata objects,
#which are used by the python package scanpy for scRNAseq
#analysis. 
export_to_anndata <- function(cds, filename){
  #first we translate the gene ids in human gene symbols
  rowData(cds)$human_symbol <-
    convert_gene_list(rowData(cds)$gene_short_name, #this is a self created function
                      convert_to = 'human')         #you find it in the script sourced above.
  #first convert it to seurat object
  seurat.obj <- as.Seurat(cds, data = NULL)
  #we normalize again
  seurat.obj <- NormalizeData(seurat.obj)
  
  #here we find out which columns are factors in our metadata..
  i <- sapply(seurat.obj@meta.data, is.factor)
  
  #...and we convert them to character columns. This is needed,
  #because of some weird things with an anndata object. 
  seurat.obj@meta.data[i] <- lapply(seurat.obj@meta.data[i], as.character)
  
  #first we have to save it as a h5Seurat object
  SaveH5Seurat(seurat.obj, file.path(vDir, 
                                     paste('3b',
                                           filename,
                                           'anndata.h5Seurat',
                                           sep = '_')),
               overwrite = TRUE)
  
  #then we convert it to a h5ad format, used by anndata
  Convert(file.path(vDir, 
                    paste('3b',
                          filename,
                          'anndata.h5Seurat',
                          sep = '_')), 
          dest = 'h5ad')
}

# Load Data ---------------------------------------------------------------
monocle.obj <- read_rds(file.path(vDir, '3a_first_annotation_full_dataset.cds'))
b_subset <- read_rds(file.path(vDir, '3a_first_annotation_b_subset.cds'))
t_subset <- read_rds(file.path(vDir, '3a_first_annotation_t_subset.cds'))
m_subset <- read_rds(file.path(vDir, '3a_first_annotation_m_subset.cds'))

# Gamma-Delta identification in T subset ----------------------------------

#Genes used for module calculation. Selection based on: 
#https://academic.oup.com/jleukbio/article/114/6/630/7223408#427715178
#
tgd_genes <- grep('^Tr[d][c|v]', rownames(t_subset), value = TRUE)
tab_genes <- grep('^Tr[a|b][c|v]', rownames(t_subset), value = TRUE)

#first convert it to seurat object
seurat.obj <- as.Seurat(t_subset, data = NULL)
seurat.obj <- NormalizeData(seurat.obj)
#calculate the alpha-beta or gamma-delta module scores 
seurat.obj <- AddModuleScore(seurat.obj, 
                             features = list(gd = tgd_genes,
                                             ab = tab_genes),
                             name = list('gamma.delta_module_score', 'alpha.beta_module_score')
)
#rename the columns
t_subset$gd.score <- seurat.obj@meta.data$gamma.delta_module_score1
t_subset$ab.score <- seurat.obj@meta.data$alpha.beta_module_score2

#assign either gamma-delta, or not-gd to cells, since this is all we tested
t_subset$ttype <- factor(case_when(t_subset$celltype == 'T_GD' ~ 'gamma-delta',
                                   t_subset$celltype == 'T_CD45.1' ~ 'not-gd',
                                   TRUE ~ NA))

#then we create an filtered object with only those that have a positive gd score
t_filtered <- t_subset[,t_subset$gd.score > 0 & !is.na(t_subset$ttype)]

#this is because we wanna calssify only based on those
ROCit_obj <- rocit(score = t_filtered$ab.score, class = t_filtered$ttype)

#The ROC curve tells us, about the predictive power of our model. In this 
#case whether we can classify well between GD and not-GD cells.
#This is in fact true. The Youden Index point (see plot) is the value
#that maximizes the true-positive-rate and minimizes the false-positive-rate.
plot(ROCit_obj)
#here we select the point that maximizes it
m <- which.max(ROCit_obj$TPR - ROCit_obj$FPR)
#then we select the value as our cutoff value for classification
t <- ROCit_obj$Cutoff[m]

#Finally, we can assign t cells as gamma delta based on that value
#Note that everything that is determined as gamma delta or CD45.1 by CITE-seq
#will not be changed, even though they might fall in a different category now.
t_subset$ttype <- factor(case_when(t_subset$celltype == 'T_GD' ~ 'gamma-delta',
                                   t_subset$celltype == 'T_CD45.1' ~ 'alpha-beta',
                                   t_subset$ab.score < t & t_subset$gd.score > 0 ~ 'gamma-delta',
                                   .default = 'alpha-beta'))

t_subset@colData %>%
    as_tibble %>% 
    ggplot(aes(ab.score, gd.score)) +
    geom_point(aes(col = ttype), size = 0.6, alpha = 0.3) +
    theme_my() + 
    facet_wrap(~celltype)+ 
    geom_vline(xintercept = ROCit_obj$Cutoff[m]) + 
    geom_hline(yintercept = 0)

ggsave(file.path(plotsDir, '3b_gamma_delta_classification.jpg'))

#Assignment of Gamma-Delta T cells based on umap cluster with enriched gamma
#delta signatures. CD45.1+ cells will be assigned by regions with bound CITE-
#seq ABs

coords <- reducedDims(t_subset)
coords <- data.frame(coords$UMAP)
colnames(coords) <- c("UMAP1", "UMAP2")
coords$ttype <- colData(t_subset)$ttype
coords$cd45x <- colData(t_subset)$cd45x

ggplot(coords) + 
  geom_hex(data = coords %>% dplyr::select(-ttype), 
           aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
  geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
  scale_fill_gradient(low = "deepskyblue2", high = "red", limits = c(0,40), oob = scales::squish)+
  facet_wrap( ~ttype, nrow = 2) +
  theme_my() 

ggsave(file.path(plotsDir, '3b_gamma_delta_classification_umap.jpg'))

ggplot(coords) + 
  geom_hex(data = coords %>% dplyr::select(-cd45x), 
           aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
  geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
  scale_fill_gradient(low = "deepskyblue2", high = "red", 
                      limits = c(0,40), 
                      oob = scales::squish) +
  facet_wrap(~cd45x, nrow = 2) +
  theme_my() 

ggsave(file.path(plotsDir, '3b_cd45_classification_umap.jpg'))

plot_cells(t_subset, group_label_size = 4)
ggsave(file.path(plotsDir, '3b_t_subset_cluster_umap'))

#visual inspection of the above plots results in following 
#preliminary assignments of whole clusters.
t_subset$celltype <- fcase(t_subset$Cluster %in% c(3, 15), 
                           'T_CD45.1',
                           t_subset$Cluster %in% c(16, 1, 4, 8, 9, 11, 5, 12, 24),
                           'T_GD',
                           default = 'T_cells')

t_subset$ct_cluster <- assign_ct_cluster(t_subset)

# Save T Subset -----------------------------------------------------------

write_rds(t_subset, file.path(vDir, '3b_first_annotation_t_subset.cds'))

# Export Subsets To AnnData for CellTypist -------------------------------------------------------

#Now we can export the subset as anndata objects in order to use them with python
#in the next script
walk2(list(t_subset, b_subset, m_subset), 
      c('t', 'b', 'm'), 
      export_to_anndata)
