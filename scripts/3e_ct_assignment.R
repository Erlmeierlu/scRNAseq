#This is after looking at clusters and ct annotations in 3d_ct_assignment.R script
library(tidyverse)
library(data.table)
library(monocle3)
library(Seurat)
library(SeuratDisk)

# Directories and Used Functions -----------------------------------------------------------
source('functions/directories_and_theme.R')
source('functions/explore_annotation.R')
source('functions/assign_ct_cluster.R')
source('functions/convert_genes.R')

#Unique Functions
clean_treatment <- function(cds){
    remove_t <- function(col){
        case_when(cds[[col]] == 'NoT' & grepl('CD45.1', cds$celltype) ~
                      'Undefined',
                  TRUE ~ cds[[col]]) %>% 
            factor(levels = c('NoT',
                              'HDAC_WT',
                              'HDAC_cKO',
                              'Undefined'))
    }
    cds$treatment <- remove_t('treatment')
    cds$treatment.agg <- remove_t('treatment.agg')
    cds
}

export_to_anndata <- function(cds, filename, remove_pre_file = TRUE){
    rowData(cds)$human_symbol <-
        convert_gene_list(rowData(cds)$gene_short_name,
                          convert_to = 'human')
    
    seurat.obj <- as.Seurat(cds, data = NULL)
    seurat.obj <- NormalizeData(seurat.obj)
    
    i <- sapply(seurat.obj@meta.data, is.factor)
    
    seurat.obj@meta.data[i] <- lapply(seurat.obj@meta.data[i], as.character)
    
    file <- file.path(dataDir, 
                      paste('3e',
                            filename,
                            'anndata.h5Seurat',
                            sep = '_'))
                      
    SaveH5Seurat(seurat.obj, file, overwrite = TRUE)
    
    Convert(file, 
            dest = 'h5ad',
            overwrite = TRUE)
    if(isFALSE(remove_pre_file)) return(NULL)
    file.remove(file)
}

# Load Data ---------------------------------------------------------------
monocle.obj <- read_rds(file.path(vDir, '3a_first_annotation_full_dataset.cds'))
b_subset <- read_rds(file.path(vDir, '3a_first_annotation_b_subset.cds'))
t_subset <- read_rds(file.path(vDir, '3b_first_annotation_t_subset.cds'))
m_subset <- read_rds(file.path(vDir, '3a_first_annotation_m_subset.cds'))

# Refine Assignment -------------------------------------------------------
## T Subset ----------------------------------------------------------------
#This assignment is now refined by looking at the first annotation plots.
#Double check ?
t_subset$celltype <- case_when(t_subset$Cluster %in% c(15, 3) ~ 'T_CD45.1_1',
                               t_subset$Cluster == 6 ~ 'T_CD45.1_2',
                               t_subset$Cluster == 2 ~ 'T_CD45.1_3',
                               t_subset$Cluster %in% c(1, 7) ~ 'T_cell_2',
                               t_subset$Cluster %in% c(8, 9, 11) ~ 'T_gd_2',
                               t_subset$Cluster %in% c(5, 12) ~ 'T_gd_1',
                               t_subset$Cluster == 16 ~ 'T_gd_4',
                               t_subset$Cluster == 18 ~ 'T_gd_5',
                               t_subset$Cluster == 4 ~ 'T_gd_3',
                               t_subset$Cluster %in% c(13, 17) ~ 'T_cell_1',
                               .default = 'Undefined_T')

t_subset$ct_cluster <- assign_ct_cluster(t_subset)

## B Subset -----------------------------------------------
#trying to define larger clusters and identify differences between them
b_subset$cluster_big <-
    case_when(
        b_subset$Cluster %in% c(3, 29, 33, 28, 35, 32) ~ 'left',
        b_subset$Cluster %in% c(15, 21, 43, 49) ~ 'bottom',
        b_subset$Cluster %in% c(14, 31, 24, 6, 48) ~ 'right',
        b_subset$Cluster %in% c(4, 30) ~ 'bottom_right',
        b_subset$Cluster %in% c(45, 46, 50, 51, 52, 53, 54, 55) ~ 'small',
        TRUE ~ 'B_follicular'
    )

markers <- top_markers(b_subset, group_cells_by = 'cluster_big', 
                       reference_cells = 20000)

top_markers <- markers %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(5, pseudo_R2)

marker_genes <- unique(top_markers %>% pull(gene_id))

plot_genes_by_group(b_subset, 
                    group_cells_by = 'cluster_big', 
                    markers = marker_genes, 
                    max.size = 3, 
                    ordering_type = 'max')
ggsave(file.path(plotsDir, '3e_b_subset_top_markers_large_clusters.pdf'))

plot_cells(b_subset, 
           color_cells_by = 'cluster_big',
           group_label_size = 3,
           label_groups_by_cluster = F)
ggsave(file.path(plotsDir, '3e_b_subset_large_cluster_umap.jpg'))

#Assignment based on some differences in top marker genes
#Please note that this is VERY simplified
b_subset$celltype <- case_when(b_subset$cluster_big == 'left' ~ 'B_follicular_Fth_neg',
                               b_subset$cluster_big == 'small' ~ 'Unassigned_B',
                               b_subset$cluster_big == 'bottom' ~ 'B_follicular_Mif+',
                               b_subset$cluster_big == 'right' ~ 'B_naive_mem',
                               .default = 'B_follicular')

b_subset$cluster_big <- NULL
b_subset$ct_cluster <- assign_ct_cluster(b_subset)

## M Subset ----------------------------------------------------------------
#This assignment is now refined by looking at the first annotation plots.
#Double check ?
m_subset$celltype <- case_when(m_subset$Cluster %in% c(1,2,3,4,7,11) ~
                                   'Macrophages_1',
                               m_subset$Cluster %in% c(5,9,16) ~
                                   'Macrophages_2',
                               m_subset$Cluster %in% c(6, 8, 10, 12) ~
                                   'Macrophages_3',
                               m_subset$Cluster %in% c(15, 13, 18, 20) ~
                                   'Undefined_M',
                               m_subset$Cluster == 17 ~
                                   'Plasma_cells',
                               m_subset$Cluster == 19 ~
                                   'Neutrophils',
                               m_subset$Cluster == 14 ~
                                   'Monocytes')

markers <- top_markers(m_subset, group_cells_by = 'celltype', 
                       reference_cells = 20000)

top_markers <- markers %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(5, pseudo_R2)

marker_genes <- unique(top_markers %>% pull(gene_id))

plot_genes_by_group(m_subset, 
                    group_cells_by = 'celltype', 
                    markers = marker_genes, 
                    max.size = 3, 
                    ordering_type = 'max')
ggsave(file.path(plotsDir, '3e_m_subset_top_markers.pdf'))

plot_cells(m_subset, 
           color_cells_by = 'celltype',
           group_label_size = 3,
           label_groups_by_cluster = F)
ggsave(file.path(plotsDir, '3e_m_subset_celltype_umap.jpg'))

m_subset$ct_cluster <- assign_ct_cluster(m_subset)

## Struc Subset --------------------------

#Structural Cells got labelled as B cells somehow

# plot_cells(monocle.obj,
#            color_cells_by = 'celltype',
#            group_label_size = 3,
#            label_groups_by_cluster = F
#            )
# kc_b_cluster <- c(5, 8, 18)

# res <- load_singler(file.path(tablesDir, 'cell_types_all_ref.csv'), 
#                     subset = monocle.obj[,monocle.obj$Cluster %in% kc_b_cluster])
# res <- prepare_res(ds = res, 
#                    subset = monocle.obj[,monocle.obj$Cluster %in% kc_b_cluster])
# 
# res[res[, by = .(label, cluster), frank(freq)]$V1 < 2]
# 
# #dont know how those got not assigned as KCs.. all seem to be KC. Changed Label
# 
# monocle.obj$celltype <- case_when(monocle.obj$Cluster %in% kc_b_cluster ~ 'Keratinocytes',
#                                   TRUE ~ monocle.obj$celltype)

subset <- monocle.obj[, monocle.obj$celltype %in% c('Fibroblasts', 'Keratinocytes')]
subset@colData <- subset@colData %>% droplevels()

subset$ct_cluster <- assign_ct_cluster(subset)

kc_subset <- subset[,subset$celltype == 'Keratinocytes']
markers <- top_markers(kc_subset, group_cells_by = 'ct_cluster', 
                       reference_cells = 20000)

top_markers <- markers %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(5, pseudo_R2)

marker_genes <- unique(top_markers %>% pull(gene_id))

plot_genes_by_group(kc_subset, 
                    group_cells_by = 'ct_cluster', 
                    markers = marker_genes,
                    max.size = 3,
                    ordering_type = 'max')
ggsave(file.path(plotsDir, '3e_kc_top_markers.pdf'), height = 7)

fb_subset <- subset[, subset$celltype == 'Fibroblasts']
markers <- top_markers(fb_subset, group_cells_by = 'ct_cluster', 
                       reference_cells = 20000)

top_markers <- markers %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(5, pseudo_R2)

marker_genes <- unique(top_markers %>% pull(gene_id))

plot_genes_by_group(fb_subset, 
                    group_cells_by = 'ct_cluster', 
                    markers = marker_genes,
                    max.size = 3, 
                    ordering_type = 'max')
ggsave(file.path(plotsDir, '3e_fb_top_markers.pdf'))

subset$celltype <- as_tibble(subset@colData) %>%
    group_by(celltype) %>%
    mutate(celltype = paste(celltype, dense_rank(Cluster), sep = '_')) %>% 
    pull(celltype)

# Assign to Full Object ---------------------------------------------------

plot_cells(monocle.obj)
ggsave(file.path(plotsDir, '3e_full_dataset_cluster_umap.jpg'))

monocle.obj$celltype[match(colnames(t_subset), 
                           colnames(monocle.obj))] <- t_subset$celltype
monocle.obj$celltype[match(colnames(b_subset), 
                           colnames(monocle.obj))] <- b_subset$celltype
monocle.obj$celltype[match(colnames(m_subset), 
                           colnames(monocle.obj))] <- m_subset$celltype
monocle.obj$celltype[match(colnames(subset), 
                           colnames(monocle.obj))] <- subset$celltype

#Cluster 93 looks might be melanocytes based on expression of
monocle.obj$celltype <- case_when(monocle.obj$Cluster == 93 ~ 'Melanocytes',
                                  TRUE ~ monocle.obj$celltype)

#Raw Cell Type to avoid many cts in full object
#Basically every cell, that is in one of the subsets is named after the subset.
#Undefined_X are cells that were annotated as a certain cell type of a subset, 
#but are not included in one of the subsets (the subsets were selected based
#on position on UMAP). 
#
monocle.obj$raw_celltype <- case_when(colnames(monocle.obj) %in% colnames(t_subset) ~
                                          'T_cells',
                                      colnames(monocle.obj) %in% colnames(b_subset) ~
                                          'B_cells',
                                      colnames(monocle.obj) %in% colnames(m_subset) ~
                                          'Myeloid_cells',
                                      grepl('^K', monocle.obj$celltype) ~ 'Keratinocytes',
                                      grepl('^F', monocle.obj$celltype) ~ 'Fibroblasts',
                                      grepl('^T', monocle.obj$celltype) ~ 'Undefined_T',
                                      grepl('^B', monocle.obj$celltype) ~ 'Undefined_B',
                                      monocle.obj$celltype == 'DC' ~ 'Undefined',
                                      TRUE ~ monocle.obj$celltype)

plot_cells(monocle.obj, color_cells_by = 'raw_celltype')

monocle.obj$celltype <- case_when(grepl('^removed|Undefined', monocle.obj$raw_celltype) ~ 
                                      monocle.obj$raw_celltype, 
                                  TRUE ~ monocle.obj$celltype)


#also for ct_cluster
monocle.obj$ct_cluster <- assign_ct_cluster(monocle.obj)

# backup for monocle object -----------------------------------------------
#because cd45.1 that are in NoT will be removed here.. 
write_rds(monocle.obj, file.path(vDir, '3_annotated_monocle_backup.cds'))
write_rds(t_subset, file.path(vDir, '3_t_subset_backup.cds'))
write_rds(b_subset, file.path(vDir, '3_b_subset_backup.cds'))
write_rds(m_subset, file.path(vDir, '3_m_subset_backup.cds'))

monocle.obj <- clean_treatment(monocle.obj)

t_subset <- clean_treatment(t_subset)

# Export to AnnData -------------------------------------------------------
#Full Dataset
export_to_anndata(monocle.obj, 'full')
# save monocle object ------------------------------------------------

write_rds(monocle.obj, file.path(dataDir, '3_annotated_monocle.cds'))
write_rds(t_subset, file.path(dataDir, '3_t_subset.cds'))
write_rds(b_subset, file.path(dataDir, '3_b_subset.cds'))
write_rds(m_subset, file.path(dataDir, '3_m_subset.cds'))

