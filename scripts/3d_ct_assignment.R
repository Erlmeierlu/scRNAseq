#For this script it is required that script 3c is run on python..
library(tidyverse)
library(data.table)
library(patchwork)
library(monocle3)
library(readxl)

# Directories and Used Functions -----------------------------------------------------------
source('functions/directories_and_theme.R')
#Check the file for detailed functions
source('functions/explore_annotation.R')

# Load Data ---------------------------------------------------------------
monocle.obj <- read_rds(file.path(vDir, '3a_first_annotation_full_dataset.cds'))
b_subset <- read_rds(file.path(vDir, '3a_first_annotation_b_subset.cds'))
t_subset <- read_rds(file.path(vDir, '3b_first_annotation_t_subset.cds'))
m_subset <- read_rds(file.path(vDir, '3a_first_annotation_m_subset.cds'))

#T Sub Type Markers
#These reference genes were collected by Melanie
new_t_markers <- read_xlsx(file.path(tablesDir, '3d_MarkerGenes_Thelper_ML.xlsx'))

# Exploring Subset Clusters ------------------------------------------------------

## T subset ----------------------------------------------------------------
#Most functions here can be found in the script mentioned at the start
#check them out for more details

#first load celltypist annotations
t_celltypist <- load_celltypist(file.path(tablesDir, 
                                          '3c_lymphoid_res_celltypist.csv'))

#then we load SingleR results
t_singler <- load_singler(file.path(tablesDir, 
                                    c('3a_SingleR_res_t_subset.csv',
                                      '3a_SingleR_res_full_dataset.csv')),
                          t_subset)

t_res <- prepare_res(t_singler, t_celltypist, t_subset)

t_res <- assign_factor_levels(t_res, order_by = 'CellTypist_label.fine')
fc_levels <- t_res$ct_cluster %>% levels

#Automated Annotation
p1 <- plot_automated_anno(t_res, ignore_ref = c('immgen_label.main',
                                                'struc_label.fine',
                                                'th_rest_label.fine'))

#CITE-seq distribution
t_abdat <- generate_abdata(t_subset, 
                           fc_levels, 
                           abgroups = c('cd4cd8', 'cd45x'))
p2 <- plot_abs(t_abdat)

#Organ
t_org <- prepare_organ_data(t_subset, fc_levels)
p3 <- plot_org(t_org)

#Top Marker genes
marker_genes <- get_top_markers(t_subset)

p4 <- plot_genes(t_subset,
                 cell_factor_levels = fc_levels,
                 markers = setNames(list(marker_genes), 'top_markers'),
                 ordering_type = 'maximal_on_diag'
) + theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none') +
  ggtitle('Top Marker Genes')

new_t_markers <- new_t_markers %>% lapply(na.omit)

names(new_t_markers) <- names(new_t_markers) %>% 
  str_remove_all(' marker genes| cells')

p5 <- plot_genes(t_subset,
                 cell_factor_levels = fc_levels,
                 markers = new_t_markers) + theme(axis.text.y = element_blank(),
                                                  axis.ticks.y = element_blank(),
                                                  axis.title.y = element_blank()) +
    ggtitle('T subtype marker gene expression') 

#Patchwork plots together
p1+p2+p3+p4+p5+
  plot_layout(widths = c(5.5, 2, 0.5, 6, 14.5))+
  plot_annotation(tag_levels = 'A')

ggsave(file.path(plotsDir, '3d_t_first_annotation_plot.jpg'), 
       height = 10, 
       width = 38)

plot_cells(t_subset, group_label_size = 4)
ggsave(file.path(plotsDir, '3d_umap_t_subset_cluster.jpg'))
## B subset ----------------------------------------------------------------
b_celltypist <- load_celltypist(file.path(tablesDir, '3c_b_res_celltypist.csv'))

b_singler <- load_singler(c(file.path(tablesDir, '3a_SingleR_res_b_subset.csv'),
                            file.path(tablesDir, '3a_SingleR_res_full_dataset.csv')),
                          b_subset)
b_res <- prepare_res(b_singler, b_celltypist, b_subset)
b_res <- assign_factor_levels(b_res, order_by = 'CellTypist_label.fine')

fc_levels <- b_res$ct_cluster %>% levels

p1 <- plot_automated_anno(b_res, use_ref = c('CellTypist_label.fine',
                                             'b_data_label.fine',
                                             'immgen_label.fine',
                                             'struc_label.fine'))
#Organ
b_org <- prepare_organ_data(b_subset, fc_levels)
p2 <- plot_org(b_org)

#Top Marker genes
marker_genes <- get_top_markers(b_subset)

p3 <- plot_genes(b_subset,
                 cell_factor_levels = fc_levels,
                 markers = setNames(list(marker_genes), 'top_markers'),
                 ordering_type = 'maximal_on_diag'
) + theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none') +
  ggtitle('Top Marker Genes')


p1+p2+p3+
  plot_layout(widths = c(13,1.5,12))+
  plot_annotation(tag_levels = 'A')
ggsave(file.path(plotsDir, '3d_b_first_annotation_plot.jpg'), height = 12, width = 21)

plot_cells(b_subset, group_label_size = 4)
ggsave(file.path(plotsDir, '3d_umap_b_subset_cluster.jpg'))

## M subset ----------------------------------------------------------------
m_celltypist <- load_celltypist(file.path(tablesDir, '3c_myeloid_res_celltypist.csv'))

m_singler <- load_singler(file.path(tablesDir, '3a_SingleR_res_full_dataset.csv'),
                          m_subset)
m_res <- prepare_res(m_singler, m_celltypist, m_subset)
m_res <- assign_factor_levels(m_res, order_by = 'CellTypist_label.fine')

fc_levels <- m_res$ct_cluster %>% levels

p1 <- plot_automated_anno(m_res,
                          use_ref = c('CellTypist_label.fine',
                                      'immgen_label.fine',
                                      'struc_label.fine'))
#Organ
m_org <- prepare_organ_data(m_subset, fc_levels)
p2 <- plot_org(m_org)

#Top Marker genes
marker_genes <- get_top_markers(m_subset)

p3 <- plot_genes(m_subset,
                 cell_factor_levels = fc_levels,
                 markers = setNames(list(marker_genes), 'top_markers'),
                 ordering_type = 'maximal_on_diag'
) + theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none') +
  ggtitle('Top Marker Genes')


p1+p2+p3+
  plot_layout(widths = c(11,1.5,10))+
  plot_annotation(tag_levels = 'A')
ggsave(file.path(plotsDir, '3d_m_first_annotation_plot.jpg'), height = 8, width = 21)

plot_cells(m_subset, group_label_size = 4)
ggsave(file.path(plotsDir, '3d_umap_m_subset_cluster.jpg'))

## FB/KC subset ----------------------------------------------------------------
s_subset <- monocle.obj[,grepl('^Fib|^Ker', 
                               monocle.obj$celltype)]
s_subset@colData <- s_subset@colData %>% droplevels

s_singler <- load_singler(file.path(tablesDir, '3a_SingleR_res_full_dataset.csv'),
                          s_subset)

s_res <- prepare_res(s_singler, subset = s_subset)
s_res <- assign_factor_levels(s_res, order_by = 'immgen_label.fine')

fc_levels <- s_res$ct_cluster %>% levels

#Automated Annotation
p1 <- plot_automated_anno(s_res,
                          use_ref = c('immgen_label.fine',
                                      'old_label.fine'))

#CITE-seq Abs
s_abdat <- generate_abdata(s_subset, levels = fc_levels, abgroups = 'hashid')
p2 <- plot_abs(s_abdat)

#Organ
s_org <- prepare_organ_data(s_subset, fc_levels)
p3 <- plot_org(s_org)

#Top Marker genes
marker_genes <- get_top_markers(s_subset)

p4 <- plot_genes(s_subset,
                 cell_factor_levels = fc_levels,
                 markers = setNames(list(marker_genes), 'top_markers'),
                 ordering_type = 'maximal_on_diag'
) + theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none') +
  ggtitle('Top Marker Genes')


p1+p2+p3+p4+
  plot_layout(widths = c(11,3, 1.5,10))+
  plot_annotation(tag_levels = 'A')
ggsave(file.path(plotsDir, '3d_structural_first_annotation_plot.jpg'), 
       height = 8, width = 21)

plot_cells(s_subset, group_label_size = 4)
ggsave(file.path(plotsDir, '3d_umap_structural_subset_cluster.jpg'))
