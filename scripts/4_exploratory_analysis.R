library(tidyverse)
library(monocle3)
library(data.table)
library(patchwork)
library(Seurat)
library(SeuratObject)
library(fgsea)
library(loupeR)
library(readxl)

# Directories and Used Functions -----------------------------------------------------------
source('functions/directories_and_theme.R')
source('functions/explore_annotation.R')
source('functions/percent_celltype.R')

#Unique Functions
create_cloupe <- function(object, filename){
  seurat.obj <- SeuratObject::as.Seurat(object, data = NULL)
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj@meta.data$organ_experiment <- do.call(paste, list(seurat.obj$organ,
                                                               seurat.obj$experiment,
                                                               sep = '_'))
  seurat.obj@meta.data$exp_treat <- do.call(paste, list(seurat.obj$experiment,
                                                        seurat.obj$treatment.agg,
                                                        sep = '_'))
  seurat.obj@meta.data$cd45_organ <- do.call(paste, list(seurat.obj$cd45x,
                                                         seurat.obj$organ,
                                                         sep = '_'))
  
  seurat.obj@meta.data <- seurat.obj@meta.data[,c('celltype', 'organ', 'Cluster',
                                                  'ct_cluster', 'treatment.agg',
                                                  'experiment', 'organ_experiment',
                                                  'exp_treat', 'cd45_organ')]
  
  create_loupe_from_seurat(seurat.obj, output_name = filename, force = TRUE)
  invisible(TRUE)
}

extract_coords <- function(cds, cols){
  coords <- as.data.table(reducedDim(cds, 'UMAP'))
  setnames(coords, c('UMAP1', 'UMAP2'))
  meta <- as.data.table(cds@colData)
  meta <- meta[,..cols]
  cbind(coords, meta)
}

plot_coords <- function(dt, col){
  dt %>% 
    ggplot() +
    geom_hex(data = dt[, !..col],
             aes(x = UMAP1, y = UMAP2), fill = 'ivory3', bins = 200) +
    geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) +
    scale_fill_gradient(low = "deepskyblue2", 
                        high = "red", 
                        limits = c(0,40), 
                        oob = scales::squish) +
    facet_wrap(col, nrow = 2) +
    theme_my(axis.ticks = element_blank(),
             panel.grid.major = element_blank(),
             axis.text = element_blank(),
             axis.text.x = element_blank()
    )
}

# Load Data ---------------------------------------------------------------
LUT <- read_rds(file.path(dataDir, '0_doublet_LUT.rds'))
setDT(LUT)

monocle.obj <- read_rds(file.path(dataDir, '4_full_dataset.cds'))
t_subset <- read_rds(file.path(dataDir, '4_t_subset.cds'))
b_subset <- read_rds(file.path(dataDir, '4_b_subset.cds'))
m_subset <- read_rds(file.path(dataDir, '4_m_subset.cds'))

#T Sub Type Markers
new_t_markers <- read_xlsx(file.path(tablesDir, '3d_MarkerGenes_Thelper_ML.xlsx'))
#M Sub Type Markers
new_m_markers <- read_xlsx(file.path(tablesDir, '4_MarkerGenes_Myeloid.xlsx'))

# Plots with new CT assignment -------------------------------------------------------

## T Subset ----------------------------------------------------------------
#first load celltypist annotations
t_celltypist <- load_celltypist(file.path(tablesDir, 
                                          '3c_lymphoid_res_celltypist.csv'))

#then we load SingleR results
t_singler <- load_singler(file.path(tablesDir, 
                                    c('3a_SingleR_res_t_subset.csv',
                                      '3a_SingleR_res_full_dataset.csv')),
                          t_subset)

t_res <- prepare_res(t_singler, t_celltypist, t_subset)

fc_levels <- fct_inorder(t_res[order(celltype), celltype]) %>% levels()

t_res[, celltype := factor(celltype,
                             levels = fc_levels)]

p1 <- plot_automated_anno(t_res, ignore_ref = c('immgen_label.main',
                                                'struc_label.fine',
                                                'th_rest_label.fine'),
                            group = 'celltype')

#CITE-seq distribution
t_abdat <- generate_abdata(t_subset, 
                           fc_levels, 
                           c('cd4cd8', 'cd45x'),
                           group = 'celltype')

p2 <- plot_abs(t_abdat, group = 'celltype')

#Organ
t_org <- prepare_organ_data(t_subset, fc_levels, group = 'celltype')
p3 <- plot_org(t_org, group = 'celltype')

#Top Marker genes
marker_genes <- get_top_markers(t_subset, group_by = 'celltype')

p4 <- plot_genes(t_subset,
                 group_cells_by = 'celltype',
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
                 group_cells_by = 'celltype',
                 cell_factor_levels = fc_levels,
                 markers = new_t_markers) + theme(axis.text.y = element_blank(),
                                                  axis.ticks.y = element_blank(),
                                                  axis.title.y = element_blank()) +
  ggtitle('T subtype marker gene expression') 

#Patchwork plots together
p1+p2+p3+p4+p5+
  plot_layout(widths = c(5.5, 2, 0.5, 4, 14.5))+
  plot_annotation(tag_levels = 'A')

ggsave(file.path(plotsDir, '4_t_annotation_plot.jpg'), 
       height = 6.2, 
       width = 38)

## B Subset ----------------------------------------------------------------
b_celltypist <- load_celltypist(file.path(tablesDir, '3c_b_res_celltypist.csv'))

b_singler <- load_singler(c(file.path(tablesDir, '3a_SingleR_res_b_subset.csv'),
                            file.path(tablesDir, '3a_SingleR_res_full_dataset.csv')),
                          b_subset)

b_res <- prepare_res(b_singler, b_celltypist, b_subset)

fc_levels <- fct_inorder(b_res[order(celltype), celltype]) %>% levels()

p1 <- plot_automated_anno(b_res, use_ref = c('CellTypist_label.fine',
                                             'b_data_label.fine',
                                             'immgen_label.fine',
                                             'old_label.fine'),
                          group = 'celltype')
#Organ
b_org <- prepare_organ_data(b_subset, fc_levels, group = 'celltype')
p2 <- plot_org(b_org, group = 'celltype')

#Top Marker genes
marker_genes <- get_top_markers(b_subset, group_by = 'celltype')

p3 <- plot_genes(b_subset,
                 group_cells_by = 'celltype',
                 cell_factor_levels = fc_levels,
                 markers = setNames(list(marker_genes), 'top_markers'),
                 ordering_type = 'maximal_on_diag'
) + theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none') +
  ggtitle('Top Marker Genes')

p1+p2+p3+
  plot_layout(widths = c(13,1.5,6))+
  plot_annotation(tag_levels = 'A')
ggsave(file.path(plotsDir, '4_b_annotation_plot.jpg'), height = 6, width = 15)

## M Subset ----------------------------------------------------------------
m_celltypist <- load_celltypist(file.path(tablesDir, '3c_myeloid_res_celltypist.csv'))

m_singler <- load_singler(file.path(tablesDir, '3a_SingleR_res_full_dataset.csv'),
                          m_subset)
m_res <- prepare_res(m_singler, m_celltypist, m_subset)

fc_levels <- fct_inorder(m_res[order(celltype), celltype]) %>% levels()

p1 <- plot_automated_anno(m_res,
                          group = 'celltype',
                          use_ref = c('CellTypist_label.fine',
                                      'immgen_label.fine',
                                      'old_label.fine'))
#Organ
m_org <- prepare_organ_data(m_subset, fc_levels, group = 'celltype')
p2 <- plot_org(m_org, group = 'celltype')

#Top Marker genes
marker_genes <- get_top_markers(m_subset, group_by = 'celltype')

p3 <- plot_genes(m_subset,
                 group_cells_by = 'celltype',
                 cell_factor_levels = fc_levels,
                 markers = setNames(list(marker_genes), 'top_markers'),
                 ordering_type = 'maximal_on_diag'
) + theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none') +
  ggtitle('Top Marker Genes')


new_m_markers <- new_m_markers %>% lapply(na.omit)

p4 <- plot_genes(m_subset,
                 group_cells_by = 'celltype',
                 cell_factor_levels = fc_levels,
                 markers = new_m_markers) + theme(axis.text.y = element_blank(),
                                                  axis.ticks.y = element_blank(),
                                                  axis.title.y = element_blank()) +
  ggtitle('M subtype marker gene expression') 

p1+p2+p3+p4+
  plot_layout(widths = c(11,1.5,7,8))+
  plot_annotation(tag_levels = 'A')
ggsave(file.path(plotsDir, '4_m_annotation_plot.jpg'), height = 8, width = 30)

# Cloupe Files ------------------------------------------------------------
create_cloupe(monocle.obj, file.path(gfsDir, 'cloupe_files/full_dataset'))
create_cloupe(t_subset, file.path(gfsDir, 'cloupe_files/t_subset'))
create_cloupe(b_subset, file.path(gfsDir, 'cloupe_files/b_subset'))
create_cloupe(m_subset, file.path(gfsDir, 'cloupe_files/m_subset'))

# Percent Celltype Analysis--------------------------------------------------------

## Prepare & Plot Data ---------------------------------------------------------------
pc_dat <- count_groups(monocle.obj)
pc_dat <- expand_groups(pc_dat)
pc_dat <- calculate_percentages(pc_dat)
stats <- calculate_statistics(pc_dat, group_by = c('raw_celltype',
                                                   'organ',
                                                   'experiment'))
stats_sig <- stats[label != 'NS']

#No Sig Hits
stats_sig

plot_percent(pc_dat, 'raw_celltype')
ggsave(file.path(plotsDir, "4_percent_celltype_broad_celltypes.pdf"), 
       height = 5, width = 9, scale = 3)


# UMAPS - celltype -------------------------------------------------------
plot_cells(monocle.obj, 
           color_cells_by = 'celltype',
           group_label_size = 3,
           label_groups_by_cluster = F)
ggsave(file.path(plotsDir, '4_full_dataset_umap_celltype.jpg'))

#broad celltype assignment
plot_cells(monocle.obj, 
           color_cells_by = 'raw_celltype',
           group_label_size = 3,
           label_groups_by_cluster = F)
ggsave(file.path(plotsDir, '4_full_dataset_umap_broad_celltype.jpg'))

plot_cells(b_subset, 
           color_cells_by = 'celltype',
           group_label_size = 3,
           label_groups_by_cluster = F)
ggsave(file.path(plotsDir, '4_b_umap_celltype.jpg'))

plot_cells(t_subset, 
           color_cells_by = 'celltype',
           group_label_size = 3,
           label_groups_by_cluster = F)
ggsave(file.path(plotsDir, '4_t_umap_celltype.jpg'))

plot_cells(m_subset, 
           color_cells_by = 'celltype',
           group_label_size = 3,
           label_groups_by_cluster = F)
ggsave(file.path(plotsDir, '4_m_umap_celltype.jpg'))

# UMAPS - CITE seq --------------------------------------------------------

## Full Dataset ------------------------------------------------------------
coords <- extract_coords(monocle.obj, c('experiment',
                                        'cd45x',
                                        'cd4cd8',
                                        'hashid',
                                        'gdcd4',
                                        'treatment.agg'))

plot_coords(coords, 'cd45x')
ggsave(file.path(plotsDir, "4_full_dataset_umap_cd45.jpg"))

plot_coords(coords, 'cd4cd8')
ggsave(file.path(plotsDir, "4_full_dataset_umap_cd4cd8.jpg"))

plot_coords(coords, 'treatment.agg')
ggsave(file.path(plotsDir, "4_full_dataset_umap_treatment_agg.jpg"))

plot_coords(coords, 'gdcd4')
ggsave(file.path(plotsDir, "4_full_dataset_umap_gdcd4.jpg"))

## T Subset ----------------------------------------------------------------
coords_t <- extract_coords(t_subset,
                           c('experiment',
                             'cd45x',
                             'cd4cd8',
                             'gdcd4',
                             'treatment.agg',
                             'organ'
                           ))

plot_coords(coords_t, 'cd4cd8')
ggsave(file.path(plotsDir, "4_t_subset_umap_cd4cd8.jpg"))

plot_coords(coords_t, 'cd45x')
ggsave(file.path(plotsDir, "4_t_subset_umap_cd45x.jpg"))

plot_coords(coords_t, 'experiment')
ggsave(file.path(plotsDir, "4_t_subset_umap_exp.jpg"))

plot_coords(coords_t, 'treatment.agg')
ggsave(file.path(plotsDir, "4_t_subset_umap_treatment_agg.jpg"))

plot_coords(coords_t, 'organ')
ggsave(file.path(plotsDir, "4_t_subset_umap_organ.jpg"))


## M Subset ----------------------------------------------------------------
coords_m <- extract_coords(m_subset,
                           c('experiment',
                             'cd45x',
                             'cd4cd8',
                             'gdcd4',
                             'treatment.agg',
                             'organ'
                           ))
plot_coords(coords_m, 'organ')
ggsave(file.path(plotsDir, "4_m_subset_umap_organ.jpg"))

## B Subset ----------------------------------------------------------------
coords_b <- extract_coords(b_subset,
                           c('experiment',
                             'cd45x',
                             'cd4cd8',
                             'gdcd4',
                             'treatment.agg',
                             'organ'
                           ))
plot_coords(coords_b, 'cd4cd8')
ggsave(file.path(plotsDir, "4_b_subset_umap_cd4cd8.jpg"))
