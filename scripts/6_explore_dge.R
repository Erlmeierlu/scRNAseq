library(tidyverse)
library(data.table)
library(patchwork)
library(fgsea)

# Directories and Used Functions -----------------------------------------------------------
source('functions/directories_and_theme.R')
source('functions/percent_celltype.R')
source('functions/convert_genes.R')

#Unique Functions
ccc <- function(x, y){
  p <- cor(x, y, use = 'pairwise')
  vx <- var(x, na.rm = T)
  vy <- var(y, na.rm = T)
  sx <- sd(x, na.rm = T)
  sy <- sd(y, na.rm = T)
  mx <- mean(x, na.rm = T)
  my <- mean(y, na.rm = T)
  
  pc <- (2*p*sx*sy)/(vx+vy+(mx-my)**2)
  pc
}

my_fgsea <- function(stats, pathways) {
  y <- lapply(stats,
              function(x) {
                lapply(
                  pathways,
                  fgsea,
                  stats = x,
                  minSize = 1,
                  nproc = 1
                )
              })
  rbindlist(lapply(y,
                   rbindlist,
                   idcol = "Database"),
            idcol = "cell")
}

# Load Data ------------------------------------------------------------
monocle.obj <- read_rds(file.path(dataDir, "3_annotated_monocle.cds"))
t_subset <- read_rds(file.path(dataDir, '3_t_subset.cds'))
b_subset <- read_rds(file.path(dataDir, '3_b_subset.cds'))
m_subset <- read_rds(file.path(dataDir, '3_m_subset.cds'))

#DGE results
res <- read_rds(file.path(resDir, "res.rds"))
setDT(res)

#psoriasis data
psoriasis_ref <- fread(file.path(tablesDir, 'new_ludwig_stats_my.csv'))[,-1]
psoriasis_markers <- readxl::read_xlsx(file.path(tablesDir, 
                                                 'MarkerGenes_Psoriasis_ML.xlsx'))

#enrichment psoriasis datasets
files <- list.files(tablesDir, pattern = '^GeneSet')

datasets <- sapply(files, function(x) {
  fread(file.path(tablesDir, x), header = FALSE)
})

# Percent Celltype with DGEs--------------------------------------------------------
## HDAC2 WT vs NoT ------------------------------------------------------------

pc_dat <- count_groups(monocle.obj, group_col = 'celltype')
pc_dat <- expand_groups(pc_dat, group_col = 'celltype')
pc_dat <- calculate_percentages(pc_dat, group_col = 'celltype')

pc_dat <- pc_dat[experiment == 'HDAC2' & treatment.agg != 'HDAC_cKO']

stats <- calculate_statistics(pc_dat, group_by = c('celltype',
                                                   'organ',
                                                   'experiment'))

stats_sig <- stats[label != 'NS']

p1 <- plot_percent(pc_dat, group_by = 'celltype') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

res_percent <- copy(res)

res_percent <- res_percent[,
                           keyby = .(celltype, organ, experiment, coef),
                           direction != 'NS'
][,
  keyby = .(celltype, organ, experiment, coef),
  sum(V1)]

res_percent <-
  res_percent[CJ(celltype = res_percent$celltype, organ, experiment, coef, unique = TRUE),
              on = .(celltype, organ, experiment, coef)] %>%
  droplevels 

setnafill(res_percent, fill = 0, cols = 'V1')

res_percent[,
            celltype := forcats::fct_inorder(celltype)]

res_percent <- res_percent[!grepl('^Un|remo', celltype)
][
  celltype %in% pc_dat$celltype
][
  experiment == 'HDAC2'
][
  coef == 'HDAC_WT_vs_NoT']

p2 <- res_percent %>% 
  ggplot(aes(x = celltype, y = V1), drop = F) +
  geom_line(group = 1, col = 'red') +
  facet_wrap(organ ~ experiment, ncol = 2) +
  theme_my(strip.background = element_blank(),
           strip.text.x = element_blank(),
           axis.text.x = element_text(size = 20),
           axis.text.y = element_text(size = 20),
           legend.title = element_text(size = 20),
           legend.text = element_text(size = 20)
  ) +
  scale_y_log10() +
  xlab("Celltype") +
  ylab("n DEGs")  

(p1 / p2) +
  plot_layout(guides = "collect",
              heights = c(10, 3)) & 
  theme(legend.position = 'right')

ggsave(file.path(plotsDir, 'perc_ct_dge_HDAC2_WT_NoT.pdf'), height = 5, width = 12, scale = 3)

# Cor HDAC1 HDAC2 ---------------------------------------------------------

cors <- res[treatment == "WT_vs_ctrl",
            .(celltype, organ, experiment, treatment, logFC, t, rn)
][,
  dcast(.SD, 
        organ + rn + celltype ~ experiment, 
        value.var = c("logFC", "t"))]

cors[, 
     keyby = .(organ, celltype), .(
       logFC = ccc(logFC_HDAC1, logFC_HDAC2),
       t = ccc(t_HDAC1, t_HDAC2))] %>% 
  na.omit() %>% 
  melt(id = c(1, 2),
       value.name = "cor") %>%
  ggplot() +
  geom_col(aes(celltype, cor), col = "black", fill = "orchid4") +
  facet_grid(variable ~ organ, scales = "free") +
  theme_my() +
  ggtitle("HDAC1/2 WT Concordance Correlation Coefficient")

ggsave(file.path(plotsDir, "cor_hdacs_wt.pdf"))

c <- ccc(cors$logFC_HDAC1, cors$logFC_HDAC2) %>% round(2)

maxi <- max(cors$logFC_HDAC1, cors$logFC_HDAC2, na.rm = T)
mod <- lm(logFC_HDAC2 ~ logFC_HDAC1, data = cors)
intercept <- mod$coefficients[[1]] %>% round(2)
slope <- mod$coefficients[['logFC_HDAC1']] %>% round(2)

ggplot(cors,
       aes(logFC_HDAC1, logFC_HDAC2)) +
  coord_fixed(xlim = c(-maxi, maxi),
              ylim = c(-maxi, maxi)) + 
  geom_abline(aes(linetype = c('dashed', 'solid'),
                  intercept = c(0, intercept),
                  slope = c(1, slope)),
              linewidth = 0.2,
              data = cors[1:2,]) +
  geom_point(col = 'grey75',
             fill = NA,
             size = 0.2,
             alpha = 0.1) +
  theme_my(legend.position = "bottom",
           axis.text.x = element_text(angle = 0,
                                      hjust = 0.5,
                                      vjust = 1)) + 
  annotate('text', 
           label = paste0('CCC=', c),
           x=Inf,
           y=-Inf,
           hjust=1.05,
           vjust=-0.4) +
  scale_linetype_manual(values = c('dashed', 'solid'), 
                        labels = c('y=x', 
                                   paste0('fit: y=', slope,'x-', abs(intercept)))) +
  labs(linetype = element_blank(),
       x = expression(paste(log[2], '(FC) HDAC1-WT')),
       y = expression(paste(log[2], '(FC) HDAC2-WT')))

ggsave(file.path(plotsDir, "hdac_wt_ccc.jpg"))
# Comparison to Psoriasis Data --------------------------------------------

psoriasis_ref <- psoriasis_ref[, melt(.SD, 
                                      id = 'cluster', 
                                      measure(
                                        treatment,
                                        value.name,
                                        sep = '_'))
][,
  .(celltype = cluster,
    condition = treatment,
    rn = str_to_title(n),
    adj_p = p)]

psoriasis_ref <- na.omit(psoriasis_ref)

psoriasis_ref[, ct_rn := do.call(paste, .(celltype, rn, sep = '_'))]

psoriasis_ref <- psoriasis_ref[condition == 'Psoriasis']

psoriasis_kc_genes <- psoriasis_markers$KCs %>% na.omit()
psoriasis_fb_genes <- psoriasis_markers$Fibs %>% na.omit()

plot_genes_by_group(monocle.obj[,monocle.obj$raw_celltype == 'Fibroblasts'],
                    psoriasis_fb_genes,
                    max.size = 3,
                    group_cells_by = 'ct_cluster') + 
  ggtitle('Fibroblasts Psoriasis Markers')

ggsave(file.path(plotsDir, 'FB_Psoriasis_markers.pdf'))

plot_genes_by_group(monocle.obj[,monocle.obj$raw_celltype == 'Keratinocytes'],
                    psoriasis_kc_genes,
                    max.size = 3,
                    group_cells_by = 'ct_cluster') + 
  ggtitle('Keratinocytes Psoriasis Markers')

ggsave(file.path(plotsDir, 'KC_Psoriasis_markers.pdf'))

## Enrichment of models ----------------------------------------------------
#convert symbols to mouse symbols by custom function
datasets <- lapply(datasets, convert_gene_list)
names(datasets) <- str_remove_all(names(datasets), '^GeneSet|.txt.V1$')

res_fgsea <- copy(res)[treatment == 'WT_vs_ctrl' & grepl('Kera|Fibr', celltype)]
res_fgsea <- res_fgsea[, .(data = .(.SD)),
                       by = .(celltype,
                              organ,
                              experiment)]

res_fgsea[, data := lapply(res_fgsea$data, function(x) setNames(x[,t], x[,rn]))]

res_fgsea <- setNames(res_fgsea$data,
                      paste(res_fgsea$celltype,
                            res_fgsea$organ,
                            res_fgsea$experiment,
                            sep = "-"))

res_fgsea <-
  my_fgsea(res_fgsea, list(datasets))

res_fgsea[, c('celltype',
              'organ',
              'experiment') := tstrsplit(cell, "-")]

fgsea_filtered <- res_fgsea[padj < 0.05]

res_fgsea <- res_fgsea[pathway %in% fgsea_filtered$pathway]

res_fgsea[, leadingEdge := vapply(leadingEdge,
                                  paste, 
                                  character(1),
                                  collapse = ",")]

res_fgsea[, let(
  database = fifelse(grepl('KC', pathway), 'KCs', 'FBs'),
  raw_celltype = fifelse(grepl('Kerat', celltype), 'KCs', 'FBs')
)]

p1 <- res_fgsea[raw_celltype == 'KCs' & database == 'KCs' & organ == 'Skin'] %>% 
  ggplot( aes(celltype, pathway)) +
  geom_point(aes(col = NES, size = pmin(5,-log10(padj)))) +
  theme_my() +
  scale_color_gradient2(
    low = "#005AB5",
    mid = "#EEEEE0",
    high = "#DC3220",
    midpoint = 0
  ) + ggtitle('Keratinocytes') + 
  facet_wrap(facets = 'experiment') +
  xlab('') +
  ylab('')

p2 <- res_fgsea[raw_celltype == 'FBs' & database == 'FBs' & organ == 'Skin'] %>% 
  ggplot( aes(celltype, pathway)) +
  geom_point(aes(col = NES, size = pmin(5,-log10(padj))), show.legend = F) +
  theme_my() +
  scale_color_gradient2(
    low = "#005AB5",
    mid = "#EEEEE0",
    high = "#DC3220",
    midpoint = 0
  ) + ggtitle('Fibroblasts') + 
  facet_wrap(facets = 'experiment') +
  ylab('')

p1 / p2 + plot_layout(guides = "collect") & 
  theme(legend.position = 'right')

ggsave(file.path(plotsDir, 'signature_enrichment_hdac1_hdac2.pdf'), width = 7, height = 8)


# DISCONTINUED ------------------------------------------------------------
#Tried to plot a huige heatmap here. Im not sure anymore what I tried to achieve
#Maybe you can look into it and decide whether its important (probably not)

#' files <- list.files(file.path(
#'     shinyDir,
#'     'data'), pattern = '^counts', full.names = T)
#' 
#' cts <- unique(res$celltype)
#' cts <- cts[order(cts, decreasing = TRUE)]
#' cts <- paste(cts, collapse = '|')
#' exps <- c('HDAC1|HDAC2')
#' orgs <- c('LN|Skin')
#' 
#' res_all <- sapply(files, fst::read.fst, as.data.table = TRUE, simplify = FALSE)
#' 
#' res_all <- lapply(seq_along(res_all), function(i, x) {
#'     wt_cols = grep("WT", names(x[[i]]), value = TRUE)
#'     ko_cols = grep("cKO", names(x[[i]]), value = TRUE)
#'     ctrl_cols = grep("NoT", names(x[[i]]), value = TRUE)
#'     
#'     ct <- names(x)[[i]] %>% str_extract(cts)
#'     exp <- names(x)[[i]] %>% str_extract(exps)
#'     org <- names(x)[[i]] %>% str_extract(orgs)
#'     
#'     cn <- do.call(paste,
#'                   list(ct,
#'                        org,
#'                        exp,
#'                        c('HDAC_WT', 'HDAC_cKO', 'NoT'),
#'                        sep = '_'))
#'     
#'     x[[i]][, c(cn) := .(rowMeans(.SD[, ..wt_cols], na.rm = TRUE),
#'                         rowMeans(.SD[, ..ko_cols], na.rm = TRUE),
#'                         rowMeans(.SD[, ..ctrl_cols], na.rm = TRUE))]
#'     
#'     y <- x[[i]][, c('rn', ..cn)]
#'     y
#' }, x = res_all)
#' 
#' res_all <- Reduce(function(x, y) merge(x,y, by = 'rn', all = T), res_all)
#' res_all <- as.data.frame(res_all)
#' rownames(res_all) <- res_all[,'rn']
#' res_all <- res_all[,-1]
#' 
#' # Complex Heatmap ---------------------------------------------------------
#' perc <- make_pseudobulk(monocle.obj,
#'                         cell_groups = c('celltype', 'organ', 'experiment', 'treatment.agg'),
#'                         func = 'mean',
#'                         return = 'percentage')
#' res_all <- res_all[colnames(res_all) %in% colnames(perc)]
#' 
#' 
#' use <- res_all[str_to_upper(rownames(res_all)) %in% str_to_upper(p_genes), grepl('HDAC1', colnames(res_all))]
#' mat <- t(scale(t(use)))
#' 
#' use_pc <- perc[str_to_upper(rownames(perc)) %in% str_to_upper(rownames(mat)),]
#' use_pc <- use_pc[, colnames(use_pc) %in% colnames(mat)]
#' use_pc <- use_pc[,match(colnames(mat), colnames(use_pc))]
#' use_pc <- use_pc[match(rownames(mat), rownames(use_pc)),]
#' 
#' stopifnot(dim(mat) == dim(use_pc))
#' stopifnot(all(colnames(mat) == colnames(use_pc)))
#' stopifnot(all(rownames(mat) == rownames(use_pc)))
#' 
#' dis_row <- dist(mat)
#' dis_col <- dist(t(mat))
#' cl <- tryCatch(hclust(dis_row),
#'                error = function(e) {
#'                    message(e)
#'                    # needs refinement if needed in future.
#'                    b <- as.matrix(dis)
#'                    problems <- which(is.na(b), arr.ind = TRUE)
#'                    tab <- problems[,'row'] %>% unname %>% table
#'                    # some way to remove problematic rows
#'                })
#' cr <- tryCatch(hclust(dis_col),
#'                error = function(e) {
#'                    message(e)
#'                    # needs refinement if needed in future.
#'                    b <- as.matrix(dis_col)
#'                    problems <- which(is.na(b), arr.ind = TRUE)
#'                    tab <- problems[,'row'] %>% unname %>% table
#'                    rem <- tab[tab == max(tab)] %>% names %>% as.integer()
#'                    mat <<- mat[,-rem]
#'                    NULL
#'                    # some way to remove problematic rows
#'                })
#' mat1 <- mat[,grepl('LN', colnames(mat))]
#' ord_fac <- factor(str_extract(colnames(mat1), 'NoT|WT|cKO'), levels = c('NoT', 'WT', 'cKO'))
#' mat1 <- mat1[, order(ord_fac)]
#' n1 <- str_extract(pattern = cts, colnames(mat1))
#' test1 <- cluster_between_groups(mat1, n1)
#' s1 <- str_extract(pattern = cts, colnames(mat1))
#' n_s1 <- s1 %>% unique %>% length
#' pc1 <- use_pc[,colnames(use_pc) %in% colnames(mat1)]
#' pc1 <- pc1[,match(colnames(mat1), colnames(pc1))]
#' stopifnot(all(all(dimnames(pc1)[[1]] == dimnames(mat1)[[1]]),
#'               all(dimnames(pc1)[[2]] == dimnames(mat1)[[2]])))
#' pc1 <- as.matrix(pc1)
#' pc1 <- pc1 + 0.01
#' pc1[is.na(mat1)] <- NA
#' 
#' col_labs1 <- setNames(colnames(mat1) %>% str_extract('[:alpha:]+$'), colnames(mat1))
#' 
#' col_fun <- circlize::colorRamp2(c(min(mat, na.rm = TRUE),
#'                                   0,
#'                                   max(mat, na.rm = TRUE)), 
#'                                 c("#005AB5", "ivory2", "#DC3220"))
#' 
#' 
#' 
#' h1 <- Heatmap(mat1, 
#'               cluster_rows = cl,
#'               column_labels = col_labs1,
#'               col = col_fun,
#'               row_split = 10,
#'               border = FALSE,
#'               row_gap = unit(3, "mm"),
#'               column_gap = unit(3, "mm"),
#'               column_split = n_s1,
#'               cluster_columns = test1,
#'               rect_gp = gpar(type = "none"),
#'               cell_fun = function(j, i, x, y, width, height, fill){
#'                   if(is.na(mat1[i,j])){
#'                       grid.segments(x - 0.125 * width, y, x + 0.125 * width, y,
#'                                     gp = gpar(col = 'grey'))
#'                       grid.text('NA', x, y, gp = gpar(fontsize = 5, col = 'red'))
#'                       # grid.rect(x = x, y = y, width = width, height = height,
#'                       #           gp = gpar(col = "grey", fill = 'grey'))
#'                   }
#'                   # grid.rect(x = x, y = y, width = width, height = height,
#'                   #           gp = gpar(col = "grey70", fill = NA))
#'                   grid.circle(x, y, r = pindex(pc1, i, j) * min(unit.c(width, height)) * 0.5,
#'                               gp = gpar(fill = col_fun(pindex(mat1, i,j)), col = 'grey25'))
#'               })
#' 
#' 
#' mat2 <- mat[,!grepl('LN', colnames(mat))]
#' ord_fac <- factor(str_extract(colnames(mat2), 'NoT|WT|cKO'), levels = c('NoT', 'WT', 'cKO'))
#' mat2 <- mat2[, order(ord_fac)]
#' n2 <- str_extract(pattern = cts, colnames(mat2))
#' test2 <- cluster_between_groups(mat2, n2)
#' s2 <- str_extract(pattern = cts, colnames(mat2))
#' n_s2 <- s2 %>% unique %>% length
#' pc2 <- use_pc[,colnames(use_pc) %in% colnames(mat2)]
#' pc2 <- pc2[,match(colnames(mat2), colnames(pc2))]
#' stopifnot(all(all(dimnames(pc2)[[1]] == dimnames(mat2)[[1]]),
#'               all(dimnames(pc2)[[2]] == dimnames(mat2)[[2]])))
#' pc2 <- as.matrix(pc2)
#' 
#' h2 <- Heatmap(mat2,
#'               cluster_rows = cl,
#'               column_split = n_s2,
#'               row_split = 10,
#'               col = col_fun,
#'               cluster_columns = test2,
#'               rect_gp = gpar(type = "none"),
#'               layer_fun = function(j, i, x, y, width, height, fill){
#'                   if(is.na(pindex(mat2, i, j))){
#'                       grid.rect(x = x, y = y, width = width, height = height,
#'                                 gp = gpar(col = "grey", fill = 'grey'))
#'                   }
#'                   # grid.rect(x = x, y = y, width = width, height = height,
#'                   #           gp = gpar(col = "grey70", fill = NA))
#'                   grid.circle(x, y, r = pindex(pc2, i, j) * min(unit.c(width, height) * 0.5), 
#'                               gp = gpar(fill = col_fun(pindex(mat2, i,j)), col = 'grey25'))
#'               })
#' 
#' h1+h2
#' 
#' #' #' Code snippet from https://github.com/jokergoo/ComplexHeatmap/issues/155#issuecomment-551153988
#' #' #' explained below.
#' #' #' BEWARE: VERY VERY VERY SLOW! LOOPPS THROUGH MATRIX ROWS AND COLS SEVERAL TIMES.
#' #' #' NOT THE BEST SOLUTION FOR LARGE DATASETS. 
#' #' #' USING THIS APPROACH ONCE, BUT MIGHT CHANGE UPSTREAM ANALYSIS, ESPECIALLY
#' #' #' DGE, TO REMOVE THE PROBLEM COMPLETELY. 
#' #' #' 
#' #' #' Hclust cannot handle matrices in which for some pairs of rows and columns,
#' #' #' only 1 or fewer shared values are non-NA. This function recurrently
#' #' #' identifies the most aggravating column/row, excludes that column/row and checks
#' #' #' whether more columns/rows need to be excluded
#' #' #'
#' #' #' @param mat Matrix to investigate
#' #' #' @param min_shared_fields Minimum number of positions that are not NA in both
#' #' #' vectors in order not to flag the vector pair as problematic
#' #' #'
#' #' identify_problematic_combs <- function(mat, min_shared_fields = 1) {
#' #'   exclude_rows <- NULL
#' #'   exclude_cols <- NULL
#' #'   # stopifnot(is.matrix(mat))
#' #'   rownames <- rownames(mat)
#' #'   mat <- setDT(as.data.frame(mat))
#' #'   ## Loop over candidate removals
#' #'   for (k in 1:nrow(mat)) {
#' #'     message('row ', k)
#' #'     candidate_rows <- setdiff(1:nrow(mat), exclude_rows)
#' #'     problem_row_combs <- NULL
#' #'     for (i in candidate_rows) {
#' #'       i_idx <- which(candidate_rows == i)
#' #'       message('candidate ', i)
#' #'       for (j in candidate_rows[i_idx:length(candidate_rows)]) {
#' #'         if (sum(!is.na(mat[i, ]) & !is.na(mat[j, ])) <= min_shared_fields) {
#' #'           message('candidate is ', j)
#' #'           problem_row_combs <- rbind(problem_row_combs, c(i, j))
#' #'         }
#' #'       }
#' #'     }
#' #'     if (is.null(problem_row_combs)) break
#' #'     exclude_rows <- c(exclude_rows,
#' #'                       as.integer(names(which.max(table(problem_row_combs)))))
#' #'   }
#' #'   apply(mat,2, function(x){
#' #'    p <- !is.na(x)
#' #'    sum(p)
#' #'   })
#' #'   apply(mat,1, function(x){
#' #'     p <- !is.na(x)
#' #'     sum(p)
#' #'   })
#' #'   
#' #'   for (k in 1:ncol(mat)) {
#' #'     message('column ', k)
#' #'     candidate_cols <- setdiff(1:ncol(mat), exclude_cols)
#' #'     problem_col_combs <- NULL
#' #'     for (i in candidate_cols) {
#' #'       i_idx <- which(candidate_cols == i)
#' #'       for (j in candidate_cols[i_idx:length(candidate_cols)]) {
#' #'         if (sum(!is.na(mat[, i]) & !is.na(mat[, j])) <= min_shared_fields) {
#' #'           problem_col_combs <- rbind(problem_col_combs, c(i, j))
#' #'         }
#' #'       }
#' #'     }
#' #'     if (is.null(problem_col_combs)) break
#' #'     exclude_cols <- c(exclude_cols,
#' #'                       as.integer(names(which.max(table(problem_col_combs)))))
#' #'   }
#' #'   
#' #'   return(list('row' = exclude_rows, 'column' = exclude_cols))
#' #' }
#' #' 
#' #' 