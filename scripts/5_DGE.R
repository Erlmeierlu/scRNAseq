library(tidyverse)
library(monocle3)
library(data.table)
library(limma)
library(edgeR)
library(fst)
library(Matrix.utils)

# Directories and Used Functions -----------------------------------------------------------
source('functions/directories_and_theme.R')

#Unique Functions
make_pseudobulk <- function(cds, 
                            cell_groups, 
                            func = c('mean', 'sum'),
                            return = c('both', 'pseudobulk', 'percentage')) {
  
  func <- match.arg(func)
  return <- match.arg(return)
  if(return == 'percentage') func <- 'sum'
  
  cds <- cds[, cds$treatment.agg != 'Undefined']
  groups <- colData(cds)[, cell_groups] %>% 
    droplevels
  
  dge_groups <- base::do.call(paste, c(as.data.frame(groups), sep = '_'))
  
  c <- counts(cds)
  
  #The next 10 lines are taken from edgeR::filterByExpr and edgeR::cpm
  #Those functions are very slow due to the conversion of sparse count matrices
  #to 'normal' matrices. Since I will not use log transformation here the results
  #are the same. This is way faster now.
  lib.size <- colSums(c)
  dge_groups <- as.factor(dge_groups)
  n <- tabulate(dge_groups)
  MinSampleSize <- min(n[n > 0L])
  MedianLibSize <- median(lib.size)
  CPM.Cutoff <- 10/MedianLibSize * 1e+06
  CPM <- t(t(c)/lib.size*1e6)
  tol <- 1e-14
  keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol)
  keep.TotalCount <- (rowSums(c) >= 15 - tol)
  gfilter <- keep.CPM & keep.TotalCount
  
  c <- c[gfilter, ]
  
  c <- t(c)
  
  pb <- aggregate.Matrix(c,
                         groupings = groups)
  if(return == 'pseudobulk' & func == 'sum') return(t(pb))
  
  pb_c <- aggregate.Matrix(c, 
                           groupings = groups, 
                           fun = 'count') 
  if(func == 'mean') {
    pb@x <- pb@x/pb_c@x
  }
  
  if(return == 'pseudobulk') return(t(pb))
  
  size <- groups %>% 
    as.data.frame() %>% 
    unite(groupings, everything(), sep = '_') %>% 
    group_by(groupings) %>%  
    summarize(n=n())
  
  size <- size[match(rownames(pb_c), size$groupings),]
  stopifnot(all(rownames(pb_c) == size$groupings))
  pc <- pb_c/size$n
  
  pb@x[is.na(pb@x)] <- 0
  if(return == 'percentage') return(t(pc))
  
  list(expression = t(pb), percentage = t(pc))
}

make_metadata <- function(cds, counts, use_ct_cluster = FALSE) {
  meta <- cds@colData %>% as.data.table()
  if(isTRUE(use_ct_cluster)) meta[, celltype := ct_cluster]
  meta <- meta[treatment.agg != 'Undefined'] %>% droplevels()
  meta <-
    meta[,.(celltype, treatment.agg, sample, sex, organ, experiment)
    ][, 
      .N, 
      keyby = .(celltype, treatment.agg, sample, sex, organ, experiment)
    ]
  
  meta[, rn := paste(sample, celltype, treatment.agg, sep = "_")]
  meta <- meta[N > 3]

  meta <- meta[match(colnames(counts), meta$rn),]
  meta <- na.omit(meta)
  setnames(meta, 'treatment.agg', 'treatment')
  meta
}

perform_dge <- function(metadata, 
                        counts, 
                        group_by = NULL, 
                        save_data_for_shiny = FALSE,
                        model_term = 'treatment'){
  
  apply_limma <- function(x, by = NULL){
    if(nrow(x) < 2) return(data.table())
    if(uniqueN(x[,t, env = list(t = model_term)]) < 2) return(data.table())
    
    counts_subset <- counts[, x$rn]
    
    model1 <- formula(paste0('~0+', model_term, '+sex'))
    model2 <- formula(paste0('~0+', model_term))
    
    frame <- model.frame(model1, x, 
                         drop.unused.levels = TRUE)
  
    tryCatch(
      designmat <- model.matrix(model1, frame),
      error = function(e){
        designmat <<- model.matrix(model2, frame)
      })
    
    colnames(designmat) <-
      colnames(designmat) %>%
      str_remove_all("\\(|\\)|celltype|organ|experiment|treatment|\\/|[0-9]\\/[0-9]|sex") %>%
      str_replace_all("\\-", "_")
    
    counts_subset <- DGEList(counts_subset)
    counts_subset <- calcNormFactors(counts_subset)
    
    voomx <- voomWithQualityWeights(counts_subset,
                                    design = designmat,
                                    plot = FALSE,
                                    save.plot = FALSE)
    
    if(save_data_for_shiny == TRUE){
      out <- as.data.table(voomx$E, keep.rownames = 'rn')
      
      filename <- paste(lapply(by, as.character), collapse = '_')
      write_fst(out, file.path(
        shinyDir,
        'data',
        paste0('5_counts_', filename, '.fst')
      ),
      compress = 0)
    }
    
    contrast.matrix <-
      makeContrasts(contrasts = combn(rev(colnames(designmat[,colnames(designmat) != 'M'])),
                                      2,
                                      FUN = paste,
                                      collapse = "-"),
                    levels = designmat)
    
    tryCatch(
      fitx <-
        voomx %>%
        lmFit(design = designmat) %>%
        contrasts.fit(contrast.matrix) %>%
        eBayes(
        ), error = function(e) fitx <<- FALSE
    )
    
    if(isFALSE(fitx)) return(data.table())
    
    map(set_names(colnames(contrast.matrix)), function(x) {
      topx <- topTable(fitx, number = Inf, coef = x)
      as.data.table(topx, keep.rownames = 'rn')
    }) %>%
      rbindlist(idcol = 'coef')
  }
  
  metadata <- metadata[, keyby = group_by, apply_limma(.SD, .BY)]
  metadata[,coef := str_replace_all(coef, '-', '_vs_')]
  
  metadata[, direction := .(fcase(logFC < -1 & adj.P.Val < 0.05, "down",
                                  logFC > 1 & adj.P.Val < 0.05, "up",
                                  default = "NS"))]
  
  metadata[adj.P.Val == 0, adj.P.Val := 5e-324]
  metadata[P.Value == 0, P.Value := 5e-324]
  metadata[]
}
# Load Data ---------------------------------------------------------------

monocle.obj <- read_rds(file.path(dataDir, "4_full_dataset.cds"))
t_subset <- read_rds(file.path(dataDir, '4_t_subset.cds'))
b_subset <- read_rds(file.path(dataDir, '4_b_subset.cds'))
m_subset <- read_rds(file.path(dataDir, '4_m_subset.cds'))

# Full Dataset -------------------------------------------------

## Create Pseudo Bulk ------------------------------------------------------
pb_res <- make_pseudobulk(monocle.obj, 
                          c("sample", "celltype", "treatment.agg"),
                          func = 'sum',
                          return = 'pseudobulk')

if(!is.list(pb_res)){
  countsx <- pb_res
}else{
  countsx <- pb_res$expression
}

metax <- make_metadata(monocle.obj, countsx)

#We probably need those in the shiny app
write_rds(countsx, file = file.path(dataDir, "5_counts_pseudobulk.rds"))
write_rds(metax, file = file.path(shinyDir, "data/5_design_pseudobulk.rds"))

## Res -------------------------------------------------------------------

# metax <- read_rds(file.path(shinyDir, 'data/5_design_pseudobulk.rds'))
# countsx <- read_rds(file.path(dataDir, '5_counts_pseudobulk.rds'))

res <- perform_dge(metax, 
                   countsx, 
                   c('celltype', 'organ', 'experiment'),
                   save_data_for_shiny = TRUE)

res[, treatment := .(factor(
  fcase(
    coef %in% "HDAC_WT_vs_NoT", "WT_vs_ctrl", 
    coef %in% "HDAC_cKO_vs_NoT", "KO_vs_ctrl",
    coef %in% "HDAC_cKO_vs_HDAC_WT", "KO_vs_WT"
  ), 
  levels = c("WT_vs_ctrl", "KO_vs_ctrl", "KO_vs_WT")
))]

write_rds(res, file.path(resDir, "5_dge_full_dataset.rds"))
fwrite(res, file.path(resDir, '5_dge_full_dataset.csv'))

res <- read_rds(file.path(resDir, '5_dge_full_dataset.rds'))
res[, ':='(start = fcase(grepl('ctrl', treatment), 'NoT',
                         default = 'WT'),
           end = fcase(grepl('KO', treatment), 'cKO',
                       default = 'WT'))]

cols <- c('coef', 't', 'B', 'direction')
# filename <- paste(lapply(by, as.character), collapse = '_')
res[, 
    keyby = .(celltype, organ, experiment), 
    write_fst(.SD, 
              path = file.path(shinyDir, 'data', paste0(
                '5_',
                paste(.BY, collapse = '_'), 
                '.fst')
              ),
              compress = 0),
    .SDcols = !cols
]

list <- res[,.(rn, celltype, organ, experiment, 'direction2' = direction,
               treatment, t)]
write.fst(list, file.path(shinyDir, 'data', '5_list_for_shiny.fst'))

design <- res[,
              keyby = .(celltype, organ, experiment),
              .(AveMin = min(round(AveExpr)),
                AveMax = max(round(AveExpr)))
][
  , AveVal := AveMin
  ]

write_rds(design, file.path(shinyDir, 'data/5_design.rds'))

# Organs in HDAC2_WT --------------------------------------------------------------

## Res -------------------------------------------------------------------

# metax <- read_rds(file.path(shinyDir, "data/5_design_pseudobulk.rds"))
# countsx <- read_rds(file.path(dataDir, "5_counts_pseudobulk.rds"))

#wt only skin vs ln in hdac2
organ_meta <- metax[experiment == 'HDAC2' & treatment == 'HDAC_WT']
res_organ <- perform_dge(organ_meta, 
                         countsx,
                         group_by = 'celltype',
                         model_term = 'organ')


write_rds(res_organ, file.path(resDir, "5_dge_organ.rds"))

# res_organ <- read_rds(file.path(resDir, "5_dge_organ.rds"))

# CD45.1 Clusters -----------------------------------------------------
## Create Subset -----------------------------------------------------------
cd45_sub <- t_subset[,grepl('CD45', t_subset$celltype) & t_subset$organ == 'LN']
cd45_sub@colData <- cd45_sub@colData %>% droplevels

## Create Pseudobulk ---------------------------------------------------------
pb_sub <- make_pseudobulk(cd45_sub, 
                          c("sample", "celltype", "treatment.agg"),
                          func = 'sum',
                          return = 'pseudobulk')

if(!is.list(pb_res)){
  countsx_sub <- pb_sub
}else{
  countsx_sub <- pb_sub$expression
}

meta_sub <- make_metadata(cd45_sub, countsx_sub)

## Res -------------------------------------------------------------------
res_sub <- perform_dge(meta_sub, countsx_sub, model_term = 'celltype')
write_rds(res_sub, file.path(resDir, "5_dge_cd45.rds"))

