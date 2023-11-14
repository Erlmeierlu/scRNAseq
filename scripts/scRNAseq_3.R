library(dplyr)
library(stringr)
library(ggplot2)
library(monocle3)
library(data.table)
library(limma)
library(edgeR)

#set seed
set.seed(42)
#personal theme
theme_my <- function() {
  theme(
    panel.grid.major = element_line(colour = "lightgray"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    ),
    panel.spacing.x = unit(10, "mm"),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.key = element_blank(),
    text = element_text(size = 15),
    # strip.text.x = element_text(size = 10, margin = margin(b = 2, t = 2)),
    strip.background = element_rect(
      fill = "lightgray",
      colour = "black",
      size = 1
    ),
    axis.text.x = element_text(
      angle = 90,
      size = 8,
      hjust = 1,
      vjust = 0.5
    )
  )
}

#setting up directories
baseDir <- getwd()
rawDir <- ("/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin/raw_data/scRNA_from_BSF/COUNT")
plotsDir <- file.path(baseDir, "plots/")
tablesDir <- file.path(baseDir, "tables/")
dataDir <- file.path(baseDir, "data/")
bulkDir <- "~/Desktop/MasterThesis/myTables/"

# Load Data ---------------------------------------------------------------

monocle.obj <- readRDS(file.path(dataDir, "old/scRNAseq_3_monocle_more_hash_cutoff.cds"))
# monocle.obj_old <- readRDS(file.path("scRNAseq_3_monocle.cds"))

# monocle.obj@colData$celltype <- monocle.obj@colData$celltype %>% str_replace_all(" |\\-", "_")

# DGE Analysis -------------------------------------------------

## Pseudo Bulk -------------------------------------------------------------

### Create Pseudo Bulk ------------------------------------------------------

monocle.obj@colData$tmp <- monocle.obj@colData$Size_Factor
monocle.obj@colData$Size_Factor <- 1
group_df <-
  data.frame(monocle.obj@colData[, c("sample", "treatment.agg", "celltype")] %>% subset(treatment.agg != "Undefined")) %>% droplevels()
setDT(group_df, keep.rownames = "rn")
group_df <- group_df[, .(rn, groups = do.call(paste, c(.SD, sep = "_"))), .SDcols = -"rn"]

countsx <-
        aggregate_gene_expression(
          monocle.obj,
                cell_group_df = group_df,
                norm_method = "size_only",
                scale_agg_values = F,
                cell_agg_fun = "sum"
        )

monocle.obj@colData$Size_Factor <- monocle.obj@colData$tmp
monocle.obj@colData$tmp <- NULL


designx <- monocle.obj@colData %>% as.data.table()
designx <- designx[treatment.agg != "Undefined"] %>% droplevels()
designx <- designx[,.(celltype, treatment.agg, sample, sex, organ, experiment)]
designx <- designx[, .N, keyby = .(celltype, treatment.agg, sample, organ, experiment)]

designx[, rn := paste(sample, treatment.agg, celltype, sep = "_")]

designx <- as.data.frame(designx)
rownames(designx) <- designx[,"rn"]
designx$rn <- NULL

designx <- designx[match(colnames(countsx), rownames(designx)),]
colnames(designx)[2] <- "treatment"

saveRDS(countsx, file = file.path(dataDir, "counts_pseudobulk.rds"))
saveRDS(designx, file = file.path(dataDir, "design_pseudobulk.rds"))

### Limma -------------------------------------------------------------------

designx <- readRDS(file.path(dataDir, "design_pseudobulk.rds"))
countsx <- readRDS(file.path(dataDir, "counts_pseudobulk.rds"))

results.DGE.pseudo <- list()

for (ct in unique(designx$celltype)) {
  message(ct, "-------------------------")
  for (organx in unique(designx$organ)) {
    for (expx in unique(designx$experiment)) {
      subsetx <- designx %>%
        subset(celltype == ct &
                 organ == organx &
                 experiment == expx &
                 N > 3)
      
      if (nrow(subsetx) < 2)
        next
      
      countsx.pseudo <- countsx[, rownames(subsetx)]
      
      if (length(unique(subsetx$treatment)) < 2)
        next
      
      if (length(unique(subsetx$sex)) < 2) {
        frame <-
          model.frame( ~ 0 + treatment,
                       subsetx, drop.unused.levels = T)
        designmat <- model.matrix( ~ 0 + treatment, frame)
        
      } else {
        frame <-
          model.frame( ~ 0 + treatment + sex,
                       subsetx, drop.unused.levels = T)
        designmat <-
          model.matrix( ~ 0 + treatment + sex, frame)
        
      }
      
      if (ncol(designmat) >= nrow(designmat)) {
        # message(pwd, "_", ct, "_", organx, "_", expx,  "---ncol(designmat) >= nrow(designmat)---")
        next
      }
      
      colnames(designmat) <-
        colnames(designmat) %>%
        str_remove_all("\\(|\\)|sample|treatment|\\/|[0-9]\\/[0-9]|sex") %>%
        str_replace_all("\\-", "_")
      
      dge <- DGEList(countsx.pseudo)
      
      # gfilter <-
      #         filterByExpr(
      #                 dge,
      #                 design = designmat,
      #                 large.n = 10,
      #                 min.count = 7,
      #                 min.prop = 0.7,
      #                 min.total.count = 8
      #         )
      
      gfilter <- rowSums(dge$counts) > 0
      
      
      dge <- dge[gfilter, , keep.lib.sizes = FALSE]
      
      dge <- calcNormFactors(dge)
      
      voomx <- voomWithQualityWeights(dge,
                                      design = designmat,
                                      plot = F,
                                      save.plot = F)
      
      contrast.matrix <-
        makeContrasts(contrasts = combn(rev(colnames(designmat)[colnames(designmat) != "M"]),
                                        2,
                                        FUN = paste,
                                        collapse = "-"),
                      levels = designmat)
      
      fitx <- lmFit(voomx, design = designmat)
      fitx <- contrasts.fit(fitx, contrast.matrix)
      fitx <- eBayes(fitx)
      
      for (coefx in colnames(contrast.matrix)) {
        topx <- topTable(fitx, number = Inf, coef = coefx)
        topx$rn <- rownames(topx)
        results.DGE.pseudo[[ct]][[organx]][[expx]][[coefx]] <-
          data.table(topx)
      }
    }
  }
}

### Res ---------------------------------------------------------------------

res <-
  rbindlist(lapply(results.DGE.pseudo, function(l1) {
    rbindlist(lapply(l1, function(l2) {
        rbindlist(lapply(l2, function(l3) {
          rbindlist(l3, idcol = "coef")}), 
          idcol = "experiment")
    }), idcol = "organ")
  }), idcol = "celltype")

res[, treatment := .(factor(
  fcase(
    coef %in% "HDAC_WT-NoT", "WT_vs_ctrl", 
    coef %in% "HDAC_cKO-NoT", "KO_vs_ctrl",
    coef %in% "HDAC_cKO-HDAC_WT", "KO_vs_WT"
  ), 
  levels = c("WT_vs_ctrl", "KO_vs_ctrl", "KO_vs_WT")
))]

res[, direction := .(fcase(logFC < -1 & adj.P.Val < 0.05, "down",
                              logFC > 1 & adj.P.Val < 0.05, "up",
                              default = "NS"))]


saveRDS(res, file.path(dataDir, "res.rds"))

# P_val Histos ------------------------------------------------------------

res <- readRDS(file.path(dataDir, "res.rds"))

res[rn == "Foxp3", .(expr = round(AveExpr), celltype, experiment, organ)] %>% unique %>% 
  ggplot() +
  geom_col(aes(celltype, expr),col = "black", fill = "ivory2") +
  facet_wrap(organ ~ experiment, scales = "free") +
  theme_my() +
  ggtitle("Foxp3 round(AveExpr)")

ggsave(file.path(plotsDir, "Foxp3_expression.pdf"))

res_old <- readRDS(file.path(dataDir, "old/list_for_shiny.rds"))
res_old$celltype <- res_old$celltype %>% str_replace_all("\\-", "_")

res_fox <- res[round(AveExpr) >= 2]

# for(ct in unique(res$celltype)) {
  for (organx in unique(res_fox$organ)) {
    for (expx in unique(res_fox$experiment)) {
      subsetx <- res_fox %>% filter(
        # celltype == ct,
                                organ == organx,
                                experiment == expx,
                                !grepl("Intercept|M", coef))
      
      if(nrow(subsetx) == 0) next
      
      subsetx %>%
        ggplot() +
        geom_histogram(aes(P.Value, fill = factor(round(AveExpr))),
                       bins = 30) +
        facet_wrap(~ treatment) +
        theme_my() +
        theme(panel.grid.major = element_blank()) +
        ggtitle(paste("ave_expr_filter_pseudo_histo",
                      # ct,
                      organx, expx, sep = "_"))
      
      ggsave(file.path(plotsDir,
                       paste0(
                         paste("ave_expr_filter_pseudo_histo",
                               # ct,
                               organx, expx, sep = "_"),
                         ".pdf"
                       )))
    }
  }
# }

# for(ct in unique(res_old$celltype)) {
#   for (organx in unique(res_old$organ)) {
#     for (expx in unique(res_old$experiment)) {
#       subsetx <- res_old %>% filter(celltype == ct,
#                                 organ == organx,
#                                 experiment == expx,
#                                 !grepl("Intercept|M", coef))
#       
#       if(nrow(subsetx) == 0) next
#       
#       subsetx %>%
#         ggplot() +
#         geom_histogram(aes(P.Value, fill = factor(round(AveExpr))),
#                        bins = 30) +
#         facet_wrap(~ treatment) +
#         theme_my() +
#         theme(panel.grid.major = element_blank()) +
#         ggtitle(paste("pseudo_histo",
#                       ct,
#                       organx, expx, sep = "_"))
#       
#       ggsave(file.path(plotsDir,
#                        paste0(
#                          paste("pval_histos/old_pseudo_histo",
#                                ct,
#                                organx, expx, sep = "_"),
#                          ".pdf"
#                        )))
#     }
#   }
# }


# Cor HDAC1 HDAC2 ---------------------------------------------------------

res[treatment == "WT_vs_ctrl",.(celltype, organ, experiment, treatment, logFC, t, rn)] %>% 
  dcast(t + logFC ~ experiment)

cors <- res[treatment == "WT_vs_ctrl",.(celltype, organ, experiment, treatment, logFC, t, rn)] %>% 
  dcast(organ + rn + celltype ~experiment, value.var = c("logFC", "t"))

cors[, keyby = .(organ, celltype), .(
  logFC = cor(logFC_HDAC1, logFC_HDAC2, use = "p", method = "s"),
  t = cor(t_HDAC1, t_HDAC2, use = "p", method = "s")
)] %>%
  na.omit() %>%
  melt(id = c(1, 2),
       value.name = "cor") %>%
    ggplot() +
      geom_col(aes(celltype, cor), col = "black", fill = "orchid4") +
      facet_grid(variable ~ organ, scales = "free") +
      theme_my() +
      ggtitle("HDAC1/2 WT cors old data")

ggsave(file.path(plotsDir, "cor_old_hdacs_wt.pdf"))

# Percent Celltype Analysis--------------------------------------------------------

  percent_ct_df <- monocle.obj@colData %>%
    as_tibble %>%
    filter(treatment.agg != "Undefined") %>% 
    group_by(experiment, organ, treatment.agg, sample,  celltype) %>% 
    summarise(n = n()) %>% 
    mutate(cells = sum(n), pc.cells = n/cells*100, sex = str_extract(sample, "^[A-z]"), batch = str_extract(sample, "\\d$")) 
  
  comp <- combn(unique(levels(percent_ct_df$treatment.agg %>%
      droplevels
  ) %>% 
    rev),
  2,
  FUN = paste,
  collapse = ";")
  
  stats_list <- list()
  for(ct in unique(percent_ct_df$celltype)){
    for(organx in unique(percent_ct_df$organ)){
      for(expx in unique(percent_ct_df$experiment)){
        for(compx in comp){
          subsetx <- percent_ct_df %>% filter(treatment.agg %in% unlist(str_split(compx, ";")),
                                   celltype == ct,
                                   organ == organx,
                                   experiment == expx)
          
            if(length(unique(subsetx$treatment.agg)) < 2) next
            else if(length(unique(subsetx$sex)) < 2){
              mod <- lm(pc.cells ~ treatment.agg, subsetx)
            }
            else mod <- lm(pc.cells ~ treatment.agg + sex, subsetx)
          
          stats <- summary(mod)$coefficient %>% as_tibble(rownames = "end")
          
          stats_fin <- stats %>% filter(grepl("WT|cKO|NoT", end))
          stats_fin$end <- stats_fin$end %>% str_remove("treatment.agg")
          colnames(stats_fin) <- colnames(stats_fin) %>% 
            str_replace_all(" ", "_") %>% 
            str_remove("\\.")
          colnames(stats_fin)[5] <- "p_value"
          
          stats_list[[ct]][[organx]][[expx]][[compx]] <- stats_fin
        }
      }
    }
  }
  
  res.stats <-
    rbindlist(lapply(stats_list, function(l1) {
      rbindlist(lapply(l1, function(l2) {
        rbindlist(lapply(l2, function(l3) {
            rbindlist(l3, idcol = "comparison")}), idcol = "experiment")
      }), idcol = "organ")
    }), idcol = "celltype")
  res.stats$start <- res.stats$comparison %>% str_extract("\\w+$")
  res.stats$comparison <- res.stats$comparison %>% str_replace(";", "_vs_")
  res.stats$adj_p <- res.stats$p_value %>% p.adjust(method = "BH")
  res.stats$label <- case_when(
    res.stats$adj_p >= 0.1 ~ "NS",
    between(res.stats$adj_p, 0.05, 0.1) ~ "*",
    between(res.stats$adj_p, 0.01, 0.049) ~ "**",
    res.stats$adj_p < 0.01 ~ "***"
  )

  
  y_values <-
    percent_ct_df %>% group_by(experiment, organ, celltype, treatment.agg) %>% slice_max(n = 1, pc.cells) %>% pivot_wider(
      names_from = "treatment.agg",
      values_from = "pc.cells",
      id_cols = c("experiment", "organ", "celltype")
    )

  res.stats <-
    left_join(res.stats, y_values, by = c("experiment", "organ", "celltype"))
                                    
  res.stats$position <- case_when(res.stats$start == "NoT" & res.stats$end == "HDAC_cKO" ~
                                    res.stats[,13:15] %>% as.matrix %>%  rowMaxs * 1.3,
                                  res.stats$start == "NoT" & res.stats$end == "HDAC_WT" ~
                                    res.stats[,13:14] %>% as.matrix %>%  rowMaxs * 1.15,
                                  TRUE ~
                                    res.stats[,14:15] %>% as.matrix %>%  rowMaxs * 1.085)
  res.stats$position <- case_when(res.stats$position < 3 ~ res.stats$position + 3,
                                  TRUE ~ res.stats$position)
  
  # stats <- res.stats %>% filter(label != "NS")
  # stats$start <- c(3-0.4/3, 3-0.4/3, 5-0.4/3, 12-0.4/3, 16-0.4/3, 16-0.4/3, 16-0.4/3)
  # stats$end <- c(3+0.4/3, 3, 5+0.4/3, 12, 16+0.4/3, 16, 16)
  # stats$position <- stats$position/10
  
  percent_ct_df %>% 
          ggplot(aes(celltype, pc.cells)) +
          geom_boxplot(aes(fill = treatment.agg), outlier.shape = NA, width = 0.4) +
          geom_point(aes(group = treatment.agg, col = batch, shape = sex), size = 2, position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.4)) + 
          geom_signif(data = stats,
                      aes(xmin = start,
                          xmax = end,
                          annotations = label,
                          y_position = position,
                          ),
                      show.legend = TRUE,
                      manual = T) +
          scale_fill_grey(start = 0.4, name = "Treatment Group") +
          scale_shape_discrete(name = "Sex") +
          scale_color_discrete(name = "Batch") +
          scale_y_log10() +
          facet_wrap(~ experiment + organ, ncol = 1) +
          theme_my() + 
          theme(strip.background = element_blank(),
                strip.text = element_text(hjust = 0,
                                          size = rel(1.5)),
                axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 20)) +
          xlab("Celltype") +
          ylab("% of Cells")
  
  ggsave("pc.pdf", height = 9, width = 7.5, scale = 3)
  
  # ggsave("laia_perc_ct.pdf", scale = 3.5)
 
  
    ggplot(percent_ct_df) +
      geom_jitter(aes(treatment.agg, pc.cells, col = batch, shape = sex), height = 0, width = 0.2) +
      facet_grid(celltype~experiment+organ, switch = "x", scales = "free_y") + 
      geom_signif(data = res.stats %>% filter(label != "NS"),
                  aes(xmin = start, 
                      xmax = end, 
                      annotations = label,
                      y_position = position),
                  show.legend = TRUE,
                  step_increase = 2,
                  vjust = 0.5,
                  manual = T) +
    theme_my() +
    theme(strip.placement = "outside")

  
ggsave(file.path(plotsDir, "perc_ct_sig.pdf"), height = 45)
  
write.csv(percent_ct_df , file.path(tablesDir,"stats_list.csv"), row.names = F)
