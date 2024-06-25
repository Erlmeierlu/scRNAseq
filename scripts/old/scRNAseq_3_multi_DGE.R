library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(monocle3)
library(ggrepel)
library(data.table)
library(limma)
# library(Matrix.utils)
library(edgeR)
library(patchwork)
library(ggsignif)
library(purrr)


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

# Analysis ----------------------------------------------------------------

marker_test_res <- top_markers(monocle.obj, group_cells_by="Cluster", 
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(4, pseudo_R2)
# 
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))


plot_genes_by_group(monocle.obj,
                    top_specific_marker_ids,
                    group_cells_by="ct_cluster",
                    ordering_type="none",
                    max.size=3) +
    theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
    )
    )

ggsave(file.path(plotsDir, "Markers.pdf"), height = 15, width = 10)

coords <- reducedDims(monocle.obj)
coords <- data.frame(coords$UMAP)
colnames(coords) <- c("UMAP1", "UMAP2")
coords$experiment <- colData(monocle.obj)$experiment
coords$cd45x <- colData(monocle.obj)$cd45x
coords$cd4cd8 <- colData(monocle.obj)$cd4cd8
coords$hashid <- colData(monocle.obj)$hashid
coords$treatment <- colData(monocle.obj)$treatment
coords$hashid.agg <- colData(monocle.obj)$hashid.agg
coords$treatment.agg <- colData(monocle.obj)$treatment.agg
coords$gdcd4 <- colData(monocle.obj)$gdcd4

UMAP_projection <- coords %>% dplyr::select(UMAP1, UMAP2)
UMAP_projection <- as_tibble(UMAP_projection, rownames = "Barcode")
colnames(UMAP_projection) <- c("Barcode", "X Coordinate", "Y Coordinate")

rename_barcodes <- function(object){
  
  object$Barcode <- case_when(grepl("fLN_40B2",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-1"),
            grepl("fLN_40B3",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-2"),
            grepl("fLN_41B1",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-3"),
            grepl("fLN_41B2",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-4"),
            grepl("fLN_41B3",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-5"),
            grepl("fLN_B1",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-6"),
            grepl("fSkin_40B2",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-7"),
            grepl("fSkin_40B3",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-8"),
            grepl("fSkin_41B1",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-9"),
            grepl("fSkin_41B2",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-10"),
            grepl("fSkin_41B3",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-11"),
            grepl("fSkin_B1",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-12"),
            grepl("mLN_40B2",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-13"),
            grepl("mLN_40B3",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-14"),
            grepl("mLN_41B1",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-15"),
            grepl("mLN_41B3",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-16"),
            grepl("mLN_41B4",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-17"),
            grepl("mLN_B1",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-18"),
            grepl("mSkin_40B2",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-19"),
            grepl("mSkin_40B3",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-20"),
            grepl("mSkin_41B1",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-21"),
            grepl("mSkin_41B3",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-22"),
            grepl("mSkin_B1",object$Barcode) ~ str_replace(object$Barcode, "\\-1.*","-23"))
  
  return(object)
}

UMAP_projection <- rename_barcodes(UMAP_projection)

write.csv(UMAP_projection, file = file.path(tablesDir, "UMAP-Projection-Ludwig4.csv"), quote = F, row.names = F)

graph_based <- monocle.obj@colData %>% as_tibble(rownames = "Barcode") %>% dplyr::select(Barcode, ct_cluster)
graph_based$Ludwig4 <- paste("Cluster", graph_based$ct_cluster)
graph_based <- graph_based %>% dplyr::select(-ct_cluster)

graph_based <- rename_barcodes(graph_based)

write.csv(graph_based, file = file.path(tablesDir, "Graph-Based-Ludwig4.csv"), quote = F, row.names = F)

treatments <- monocle.obj@colData %>% as_tibble(rownames = "Barcode") %>% dplyr::select(Barcode, treatment.agg)

treatments <- rename_barcodes(treatments)

write.csv(treatments, file = file.path(tablesDir, "Treatments-Ludwig4.csv"), row.names = F)

experiments <- monocle.obj@colData %>% as_tibble(rownames = "Barcode") %>% dplyr::select(Barcode, experiment)

experiments <- rename_barcodes(experiments)

write.csv(experiments, file = file.path(tablesDir, "experiments-Ludwig4.csv"), row.names = F)

exp_treat <- inner_join(experiments,treatments) %>% 
  mutate(exp_treat = paste(experiment, 
                           str_extract(treatment.agg, "[:alpha:]+$"), 
                           sep = "_")) %>% 
  select(Barcode, exp_treat)

write.csv(exp_treat, file = file.path(tablesDir, "exp-treat-Ludwig4.csv"), row.names = F)


ggplot(coords) + 
    geom_hex(data = coords %>% dplyr::select(-cd45x), 
               aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
    geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
    scale_fill_gradient(low = "deepskyblue2", high = "red", limits = c(0,40), oob = scales::squish)+
    facet_wrap( ~cd45x, nrow = 2) +
    theme_my() #+
    #theme(panel.background = element_rect(fill = "black"),
          #panel.grid.major = element_blank())

ggsave(file.path(plotsDir, "UMAP_cd45.jpg"))

ggplot(coords) + 
    geom_hex(data = coords %>% dplyr::select(-cd4cd8), 
             aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
    geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
    scale_fill_gradient(low = "deepskyblue2", high = "red", limits = c(0,40), oob = scales::squish)+
    facet_wrap(~cd4cd8, nrow = 2) +
    theme_my()


ggsave(file.path(plotsDir, "UMAP_cd4cd8.jpg"))

ggplot(coords) + 
    geom_hex(data = coords %>% dplyr::select(-hashid), 
             aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
    geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
    scale_fill_gradient(low = "deepskyblue2", high = "red", limits = c(0,40), oob = scales::squish) +
    facet_wrap(~hashid, ncol = 3) +
    theme_my()

ggsave(file.path(plotsDir, "UMAP_hash.jpg"))

ggplot(coords) + 
  geom_hex(data = coords %>% dplyr::select(-hashid.agg), 
           aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
  geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
  scale_fill_gradient(low = "deepskyblue2", high = "red", limits = c(0,40), oob = scales::squish) +
  facet_wrap(~hashid.agg, ncol = 2) +
  theme_my()

ggsave(file.path(plotsDir, "UMAP_hash_agg.jpg"))

ggplot(coords) + 
    geom_hex(data = coords %>% dplyr::select(-treatment), 
             aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
    geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
    scale_fill_gradient(low = "deepskyblue2", high = "red", limits = c(0,40), oob = scales::squish) +
    facet_wrap( ~treatment, nrow = 2) +
    theme_my()

ggsave(file.path(plotsDir, "UMAP_treatment.jpg"))

ggplot(coords) + 
  geom_hex(data = coords %>% dplyr::select(-treatment.agg), 
           aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
  geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
  scale_fill_gradient(low = "deepskyblue2", high = "red", limits = c(0,40), oob = scales::squish) +
  facet_wrap( ~treatment.agg, nrow = 2) +
  theme_my()


ggsave(file.path(plotsDir, "UMAP_treatment_agg.jpg"))

ggplot(coords) + 
  geom_hex(data = coords %>% dplyr::select(-gdcd4), 
           aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
  geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
  scale_fill_gradient(low = "deepskyblue2", high = "red", limits = c(0,40), oob = scales::squish) +
  facet_wrap( ~gdcd4, nrow = 2) +
  theme_my()

ggsave(file.path(plotsDir, "UMAP_gdcd4.jpg"))

plot_cells(monocle.obj, color_cells_by = "Cluster", group_label_size = 3.5)
ggsave(file.path(plotsDir, "UMAP_cluster.jpg"))

# plot with colors according to new anno

col <- read.csv("hexcodes.csv", header = T)

plot_cells(monocle.obj, group_cells_by = "cluster", 
           color_cells_by = "celltype",
           label_groups_by_cluster = F,
           group_label_size = 4) +
           scale_color_manual(values = test)
           # scale_color_manual(values = c("#D0021B", "#BA98FF", "#F5A623"), 
           #                    labels = c("GC B cells", "Plasma cells", "Pre-plasmablasts"))

        # scale_color_manual(breaks = c("GC B cells", "Plasma cells", "Pre-plasmablasts"), values = c("#D0021B", "#BA98FF", "#F5A623"), labels = c("GC B cells", "Plasma cells", "Pre-plasmablasts"))

ggsave(file.path("../","UMAP_annotated.jpg"), scale = 1.5)

monocle.obj@colData %>% as_tibble() %>%
  ggplot() +
  geom_hline(yintercept = 0.6, col = "deepskyblue2") +
  geom_violin(aes(sample, gd.ratio, fill = experiment)) +
  theme_my() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  ylab("Ratio of GDT / GDT + CD4 reads") +
  scale_fill_viridis_d(option = "G")

ggsave(file.path(plotsDir, "gd_ratio.pdf"))


q# 
# plot_cells(monocle.obj, group_cells_by="cluster", 
#            color_cells_by="cell_types_immgen_label.fine",
#            label_groups_by_cluster=F,
#            group_label_size = 4)
# 
# ggsave(file.path(plotsDir, "ref_immgen_anno.pdf"))


# Celltype Hashtag Reads --------------------------------------------------
#with hash ratio
for(experimentx in unique(monocle.obj@colData$experiment)) {
  for (ratiox in c("hash.ratio", "hash.agg.ratio")) {
    colData(monocle.obj) %>% as_tibble() %>% filter(experiment == experimentx, get(ratiox) > 0.6) %>%
      ggplot() +
      geom_violin(aes(celltype, hash.sum, fill = celltype)) +
      facet_wrap(~ sample, scales = "free_y", ncol = 3) +
      scale_y_log10() +
      theme_my() +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 6
      )) +
      labs(title = paste("Hashtag-Reads per cell type with Hash-Ratio > 60%", experimentx, ratiox))
    
    ggsave(file.path(plotsDir, paste0("Celltype_hashtag_reads_ratio_",experimentx,"_",ratiox,".pdf")))
  }
}


#no hash ratio
for(experimentx in unique(monocle.obj@colData$experiment)) {
colData(monocle.obj) %>% as_tibble() %>% filter(experiment == experimentx) %>% 
    ggplot() +
    geom_violin(aes(celltype, hash.sum, fill = celltype)) +
    facet_wrap(~sample, scales = "free_y", ncol = 2)+
    scale_y_log10() +
    theme_my() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))+
    labs(title = "Hashtag-Reads per cell type")

ggsave(file.path(plotsDir, paste0("Celltype_hashtag_reads_",experimentx,".pdf")))
}
# Reads post-cleanup ------------------------------------------------------
for (experimentx in unique(colData(monocle.obj)$experiment)) {
  for (ratiox in c("hash.ratio", "hash.agg.ratio")) {
    monocle.obj@colData %>% as_tibble() %>% filter(experiment == experimentx) %>%
      ggplot() +
      stat_bin_hex(aes(hash.sum, get(ratiox), fill = log(..count..)),
                   bins = 30,
                   col = "black") +
      facet_wrap(~ sample, ncol = 3) +
      scale_fill_gradient(low = "ivory2", high = "red") +
      scale_x_log10() +
      geom_hline(
        yintercept = 0.6,
        alpha = 0.7,
        col  = "blue",
        linetype = "dashed"
      ) +
      labs(
        x = "Hashtag-Reads",
        y = "Hashtag-Ratio",
        title = paste("Hashtag-Ratio post cleanup", experimentx, ratiox)
      ) +
      theme_my()
    
    ggsave(file.path(
      plotsDir,
      paste0(
        "Hashtag-Ratio_post_cleanup",
        experimentx,
        "_",
        ratiox,
        ".pdf"
      )
    ))
  }
}

    

# Specific celltype reads -------------------------------------------------
ct_plotlist <- list()
for (ct in unique(colData(monocle.obj)$celltype)) {
    ct_plotlist[[ct]] <-
        monocle.obj@colData %>% as_tibble() %>% filter(celltype == ct &
                                                           batch == "B2") %>%
        ggplot() +
        stat_bin_hex(aes(hash.sum, hash.ratio, fill = log(..count..)),
                     bins = 30,
                     col = "black") +
        facet_wrap(~ sample, ncol = 2) +
        scale_fill_gradient(low = "ivory2", high = "red") +
        scale_x_log10() +
        geom_hline(yintercept = 0.6, alpha = 0.7, col  = "blue", linetype = "dashed") +
        labs(x = "Hashtag-Reads",
             y = "Hashtag-Ratio",
             title = paste(ct)) +
        theme_my()
}
names(ct_plotlist) <- str_replace_all(names(ct_plotlist), " |/", "_")

for (plotx in names(ct_plotlist)){
    ggsave(file.path(plotsDir, paste0("Hashtag_Ratio_vs_reads_", plotx,".pdf")), plot = ct_plotlist[[plotx]])
}



# DGE Analysis -------------------------------------------------

## Single Cell ------------------------------------------------------------

### Limma --------------------------------------------------------------------

comp <- combn(unique(levels(
    monocle.obj[, colData(monocle.obj) %>%
                    subset(treatment != "Undefined") %>% row.names]@colData$treatment %>%
        droplevels
)),
2,
FUN = paste,
collapse = ";")

results.DGE <- list()
counts_HM <- list()
for (ct in unique(monocle.obj@colData$celltype)) {
  for (organx in unique(monocle.obj@colData$organ)) {
    for (compx in comp) {
      for (expx in unique(monocle.obj@colData$experiment)) {
        subsetx <- monocle.obj[,
                               colData(monocle.obj) %>%
                                 subset(
                                   celltype == ct &
                                     organ == organx &
                                     treatment.agg %in% unlist(str_split(compx, ";")) &
                                     experiment == expx
                                 ) %>%
                                 row.names]
        
        if (length(unique(subsetx@colData$treatment.agg)) >= 2 &
            length(unique(subsetx@colData$sex)) >= 2) {
          frame <-
            model.frame(~ treatment.agg + sex,
                        subsetx@colData %>% droplevels())
          designmat <- model.matrix(~ treatment.agg + sex, frame)
          
          if (ncol(designmat) >= nrow(designmat)){
            
            rm(list = c("subsetx",
                        "frame",
                        "designmat"))
            gc(reset = T)
            
            next
          }
          
        } else if (length(unique(subsetx@colData$treatment.agg)) >= 2 &
                   length(unique(subsetx@colData$sex)) < 2) {
          frame <-
            model.frame(~ treatment.agg, subsetx@colData %>% droplevels())
          designmat <- model.matrix(~ treatment.agg, frame)
          
          if (ncol(designmat) >= nrow(designmat)){
            
            rm(list = c("subsetx",
                        "frame",
                        "designmat"))
            gc(reset = T)
            
            next
          }
          
        } else {
          rm(list = c("subsetx"))
          gc(reset = T)
          
          next
        }
        colnames(designmat) <-
          colnames(designmat) %>%
          str_remove_all("\\(|\\)|sample|treatment.agg|\\/|[0-9]\\/[0-9]|sex") %>%
          str_replace_all("\\-", "_")
        
        countsx <- counts(subsetx)
        
        # gfilter <-
        #   filterByExpr(
        #     countsx,
        #     design = designmat,
        #     min.prop = 0.1,
        #     min.count = 1
        #   )
        
        gfilter <- rowSums(countsx) > 0
        
        voomx <-
          voom(
            countsx[gfilter,],
            design = designmat,
            plot = F,
            save.plot = T
          )
        
        #counts_HM[[ct]][[organx]][[compx]][[expx]] <- voomx$E
        
        ggplot() +
          geom_point(aes(voomx$voom.xy$x, voomx$voom.xy$y), size = 0.1) +
          geom_line(aes(voomx$voom.line$x, voomx$voom.line$y), col = "red") +
          xlab(voomx$voom.xy$xlab) +
          ylab(voomx$voom.xy$ylab) +
          ggtitle("voom: Mean-variance trend") +
          theme_my()
        
        
        ggsave(file.path(
          plotsDir,
          paste0(
            "voom_",
            str_replace_all(paste0(unlist(
              str_split(compx, ";")
            ), collapse = "_"), "-|\\/", "_"),
            "_",
            str_replace_all(ct, " |\\/", "_"),
            "_",
            organx,
            "_",
            expx,
            ".pdf"
          )
        ))
        
        fitx <- lmFit(voomx, design = designmat)
        fitx <- eBayes(fitx)
        
        for (coefx in colnames(designmat)) {
          topx <- topTable(fitx, number = Inf, coef = coefx)
          topx$rn <- rownames(topx)
          results.DGE[[ct]][[organx]][[compx]][[expx]][[coefx]] <-
            data.table(topx)
          
        }
        
        rm(
          list = c(
            "subsetx",
            "fitx",
            "designmat",
            "frame",
            "voomx",
            "countsx",
            "gfilter",
            "topx"
          )
        )
        gc(reset = T)
        
      }
      
    }
    
  }
  
}


### Res ---------------------------------------------------------------------
res <-
    rbindlist(lapply(results.DGE, function(l1) {
      rbindlist(lapply(l1, function(l2) {
        rbindlist(lapply(l2, function(l3) {
          rbindlist(lapply(l3, function(l4) {
            rbindlist(l4, idcol = "coef")}), idcol = "experiment")
            }), idcol = "treatment")
        }), idcol = "organ")
    }), idcol = "celltype")

saveRDS(res, file.path(dataDir, "DGE_raw.rds"))

res$treatment <-
  case_when(
    res$treatment %in% "NoT;HDAC_WT" ~ "WT_vs_ctrl",
    res$treatment %in% "NoT;HDAC_cKO" ~ "KO_vs_ctrl",
    res$treatment %in% "HDAC_WT;HDAC_cKO" ~ "KO_vs_WT",
    TRUE ~ res$treatment
  )

res$P.Value <- 
  case_when(
    res$P.Value == 0 ~ 5e-324,
    TRUE ~ res$P.Value
  )

res$adj.P.Val <- 
  case_when(
    res$adj.P.Val == 0 ~ 5e-324,
    TRUE ~ res$adj.P.Val
  )

res$celltype <- res$celltype %>% 
  str_replace_all(" |\\/", "_") %>% 
  str_remove(",")

res$direction <- case_when(res$logFC < 0 ~ "down",
                           res$logFC > 0 ~ "up",
                           TRUE ~ "NS")
res$direction2 <- case_when(res$logFC < -1 & res$adj.P.Val < 0.05 ~ "down",
                            res$logFC > 1  & res$adj.P.Val < 0.05 ~ "up",
                            TRUE ~ "NS")
res$direction3 <- case_when(res$logFC < -1 & res$adj.P.Val < 0.05 ~ "down.adj",
                            res$logFC < -1 & res$P.Value < 0.05 & res$adj.P.Val >= 0.05 ~ "down",
                            res$logFC > 1 & res$adj.P.Val < 0.05 ~ "up.adj",
                            res$logFC > 1 & res$P.Value < 0.05 & res$adj.P.Val >= 0.05 ~ "up",
                               TRUE ~ "NS")

saveRDS(res, file.path(dataDir, "DGE_res.rds"))

for (ct in unique(res$celltype)){
  for (treat in unique(res$treatment)){
    for (organx in unique(res$organ)){
      for (experimentx in unique(res$experiment)){
        
        sx <- res %>%
          filter(
            celltype == ct,
            treatment == treat,
            organ == organx,
            experiment == experimentx,
            adj.P.Val < 0.05,
            !coef %in% c("M", "Intercept")
          ) 
        
        if (nrow(sx) == 0) next
        
        sx %>%
          write.csv(file.path(tablesDir, paste0(
            paste("gene_list", ct, organx, experimentx, treat,
                  sep = "_"),
            ".csv"
          )), row.names = F)
        
      }
    }
  }
}




### Vulcano & Histo ---------------------------------------------------------

res <- readRDS(file.path(dataDir, "DGE_res.rds"))

for(ct in unique(res$celltype)) {
  for (organx in unique(res$organ)) {
    for (expx in unique(res$experiment)) {
      sx <- res %>% filter(celltype == ct,
                     organ == organx,
                     experiment == expx,
                     !grepl("Intercept|M", coef)) 
      
      if (nrow(sx) == 0) next
      
      xmax <- max(abs(sx$logFC)) + 1
      xmin <- -xmax
      ymax <- sx %>% filter(direction3 %in% c("up","down")) %>% group_by(treatment) %>% summarise(x = max(-log10(P.Value)))
      
      
      sx %>% 
        ggplot(aes(logFC, -log10(P.Value))) +
        geom_point(aes(fill = direction2),
                   shape = 21,
                   alpha = 0.5) +
        facet_wrap(~ treatment, scales = "free_y")+
        # geom_ribbon(data = sx %>% filter(direction3 %in% c("up","down")),
        #             aes(y = -log10(P.Value),
        #                 x = logFC,
        #                 xmin = xmin,
        #                 xmax = xmax,
        #                 group = treatment),
        #             alpha = 0.2)  +
        scale_x_continuous(limits = c(xmin - 0.1, xmax + 0.1), expand = expansion(0,0)) + 
        theme_my() +
        theme(
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(
            colour = "black",
            fill = NA,
            size = 1
          ),
          axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.key = element_blank()
        ) +
        scale_fill_manual(
          values = c("#C8466D", "#3F60AE","#E6E6E6"),
          name = "",
          limits = c("up", "down", "NS"),
          labels = c("Up", "Down", "Unaffected")
        ) +
        geom_vline(
          xintercept = 1,
          color = "black",
          linetype = "dotted",
          size = 0.5
        ) +
        geom_vline(
          xintercept = -1,
          color = "black",
          linetype = "dotted",
          size = 0.5
        ) +
        geom_hline(data = ymax,
                   aes(yintercept = x),
                   color = "black",
                   linetype = "dotted",
                   size = 0.5
        ) +
        ggtitle(paste("vulcano", ct, organx, expx, sep = "_"))
      
      ggsave(file.path(plotsDir,
                       paste0(
                         paste("vulcano", ct, organx, expx, sep = "_"),
                         ".pdf"
                       )))
      
      sx %>% 
        ggplot() +
        geom_histogram(aes(P.Value, fill = factor(floor(AveExpr))),
                       bins = 30
        ) +
        facet_wrap(~ treatment) +
        theme_my() +
        theme(panel.grid.major = element_blank())
      
      ggsave(file.path(plotsDir,
                       paste0(
                         paste("histo", ct, organx, expx, sep = "_"),
                         ".pdf"
                       )))
    }
  }
}

# ### Heatmap -----------------------------------------------------------------
# scaled_counts <- t(scale(t(voomx$E)))
# 
# 
# subsetx <- monocle.obj[,
#                        colData(monocle.obj) %>%
#                          subset(
#                            celltype == ct &
#                              organ == organx &
#                              treatment.agg %in% unlist(str_split(compx, ";")) &
#                              experiment == expx
#                          ) %>%
#                          row.names]
# 
# Heatmap(scaled_counts)

## Pseudo Bulk -------------------------------------------------------------

### Create Pseudo Bulk ------------------------------------------------------

countsx <- as(matrix(nrow = nrow(monocle.obj)), "sparseMatrix")
designx <- DataFrame()

for(samplex in unique(monocle.obj@colData$sample)) {
  for (ct in unique(monocle.obj@colData$celltype)) {
    for (treat in c("NoT", "HDAC_WT", "HDAC_cKO")) {
      
      subsetx <- monocle.obj[, monocle.obj@colData %>% subset(celltype == ct &
                                                                treatment.agg == treat &
                                                                sample == samplex) %>% rownames()]
      
      if(ncol(subsetx) == 0) next
      
      x <-  counts(subsetx) %>% rowSums() %>% as("sparseMatrix")
      
      
      countsx <- cbind(countsx, x)
      
      designx <- rbind(designx,
                       DataFrame(
                         subsetx@colData %>% 
                           as_tibble() %>% 
                           select(celltype, treatment.agg, sample, sex, organ, experiment) %>%  
                           unique()
                       )
      )
      
      colnames(countsx)[ncol(countsx)] <-
        paste(samplex, str_replace_all(ct, " |\\/", "_"), treat, sep = "_")
      
    }
  }
}

countsx <- countsx[,-1]
row.names(designx) <- colnames(countsx)

designx$treatment07 <- forcats::fct_relevel(designx$treatment07, c("NoT", "HDAC_WT", "HDAC_cKO", "Undefined"))
colnames(designx) <- colnames(designx) %>% str_remove("\\d+")

saveRDS(countsx, file = file.path("data/07", "counts_pseudobulk.rds"))
saveRDS(designx, file = file.path("data/07", "design_pseudobulk.rds"))

### Limma -------------------------------------------------------------------
# comp <- combn(unique(levels(
#         monocle.obj[, colData(monocle.obj) %>%
#                             subset(treatment != "Undefined") %>% row.names]@colData$treatment %>%
#                 droplevels
# )),
# 2,
# FUN = paste,
# collapse = ";")
gene_filter_list <- list()
for(pwd in c("data/agg", "data/55", "data/07", "data/08")) {

countsx <- readRDS(file = file.path("../Desktop/Shiny_Backup_data/",pwd, "counts_pseudobulk.rds"))
designx <- readRDS(file = file.path("../Desktop/Shiny_App/",pwd, "design_pseudobulk.rds"))

message(pwd,"limma start ---------------")
results.DGE.pseudo <- list()

for (ct in unique(designx$celltype)) {
        message(ct, "-------------------------")
  for (organx in unique(designx$organ)) {
      for (expx in unique(designx$experiment)) {
              
        subsetx <- designx %>% 
          subset(
            celltype == ct &
            organ == organx &
            experiment == expx
        )
        
        if(nrow(subsetx) < 2) next
        
        countsx.pseudo <- countsx[,rownames(subsetx)]
        
        if (length(unique(subsetx$treatment)) < 2) next
        
        if (length(unique(subsetx$sex)) < 2) {
                frame <-
                        model.frame(~ 0 + treatment,
                                    subsetx, drop.unused.levels = T)
                designmat <- model.matrix(~ 0 + treatment, frame)
                
        } else {
                frame <-
                        model.frame(~ 0 + treatment + sex,
                                    subsetx, drop.unused.levels = T)
                designmat <-
                        model.matrix(~ 0 + treatment + sex, frame)
                
        }
        
        if (ncol(designmat) >= nrow(designmat)) {
                message(pwd, "_", ct, "_", organx, "_", expx,  "---ncol(designmat) >= nrow(designmat)---")
                next
        }
        
        colnames(designmat) <-
                colnames(designmat) %>%
                str_remove_all("\\(|\\)|sample|treatment|\\/|[0-9]\\/[0-9]|sex") %>%
                str_replace_all("\\-", "_")

        dge <- DGEList(countsx.pseudo)
        
        gfilter <-
                filterByExpr(
                        dge,
                        design = designmat,
                        large.n = 10,
                        min.count = 7,
                        min.prop = 0.7,
                        min.total.count = 8
                )
        
        gfilter2 <- rowSums(dge$counts) > 0
        
        # df <- data.frame(strict_filter = sum(gfilter), old_filter = sum(gfilter2))
        
        # gene_filter_list[[ct]][[organx]][[expx]][[pwd]] <- df
        
      # }}}}

        
        dge <- dge[gfilter2,,keep.lib.sizes=FALSE]

        dge <- calcNormFactors(dge)

        voomx <- voomWithQualityWeights(dge,
                                        design = designmat,
                                        plot = F,
                                        save.plot = F)
        
        contrast.matrix <- makeContrasts(contrasts = combn(rev(colnames(designmat)[colnames(designmat) != "M"]),
                                                           2,
                                                           FUN = paste,
                                                           collapse = "-"),
                                         levels = designmat)
        
        fitx <- lmFit(voomx, design = designmat)
        fitx <- contrasts.fit(fitx, contrast.matrix)
        fitx <- eBayes(fitx)

        write.fst(data.table(voomx$E, keep.rownames = T), path = file.path(pwd,paste0("voom_counts_", ct, "_", organx, "_", expx, ".fst")), compress = 50)
        
        for (coefx in colnames(contrast.matrix)) {
          topx <- topTable(fitx, number = Inf, coef = coefx)
          topx$rn <- rownames(topx)
          results.DGE.pseudo[[ct]][[organx]][[expx]][[coefx]] <-
            data.table(topx)
        }
        
      }
      
    }
    
  }
  
# genes <-
#         res.pb <-
#         rbindlist(lapply(gene_filter_list, function(l1) {
#                 rbindlist(lapply(l1, function(l2) {
#                         rbindlist(lapply(l2, function(l3) {
#                                 rbindlist(l3, idcol = "strictness")}), 
#                         idcol = "experiment")
#                 }), idcol = "organ")
#         }), idcol = "celltype")
# 
# genes$strictness <- genes$strictness %>% str_remove("\\w+\\/")
# 
# genes <- genes %>% pivot_longer(cols = contains("filter"), names_to = "filter_type", values_to = "n_genes")
# 
# genes %>% ggplot + 
#         geom_col(aes(celltype, n_genes, fill = filter_type), position = position_dodge()) +
#         facet_wrap( ~ organ + experiment + strictness, ncol = 4) +
#         theme_my()  
# 
# ggsave("gfilter.pdf", scale = 1.5, height = 7.5, width = 15)

### Res ---------------------------------------------------------------------
message(pwd,"limma end ---------------")
res.pb <-
  rbindlist(lapply(results.DGE.pseudo, function(l1) {
    rbindlist(lapply(l1, function(l2) {
        rbindlist(lapply(l2, function(l3) {
          rbindlist(l3, idcol = "coef")}), 
          idcol = "experiment")
    }), idcol = "organ")
  }), idcol = "celltype")

# saveRDS(res.pb, file.path(dataDir, "DGE_raw_pb.rds"))

res.pb$treatment <-
  case_when(
    res.pb$coef %in% "HDAC_WT-NoT" ~ "WT_vs_ctrl",
    res.pb$coef %in% "HDAC_cKO-NoT" ~ "KO_vs_ctrl",
    res.pb$coef %in% "HDAC_cKO-HDAC_WT" ~ "KO_vs_WT",
    TRUE ~ res.pb$coef
  )

# res.pb$P.Value <- 
#   case_when(
#     res.pb$P.Value == 0 ~ 5e-324,
#     TRUE ~ res.pb$P.Value
#   )
# 
# res.pb$adj.P.Val <- 
#   case_when(
#     res.pb$adj.P.Val == 0 ~ 5e-324,
#     TRUE ~ res.pb$adj.P.Val
#   )

res.pb$celltype <- res.pb$celltype %>% 
  str_replace_all(" |\\/\\-", "_") %>% 
  str_remove(",")

res.pb$direction <- case_when(res.pb$logFC < 0 ~ "down",
                           res.pb$logFC > 0 ~ "up",
                           TRUE ~ "NS")
res.pb$direction2 <- case_when(res.pb$logFC < -1 & res.pb$adj.P.Val < 0.05 ~ "down",
                            res.pb$logFC > 1  & res.pb$adj.P.Val < 0.05 ~ "up",
                            TRUE ~ "NS")

res.pb$direction3 <- case_when(res.pb$logFC < -1 & res.pb$adj.P.Val < 0.05 ~ "down.adj",
                               res.pb$logFC < -1 & res.pb$P.Value < 0.05 & res.pb$adj.P.Val >= 0.05 ~ "down",
                               res.pb$logFC > 1 & res.pb$adj.P.Val < 0.05 ~ "up.adj",
                               res.pb$logFC > 1 & res.pb$P.Value < 0.05 & res.pb$adj.P.Val >= 0.05 ~ "up",
                               TRUE ~ "NS")



# saveRDS(res.pb, file.path("old_res.rds"))

# shiny_list <- res.pb %>% filter(!grepl("M|Inter", coef))

message(pwd,"creating files ---------------")
x <- res.pb
x$celltype <- x$celltype %>% str_replace("\\-","_")

for(ct in unique(x$celltype)){
        for(organx in unique(x$organ)){
                for(expx in unique(x$experiment)){

                        subsetx <- x %>% filter(celltype == ct,
                                                organ == organx,
                                                experiment == expx)

                        subsetx$start <- case_when(subsetx$treatment %>% str_ends("ctrl") ~ "NoT",
                                                   subsetx$treatment %>% str_ends("WT") ~ "WT")

                        subsetx$end <- case_when(subsetx$treatment %>% str_starts("WT") ~ "WT",
                                                 subsetx$treatment %>% str_starts("KO") ~ "cKO")

                        subsetx$treatment <- factor(subsetx$treatment, levels = c("WT_vs_ctrl", "KO_vs_ctrl", "KO_vs_WT"))

                        subsetx <- subsetx %>% select(-c(direction, B, t, AveExpr, coef, direction2, direction3))

                        write.fst(subsetx,path = file.path(pwd, paste0("list_",ct,"_",organx,"_",expx,".fst")), compress = 50)
                }
        }
}

x <- x %>% select(rn, celltype, organ, experiment, direction2, treatment, t)
x$treatment <- factor(x$treatment, levels = c("WT_vs_ctrl", "KO_vs_ctrl", "KO_vs_WT"))
write.fst(x, file.path(pwd,"list_for_shiny.fst"), compress = 100)

# saveRDS(res.pb, file.path(pwd,"list_for_shiny.rds"))
}


for (ct in unique(res.pb$celltype)){
  for (treat in unique(res.pb$treatment)){
    for (organx in unique(res.pb$organ)){
      for (experimentx in unique(res.pb$experiment)){
        
        sx <- res.pb %>%
          filter(
            celltype == ct,
            treatment == treat,
            organ == organx,
            experiment == experimentx,
            adj.P.Val < 0.05,
            !coef %in% c("M", "Intercept")
          ) 
        
        if (nrow(sx) == 0) next
        
        sx %>%
          write.csv(file.path(tablesDir, paste0(
            paste("gene_list", ct, organx, experimentx, treat, "pseudobulk",
                  sep = "_"),
            ".csv"
          )), row.names = F)
        
      }
    }
  }
}

### cors --------------------------------------------------------------------


for(organx in unique(designx$organ)){
  for(ct in unique(designx$celltype)){
    subsetx <- designx %>% subset(celltype == ct & organ == organx & experiment == "HDAC1")
    countsx.pseudo <- countsx[,rownames(subsetx)]
    corx <- cor(as.matrix(countsx.pseudo))
    corx <- corx %>% as_tibble(rownames = "rn") %>% tibble::column_to_rownames("rn")
    
    not <- corx %>% select(contains("NoT"))
    
    avg_not <- ifelse(grepl("NoT", rownames(not)),
                      (rowSums(not) - 1) / (ncol(not) - 1),
                      rowSums(not) / ncol(not)) %>% setNames(rownames(not))
    
    wt <- corx %>% select(contains("HDAC_WT"))
    
    avg_wt <- ifelse(grepl("HDAC_WT", rownames(wt)),
                     (rowSums(wt) - 1) / (ncol(wt) - 1),
                     rowSums(wt) / ncol(wt)) %>% setNames(rownames(wt))
    
    ko <- corx %>% select(contains("HDAC_cKO"))
    
    avg_ko <- ifelse(grepl("HDAC_cKO", rownames(ko)),
                     (rowSums(ko) - 1) / (ncol(ko) - 1),
                     rowSums(ko) / ncol(ko)) %>% setNames(rownames(ko))
    
    if(all(is.na(c(avg_ko,avg_wt, avg_not)))) next
    
    result <- data.frame(avg_not,avg_wt, avg_ko) %>% 
      as_tibble(rownames = "rn") %>% 
      mutate(soi = ifelse(grepl("f.*40B3", rn),"yes","no"))
    
    result <- result %>% 
      pivot_longer(2:4, names_to = "treat", values_to = "cor") %>% 
      mutate(type = str_extract(rn, "[:alpha:]+$")) %>%  
      arrange(treat,type)
    
    result$rn <- forcats::fct_inorder(result$rn)
    
    result %>% 
      ggplot()+
      geom_col(aes(rn, cor, fill = soi),
               show.legend = F) + 
      theme_my() +
      facet_wrap(~treat) +
      scale_fill_manual(values = c("yes"="red"))
    
    
    ggsave(file.path(plotsDir,paste0("avg_cor_",
                                     str_replace(ct," ","_"),
                                     "_",
                                     organx,
                                     "_HDAC1.pdf")),
           width = 9)
  }
}


### Vulcano & Histo ---------------------------------------------------------

res.pb <- readRDS(file.path(dataDir, "DGE_res_pb.rds"))


for(ct in unique(res.pb$celltype)) {
  for (organx in unique(res.pb$organ)) {
    for (expx in unique(res.pb$experiment)) {
      subsetx <- res.pb %>% filter(celltype == ct,
                                   organ == organx,
                                   experiment == expx,
                                   !grepl("Intercept|M", coef)) 
      if (nrow(subsetx) == 0) next

      xmax <- max(abs(subsetx$logFC)) + 1
      xmin <- -xmax
      ymax <-
        subsetx %>% filter(direction3 %in% c("up", "down")) %>% group_by(treatment) %>% summarise(x = max(-log10(P.Value)))
      
      subsetx %>% 
        ggplot(aes(logFC, -log10(P.Value))) +
        geom_point(aes(fill = direction2),
                   shape = 21,
                   alpha = 0.5) +
        facet_wrap(~ treatment, scales = "free_y")+
        # geom_ribbon(data = subsetx %>% filter(direction3 %in% c("up","down")),
        #             aes(y = -log10(P.Value),
        #                 x = logFC,
        #                 xmin = xmin,
        #                 xmax = xmax,
        #                 group = treatment),
        #             alpha = 0.2)  +
        scale_x_continuous(limits = c(xmin - 0.1, xmax + 0.1), expand = expansion(0,0)) + 
        theme_my() +
        theme(
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(
            colour = "black",
            fill = NA,
            size = 1
          ),
          axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.key = element_blank()
        ) +
        scale_fill_manual(
          values = c("#C8466D", "#3F60AE","#E6E6E6"),
          name = "",
          limits = c("up", "down", "NS"),
          labels = c("Up", "Down", "Unaffected")
        ) +
        geom_vline(
          xintercept = 1,
          color = "black",
          linetype = "dotted",
          size = 0.5
        ) +
        geom_vline(
          xintercept = -1,
          color = "black",
          linetype = "dotted",
          size = 0.5
        )  +
        geom_hline(data = ymax,
                   aes(yintercept = x),
                   color = "black",
                   linetype = "dotted",
                   size = 0.5
        ) +
        ggtitle(paste("pseudo_vulcano", ct, organx, expx, sep = "_"))
      
      ggsave(file.path(plotsDir,
                       paste0(
                         paste("pseudo_vulcano", ct, organx, expx, sep = "_"),
                         ".pdf"
                       )))
      
      subsetx %>% 
        ggplot() +
        geom_histogram(aes(P.Value, fill = factor(floor(AveExpr))),
                       bins = 30
        ) +
        facet_wrap(~ treatment) +
        theme_my() +
        theme(panel.grid.major = element_blank()) +
        ggtitle(paste("pseudo_histo", ct, organx, expx, sep = "_"))
      
      ggsave(file.path(plotsDir,
                       paste0(
                         paste("pseudo_histo", ct, organx, expx, sep = "_"),
                         ".pdf"
                       )))
      
    }
  }
}


## SC Pseudobulk Cor -------------------------------------------------------
res <- readRDS(file.path(dataDir, "DGE_res.rds"))
res.pb <- readRDS(file.path(dataDir, "DGE_res_pb.rds"))


for (ct in unique(res$celltype)){
  for (organx in unique(res$organ)){
    for (expx in unique(res$experiment)){
      
      sx <- res %>%
        filter(
          !coef %in% c("M", "Intercept"),
          organ == organx,
          experiment == expx,
          celltype == ct
          ) %>%
        select(rn, logFC, treatment)
      
      sx.pb <- res.pb %>% 
        filter(
          !coef %in% c("M", "Intercept"),
          organ == organx,
          experiment == expx,
          celltype == ct,
          rn %in% sx$rn
        ) %>%
        select(rn, logFC, treatment)
      
      sx <- sx %>% filter(rn %in% sx.pb$rn)
      
      if (any(nrow(sx) == 0, nrow(sx.pb) == 0)) next
      
      df <- inner_join(sx, sx.pb, by = c("rn", "treatment"))
      
      colnames(df) <-
        case_when(
          grepl("\\.x", colnames(df)) ~ str_replace(colnames(df), "\\.x", "_sc"),
          grepl("\\.y", colnames(df)) ~ str_replace(colnames(df), "\\.y", "_pb"),
          TRUE ~ colnames(df)
        )
      
      cors <- list()
      for(tx in unique(df$treatment)){
        dfx <- df %>% filter(treatment == tx)
        cors[[tx]] <- cor(dfx$logFC_sc, dfx$logFC_pb, method = "spearman")
      } 
      cors <- bind_rows(cors, .id = "treatment")

      cors <- cors %>% pivot_longer(cols = 1:ncol(cors), values_to = "Correlation")

      p1 <- 
        df %>%
        ggplot +
        stat_binhex(aes(logFC_sc,
                        logFC_pb,
                        fill = log10(..count..)),
                    bins = 30,
                    col = "grey15") +
        scale_fill_gradient(low = "ivory2",
                            high = "#D35FB7") +
        facet_wrap(~ treatment,
                   scales = "free",
                   nrow = 2) +
        geom_smooth(
          aes(logFC_sc,
              logFC_pb),
          method = "lm",
          se = T,
          col = "black"
        ) +
        theme_my() 
      
      p2 <- 
        cors %>% 
        mutate(x = "Correlation") %>% 
        ggplot() + 
        geom_col(aes(name, 
                     Correlation),
                 col = "black",
                 fill = "deepskyblue1") +
        facet_wrap(~x) +
        theme_my() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
              axis.title = element_blank())
      
      if (nrow(cors) == 3){
      p1 + inset_element(p2, left = 0.55, bottom = 0, right = 0.95, top = 0.40, on_top = F)
      }
      else p1 + p2
      ggsave(file.path(plotsDir, paste0("cor_sc_pb_",
                                       ct,
                                       "_",
                                       organx,
                                       "_",
                                       expx,
                                       ".pdf")))
    }
  }
}



## FB & KC Clusters --------------------------------------------------------

### Limma -------------------------------------------------------------------

results.MA <- list()
for (ct in c("Keratinocytes", "Fibroblasts")) {
  subsetx <- monocle.obj[,
                         colData(monocle.obj) %>%
                           subset(celltype == ct &
                                    organ == "Skin" &
                                    treatment.agg %in% c("NoT", "HDAC_WT")) %>%
                           row.names]
  for (clusterx in unique(subsetx@colData$Cluster)) {
    
    sx <- subsetx[, colData(subsetx) %>% subset(Cluster == clusterx) %>% rownames()]
    
    if (length(unique(sx@colData$treatment.agg)) >= 2 &
        length(unique(sx@colData$sex)) >= 2) {
      
      frame <-
        model.frame( ~ treatment.agg + sex,
                     sx@colData %>% droplevels())
      designmat <- model.matrix( ~ treatment.agg + sex, frame)
      
      if (ncol(designmat) >= nrow(designmat)) next
      
    } else if (length(unique(sx@colData$treatment.agg)) >= 2 &
               length(unique(sx@colData$sex)) < 2) {
      frame <-
        model.frame( ~ treatment.agg, sx@colData %>% droplevels())
      designmat <- model.matrix( ~ treatment.agg, frame)
      
      if (ncol(designmat) >= nrow(designmat)) next
      
      
    } else next
    
    colnames(designmat) <-
      colnames(designmat) %>%
      str_remove_all("\\(|\\)|sample|treatment.agg|\\/|[0-9]\\/[0-9]|sex") %>%
      str_replace_all("\\-", "_")
    
    
    countsx <- counts(sx)
    
    # gfilter <-
    #   filterByExpr(
    #     countsx,
    #     design = designmat,
    #     min.prop = 0.1,
    #     min.count = 1
    #   )
    
    gfilter <- rowSums(countsx) > 0
    
    voomx <-
      voom(countsx[gfilter, ],
           design = designmat,
           plot = F,
           save.plot = T)
    
    ggplot() +
      geom_point(aes(voomx$voom.xy$x, voomx$voom.xy$y), size = 0.1) +
      geom_line(aes(voomx$voom.line$x, voomx$voom.line$y), col = "red") +
      xlab(voomx$voom.xy$xlab) +
      ylab(voomx$voom.xy$ylab) +
      ggtitle("voom: Mean-variance trend") +
      theme_my()
    
    
    ggsave(file.path(
      plotsDir,
      paste0(
        "MA_voom_",
        str_replace_all(ct, " |\\/", "_"),
        "_",
        clusterx,
        ".pdf"
      )
    ))
    
    fitx <- lmFit(voomx, design = designmat)
    fitx <- eBayes(fitx)
    
    for (coefx in colnames(designmat)) {
      topx <- topTable(fitx, number = Inf, coef = coefx)
      topx$rn <- rownames(topx)
      results.MA[[ct]][[clusterx]][[coefx]] <-
        data.table(topx)
      
    }
    
  }
  
}


### Res ---------------------------------------------------------------------

res.MA <-
  rbindlist(lapply(results.MA, function(l1) {
    rbindlist(lapply(l1, function(l2) {
      rbindlist(l2, idcol = "coef")
    }), idcol = "cluster")
  }), idcol = "celltype")

saveRDS(res.MA, file.path(dataDir, "DGE_raw_MA.rds"))



res.MA$P.Value <- 
  case_when(
    res.MA$P.Value == 0 ~ 5e-324,
    TRUE ~ res.MA$P.Value
  )

res.MA$adj.P.Val <- 
  case_when(
    res.MA$adj.P.Val == 0 ~ 5e-324,
    TRUE ~ res.MA$adj.P.Val
  )

res.MA$celltype <- res.MA$celltype %>% 
  str_replace_all(" |\\/", "_") %>% 
  str_remove(",")

res.MA$direction <- case_when(res.MA$logFC < 0 ~ "down",
                              res.MA$logFC > 0 ~ "up",
                              TRUE ~ "NS")
res.MA$direction2 <- case_when(res.MA$logFC < -1 & res.MA$adj.P.Val < 0.05 ~ "down",
                               res.MA$logFC > 1  & res.MA$adj.P.Val < 0.05 ~ "up",
                               TRUE ~ "NS")
res.MA$direction3 <- case_when(res.MA$logFC < -1 & res.MA$adj.P.Val < 0.05 ~ "down.adj",
                            res.MA$logFC < -1 & res.MA$P.Value < 0.05 & res.MA$adj.P.Val >= 0.05 ~ "down",
                            res.MA$logFC > 1 & res.MA$adj.P.Val < 0.05 ~ "up.adj",
                            res.MA$logFC > 1 & res.MA$P.Value < 0.05 & res.MA$adj.P.Val >= 0.05 ~ "up",
                            TRUE ~ "NS")

saveRDS(res.MA, file.path(dataDir, "DGE_res_MA.rds"))



### Vulcano & Histo ---------------------------------------------------------

res.MA <- readRDS(file.path(dataDir, "DGE_res_MA.rds"))

for (ct in unique(res.MA$celltype)) {
  for (clusterx in unique(res.MA$cluster)) {
    subsetx <- res.MA %>% filter(celltype == ct,
                                 cluster == clusterx,
                                 !grepl("Intercept|M", coef))
    
    if (nrow(subsetx) == 0) next
    
    xmax <- max(abs(subsetx$logFC)) + 1
    xmin <- -xmax
    
    subsetx %>%
      ggplot(aes(logFC,-log10(P.Value))) +
      geom_point(aes(fill = direction2),
                 shape = 21,
                 alpha = 0.5) +
      geom_ribbon(data = subsetx %>% filter(direction3 %in% c("up","down")),
                  aes(y = -log10(P.Value),
                      x = logFC,
                      xmin = xmin,
                      xmax = xmax),
                  alpha = 0.2)  +
      scale_x_continuous(limits = c(xmin - 0.1, xmax + 0.1), expand = expansion(0,0)) + 
      #facet_wrap( ~ organ, scales = "free") +
      theme_my() +
      theme(panel.grid.major = element_blank()) +
      scale_fill_manual(
        values = c("#C8466D", "#3F60AE", "#E6E6E6"),
        name = "",
        limits = c("up", "down", "NS"),
        labels = c("Up", "Down", "Unaffected")
      ) +
      geom_vline(
        xintercept = 1,
        color = "black",
        linetype = "dotted",
        size = 0.5
      ) +
      geom_vline(
        xintercept = -1,
        color = "black",
        linetype = "dotted",
        size = 0.5
      ) +
      ggtitle(paste("MA_vulcano", ct,"cluster", clusterx, sep = "_"))
    
    ggsave(file.path(plotsDir,
                     paste0(
                       paste("MA_vulcano", ct,"cluster", clusterx, sep = "_"),
                       ".pdf"
                     )))
    
    subsetx %>% 
      ggplot() +
      geom_histogram(aes(P.Value, fill = factor(floor(AveExpr))),
                     bins = 30
      ) +
      theme_my() +
      theme(panel.grid.major = element_blank()) +
      ggtitle(paste("MA_histo", ct,"cluster", clusterx, sep = "_"))
    
    ggsave(file.path(plotsDir,
                     paste0(
                       paste("MA_histo", ct,"cluster", clusterx, sep = "_"),
                       ".pdf"
                     )))
  }
}
## Clusters Pseudo Bulk -------------------------------------------------------------

### Create Pseudo Bulk ------------------------------------------------------

pb_counts_MA <- as(matrix(nrow = nrow(monocle.obj)), "sparseMatrix")
pb_design_MA <- DataFrame()

pb_object <- monocle.obj[,
                         colData(monocle.obj) %>%
                           subset(celltype %in% c("Keratinocytes", "Fibroblasts") &
                                    organ == "Skin" &
                                    treatment.agg %in% c("NoT", "HDAC_WT") &
                                    experiment == "HDAC1") %>%
                           row.names]

for(samplex in unique(pb_object@colData$sample)) {
  for (ct in unique(pb_object@colData$celltype)) {
    for (treat in unique(pb_object@colData$treatment.agg[pb_object@colData$treatment.agg != "Undefined"])) {
      subsetx <-
        pb_object[, pb_object@colData %>% subset(celltype == ct &
                                                       treatment.agg == treat &
                                                       sample == samplex) %>% rownames()]
      if(ncol(subsetx) == 0) next
      
      for (clusterx in unique(subsetx@colData$Cluster)) {
        sx <-
          subsetx[, subsetx@colData %>% subset(Cluster == clusterx) %>% rownames()]
        
        if (ncol(sx) == 0) next
        
        x <-  counts(sx) %>% rowSums() %>% as("sparseMatrix")
        
        
        pb_counts_MA <- cbind(pb_counts_MA, x)
        
        pb_design_MA <- rbind(pb_design_MA,
                              DataFrame(
                                sx@colData %>%
                                  as_tibble() %>%
                                  select(celltype, treatment.agg, sample, sex, organ, experiment, Cluster) %>%
                                  unique()
                              ))
        
        colnames(pb_counts_MA)[ncol(pb_counts_MA)] <-
          paste("Cluster",
                clusterx,
                samplex,
                str_replace_all(ct, " |\\/", "_"),
                treat,
                sep = "_")
        
      }
    }
  }
}
pb_counts_MA <- pb_counts_MA[,-1]
row.names(pb_design_MA) <- colnames(pb_counts_MA)

saveRDS(pb_counts_MA, file = file.path(dataDir, "MA_counts_pseudobulk.rds"))
saveRDS(pb_design_MA, file = file.path(dataDir, "MA_design_pseudobulk.rds"))

### Limma -------------------------------------------------------------------

pb_counts_MA <- readRDS(file = file.path(dataDir, "MA_counts_pseudobulk.rds"))
pb_design_MA <- readRDS(file = file.path(dataDir, "MA_design_pseudobulk.rds"))


results.DGE.MA.pseudo <- list()
for (ct in unique(pb_design_MA$celltype)) {
        subsetx <- pb_design_MA %>%
                        subset(celltype == ct &
                                 organ == "Skin" &
                                 treatment.agg %in% c("NoT", "HDAC_WT"))
        
        for (clusterx in unique(subsetx$Cluster)) {
          
          sx <- subsetx %>% subset(Cluster == clusterx) 
          
          if (length(unique(sx$treatment.agg)) >= 2 &
              length(unique(sx$sex)) >= 2) {
            
            frame <-
              model.frame( ~ treatment.agg + sex,
                           sx %>% droplevels())
            designmat <- model.matrix( ~ treatment.agg + sex, frame)
            
            if (ncol(designmat) >= nrow(designmat)) next
            
          } else if (length(unique(sx$treatment.agg)) >= 2 &
                     length(unique(sx$sex)) < 2) {
            frame <-
              model.frame( ~ treatment.agg, sx %>% droplevels())
            designmat <- model.matrix( ~ treatment.agg, frame)
            
            if (ncol(designmat) >= nrow(designmat)) next
            
            
          } else next
          
          colnames(designmat) <-
            colnames(designmat) %>%
            str_remove_all("\\(|\\)|sample|treatment.agg|\\/|[0-9]\\/[0-9]|sex") %>%
            str_replace_all("\\-", "_")
          
          
          countsx <- pb_counts_MA[,rownames(sx)]
          
          # gfilter <-
          #   filterByExpr(
          #     countsx,
          #     design = designmat,
          #     min.prop = 0.1,
          #     min.count = 1
          #   )
          
          gfilter <- rowSums(countsx) > 0
          
          voomx <-
            voom(countsx[gfilter, ],
                 design = designmat,
                 plot = F,
                 save.plot = T)
          
          # ggplot() +
          #   geom_point(aes(voomx$voom.xy$x, voomx$voom.xy$y), size = 0.1) +
          #   geom_line(aes(voomx$voom.line$x, voomx$voom.line$y), col = "red") +
          #   xlab(voomx$voom.xy$xlab) +
          #   ylab(voomx$voom.xy$ylab) +
          #   ggtitle("voom: Mean-variance trend") +
          #   theme_my()
          # 
          # 
          # ggsave(file.path(
          #   plotsDir,
          #   paste0(
          #     "MA_pseudo_voom_",
          #     clusterx,
          #     str_replace_all(ct, " |\\/", "_"),
          #     "_",
          #     ".pdf"
          #   )
          # ))
          
          fitx <- lmFit(voomx, design = designmat)
          fitx <- eBayes(fitx)
          
          for (coefx in colnames(designmat)) {
            topx <- topTable(fitx, number = Inf, coef = coefx)
            topx$rn <- rownames(topx)
            results.DGE.MA.pseudo[[ct]][[clusterx]][[coefx]] <-
              data.table(topx)
            
          }
          
        }
        
}
### Res ---------------------------------------------------------------------

res.pb.MA <-
  rbindlist(lapply(results.DGE.MA.pseudo, function(l1) {
    rbindlist(lapply(l1, function(l2) {
      rbindlist(l2, idcol = "coef")
    }), idcol = "cluster")
  }), idcol = "celltype")

saveRDS(res.pb.MA, file.path(dataDir, "DGE_raw_MA_pb.rds"))

res.pb.MA$P.Value <- 
  case_when(
    res.pb.MA$P.Value == 0 ~ 5e-324,
    TRUE ~ res.pb.MA$P.Value
  )

res.pb.MA$adj.P.Val <- 
  case_when(
    res.pb.MA$adj.P.Val == 0 ~ 5e-324,
    TRUE ~ res.pb.MA$adj.P.Val
  )

res.pb.MA$celltype <- res.pb.MA$celltype %>% 
  str_replace_all(" |\\/", "_") %>% 
  str_remove(",")

res.pb.MA$direction <- case_when(res.pb.MA$logFC < 0 ~ "down",
                              res.pb.MA$logFC > 0 ~ "up",
                              TRUE ~ "NS")
res.pb.MA$direction2 <- case_when(res.pb.MA$logFC < -1 & res.pb.MA$adj.P.Val < 0.05 ~ "down",
                               res.pb.MA$logFC > 1  & res.pb.MA$adj.P.Val < 0.05 ~ "up",
                               TRUE ~ "NS")

res.pb.MA$direction3 <- case_when(res.pb.MA$logFC < -1 & res.pb.MA$adj.P.Val < 0.05 ~ "down.adj",
                               res.pb.MA$logFC < -1 & res.pb.MA$P.Value < 0.05 & res.pb.MA$adj.P.Val >= 0.05 ~ "down",
                               res.pb.MA$logFC > 1 & res.pb.MA$adj.P.Val < 0.05 ~ "up.adj",
                               res.pb.MA$logFC > 1 & res.pb.MA$P.Value < 0.05 & res.pb.MA$adj.P.Val >= 0.05 ~ "up",
                               TRUE ~ "NS")



saveRDS(res.pb.MA, file.path(dataDir, "DGE_res_MA_pb.rds"))

### Vulcano & Histo ---------------------------------------------------------


# Markers -----------------------------------------------------------------

markers <- read.table(file.path(tablesDir, "PanglaoDB_markers_27_Mar_2020.tsv"), sep = "\t", header = T, quote = '')

markers <- markers %>% filter(
    grepl("T .*cells|B .*cells|Keratinocytes|Fibroblasts", cell.type),
    grepl("Mm", species)
) %>% mutate(gene_symbol = str_to_title(official.gene.symbol)) %>%
    dplyr::select(gene_symbol, organ, cell.type)

markers <- markers %>% add_row(gene_symbol = "Btla",organ = "x",cell.type = "T follicular helper cells")
markers <- markers %>% add_row(gene_symbol = "Il21",organ = "x",cell.type = "T follicular helper cells")
markers <- markers %>% add_row(gene_symbol = "Sh2d1a",organ = "x",cell.type = "T follicular helper cells")

for (ct in unique(markers$cell.type)) {
    plot_genes_by_group(
        monocle.obj,
        markers %>% filter(cell.type == ct) %>% pull(gene_symbol),
        group_cells_by = "ct_cluster",
        ordering_type = "none",
        max.size = 3
    ) +
        theme(axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
        )) +
        ggtitle(paste(ct, "markers"))
    if (length(markers %>% filter(cell.type == ct) %>% pull(gene_symbol)) >
        150) {
        ggsave(file.path(plotsDir, paste0(str_replace_all(ct," ","_"),"_markers.pdf")), height = 19)
    } else if (length(markers %>% filter(cell.type == ct) %>% pull(gene_symbol)) >
               70) {
        ggsave(file.path(plotsDir, paste0(str_replace_all(ct," ","_"),"_markers.pdf")), height = 12)
    } else if (length(markers %>% filter(cell.type == ct) %>% pull(gene_symbol)) >
               40) {
        ggsave(file.path(plotsDir, paste0(str_replace_all(ct," ","_"),"_markers.pdf")), height = 10)
    } else
        ggsave(file.path(plotsDir, paste0(str_replace_all(ct," ","_"),"_markers.pdf")))
}

marker_list <- list()
marker_list$pan_cd4 <-
  str_split(
    "Cd3, Cd4, Cd5, Cd7, Il2ra, Cd27, Cd28, Cd44, Sell, Cd69, Il7r, Tnfrsf4, Tnfrsf9, Ctla4, Cd40lg, Btla, Icos, Prdm1",
    pattern = ", "
  ) %>% unlist

marker_list$th1_cd4 <-
  str_split(
    "Ifng, Il2, Tnf, Lta, Klrd1, Ifngr1, Cxcr3, Cxcr6, Ccr1, Ccr5, Il12rb1, Tnfsf11, Havcr2, Stat1, Stat4, Tbx21",
    pattern = ", "
  ) %>% unlist

marker_list$th2_cd4 <- str_split(
  "Areg, Il4, Il5, Il9, Il10, Il13, Il21, Cxcr4, Ccr4, Ccr8, Ptgdr2, Hvacr1, Il17rb, Il1rl1, Batf, Gata3, Irf4, Stat6",
  pattern = ", "
) %>% unlist

marker_list$th17_cd4 <- str_split(
  "Csf2, Il17a, Il17af, Il17f, Il21, Il22, Ccr4, Ccr6, Il1r1, Il12b, Il21r, Il23r, Ahr, Batf, Maf, Irf4, Rora, Rorc, Stat3",
  pattern = ", "
) %>% unlist

marker_list$tfh_cd4 <- str_split(
  "Il21, Cxcr5, Icos, Prdm1, Batf, Bcl6, Maf, Irf4, Stat3",
  pattern = ", "
) %>% unlist

marker_list$treg_cd4 <- str_split(
  "Il10, Tgfb, Il2ra, Tnfrsf9, Entpd1, Nt5e, Cd101, Il1r1, Il1r2, Ctla4, Nrp1, Tnfrsf18, Cd109, Tigit, Rel, Ikzf4, Foxp3, Ikzf2, Stat5",
  pattern = ", "
) %>% unlist

marker_list$pan_cd8 <- str_split(
  "Cd3, Cd5, Cd8a, Cd27, Cd28",
  pattern = ", "
) %>% unlist

marker_list$cmem_cd8 <- str_split(
  "Cd44, Sell, Il7r, Ccr7, Bcl6, Eomes, Tbx21",
  pattern = ", "
) %>% unlist

marker_list$emem_cd8 <- str_split(
  "Gzmb, Ifng, Il2, Prf1, Tnf, Cd44, Klrg1, Eomes, Tbx21, Prdm1",
  pattern = ", "
) %>% unlist

marker_list$eff_cd8 <- str_split(
  "Ccl3, Ccl4, Ccl5, Gzmq, Gzmb, Ifng, Il2, Prf1, Tnf, Il2ra, Tnfrsf8, Cd44, Cd69, Il2rb, Tnfrsf4, Tnfrsf9, Lag3, Btla, Icos, Pdcd1, Havcr2, Klrg1, Prdm1, Id2, Tbx21",
  pattern = ", "
) %>% unlist

marker_list$naive_cd8 <- str_split(
  "Sell, Il7r, Ccr7, Cxcr3",
  pattern = ", "
) %>% unlist

marker_list$gd_t <- str_split(
  "Gzmb, Ifng, Il4, Il5, Il17a, Prf1, Tnf, Il2ra, Cd27, Cd28, Cd44, Ptprc, Sell, Cd69, Tfrc, Itgae, Il2rb, Il7r, Slamf1, Ccr6, Tcrg, Il1r1, Klrb1c, Eomes, Rorc, Tbx21",
  pattern = ", "
) %>% unlist

marker_list$nk_t <- str_split(
  "Ifng, Il4, Il10, Il13, Il17, Il21, Il22, Tnf, Cd4, Cd8a, Cd44, Itga2, Klrd1, Cd160, Cxcr3, Cxcr4, Cxcr6, Klrk1, Klra, Klrb1c, Trav14, Zbtb16, Egr2, Tcf12, Zfp683",
  pattern = ", "
) %>% unlist

marker_list$pan_b <- str_split(
  "Ptprc, Cd19, Ms4a1, Cd22, Fcer1g, Il2ra, Cd36, Cd40, Cd69, Cd80, Cd81, Cd86, Cd180, Tnfsf4, Ighd",
  pattern = ", "
) %>% unlist

marker_list$follicular_b <- str_split(
  "Ms4a1, Cd22, Fcer1g, Il2ra, Cd36, Cd40, Cd69, Cd80, Cd81, Cd86, Cd180, Tnfsf4, Ighd, Bcl6, Ebf1, Foxo1, Ikzf1, Pax5, Cd38, Cxcr5, H2-Ab1, Sell",
  pattern = ", "
) %>% unlist

marker_list$marginal_b <- str_split(
  "Cd1d1, Cd9, Cr2, Cd36, Cd81, Cd180, S1pr1, Fcgr4, Ighm, Ebf1, Notch2, Cd22, Tcf3, Pou2f2, Pax5",
  pattern = ", "
) %>% unlist

marker_list$plasma_b <- str_split(
  "Cd9, Fcgr3, Fcgr2b, Cd81, Sdc1, Ighg, Prdm1, Irf4, Xbp1, Tnfrsf17, Cd27, Cd38, Cd93, Slc3a2, Cxcr4, Ly6k",
  pattern = ", "
) %>% unlist

marker_list$mem_b <- str_split(
  "Cd38, Ptprc, Tnfrsf13b, Tnfrsf17, Ighg, Rbpj, Spib, Nt5e, Cd80, Cr2, Cd27, Cd40, Fas, Sell, H2-Ab1, Pdcd1lg2, Pou2af1, Pax5",
  pattern = ", "
) %>% unlist

marker_list$reg_b <- str_split(
  "Cd1d1, Cd5, Cr2, Cd24a, Cd40, Ighd, Ighm, Havcr1, Il10, Tgfb, Ebf1, Tcf3, Pou2f2, Pax5",
  pattern = ", "
) %>% unlist

marker_list$pan_dc <- str_split(
  "Itgax, Cd86, Cd80, H2-Ab1",
  pattern = ", "
) %>% unlist

marker_list$conventional_DC <- str_split(
  "Ido1, Il1b, Il6, Cxcl15, Il12a, Il12b, Il15, Il23a, Cd40, Cd83, Ccr7, Ly75, Cd208, Cd209a, Pdcd1lg2, Clec9a, Clec4a4, Flt3, Tlr3, Tlr9",
  pattern = ", "
) %>% unlist

marker_list$plasmacytoid_DC <- str_split(
  "Ifna, Ifnb, Ptprc, Bst2, Siglech, Tlr7, Tlr9, Tcf4, Irf8, Sirpa",
  pattern = ", "
) %>% unlist

marker_list$pan_MP <- str_split(
  "Axl, Itgam, Fcgr3, Fcgr2b, Cd40, Fcgr1, Cd68, Lamp2, Csf1r, Tlr2, Tlr4, Adgre1, Lgals3, Tnfsf18, Gpnmb, H2-Ab1, Vsig4, Vsir",
  pattern = ", "
) %>% unlist

marker_list$M2_MP <- str_split(
  "Arg1, Il10, Tgfb, Chil3, Cd14, Csf1r, Cd163, Msr1, Mrc1, Cd209a, Fcer1, Ly6c1, Mertk, Irf4, Retnla, Stat6 ",
  pattern = ", "
) %>% unlist

marker_list$M1_MP <- str_split(
  "Ido, Ifng, Il1, Il6, Il12a, Il12b, Il23, Tnf, Cd14, Fcgr3, Fcgr2b, Fcgr1, Cd86, Cd80, Ly6c1, H2-Ab1, Cd68, Irf5, Nos2, Stat1",
  pattern = ", "
) %>% unlist

marker_list$tissue_MP <- str_split(
  "Axl, Ccr2, Cd14, Cd68, Csf1r, Mrc1, Pdcd1lg2, Nos2",
  pattern = ", "
) %>% unlist

marker_list$microglia_MP <- str_split(
  "Ptprc, Csf1r, Cxcr1, Tmem119, Sall1, Siglech",
  pattern = ", "
) %>% unlist

marker_list$monocytes <- str_split(
  "Csf1r, Ccr2, Cxcr1, Ly6c1",
  pattern = ", "
) %>% unlist

marker_list$pan_granulocytes <- str_split(
  "Itgam, Fcgr3, Fcgr2b",
  pattern = ", "
) %>% unlist

marker_list$mast_cells <- str_split(
  "Tnf, Il4, Tgfb, Ngf, Cd9, Fut4, CD24a, Cr2, Spn, Fcgr1, Csf2ra, Kit, Il3ra, Il5ra, Il6ra, Fcer1, Il1rl1",
  pattern = ", "
) %>% unlist

marker_list$neutrophils <- str_split(
  "Elane, Il6, Il12, Tnf, Il1a, Il1b, Mme, Fut4, CD24a, Cr2, Spn, Ceacam1, Ceacam3, Igha, Cd93, Nectin2, Csf3r, Csf2ra, Bst1, Cd177, Cxcr1, Tlr2, Tlr4, Tlr6, Ly6g, Tlr1, Tlr9",
  pattern = ", "
) %>% unlist

marker_list$basophils <- str_split(
  "Il4, Il13, Ccl3, Cd9, Itgal, Anpep, Il2ra, Cd33, Cd38, Spn, Il3ra, Il5ra, Cd40lg, Ccr2, Il18r1, Il18rap, Tlr2, Tlr4, Tlr6, Ptgdr2, Fcer1, Tlr1, Tlr9, Cebpa, Gata2",
  pattern = ", "
) %>% unlist

marker_list$eosinophils <- str_split(
  "Epx, Mbp, Cd9, Fut4, CD24a, Cr2, Spn, Fcgr1, Csf2ra, Il3ra, Il5ra, Il6ra, Siglecf, Ccr3, Cd244a, Fcer1",
  pattern = ", "
) %>% unlist

marker_list$nk_cells <- str_split(
  "Gzma, Gzmb, Gzmk, Ifng, Il10, Prf1, Tnf, Itgam, Fcgr3, Il2ra, Itga2, Klrd1, Il2rb, Klrb1c, Cd226, Klrk1, Ncr1, Klra, Tigit, Nfil3, Eomes, Gata3, Id2, Runx1, Tbx21, Tox",
  pattern = ", "
) %>% unlist

marker_list$endothelial_cells <- str_split(
  "Pecam1, Cd34, Icam1, Itgb3, Sele, Eng, Vcam1, Cdh5, Mcam, Tek, Kdr, Pdpn, Flt4",
  pattern = ", "
) %>% unlist

marker_list$treg_cd4 <- str_split(
  "Il10, Tgfb, Il2ra, Tnfrsf9, Entpd1, Nt5e, Cd101, Il1r1, Il1r2, Ctla4, Nrp1, Tnfrsf18, Cd109, Tigit, Rel, Ikzf4, Foxp3, Ikzf2, Stat5",
  pattern = ", "
) %>% unlist()

for(ct in names(marker_list)){

  genes <- marker_list[[ct]]
  plot_genes_by_group(
    monocle.obj,
    genes,
    group_cells_by = "ct_cluster",
    ordering_type = "none",
    max.size = 3
  ) +
    theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) +
    ggtitle(paste0(ct,"markers", sep = "_"))
  
  ggsave(file.path(plotsDir, paste0(paste(ct,"markers", sep = "_"), ".pdf")))
  
}


tregs <-
  c(
    "Icos",
    "Tnfrsf18",
    "Tnfrsf9",
    "Klrg1",
    "Il1r1",
    "Satb1",
    "Nr4a",
    "Foxo",
    "AP-1",
    "Nfat",
    "Smad3",
    "Runx" ,
    "Foxp3",
    "Ets-1",
    "Creb",
    "Stat5",
    "c-Rel",
    "Ctla-4",
    "Lag3",
    "Havcr2",
    "Itgal",
    "Ccr7",
    "Ccr4",
    "Ccr6",
    "Ccr10",
    "Il1rl"
  ) %>% unique()

plot_genes_by_group(
    monocle.obj,
    tregs,
    group_cells_by = "ct_cluster",
    ordering_type = "none",
    max.size = 3
) +
    theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
    )) +
    ggtitle("T_reg_markers")
ggsave(file.path(plotsDir, "Treg_markers_Melanie.pdf"))

innate <- c("Cd80", "Cd86", "Il33")

plot_genes_by_group(
    monocle.obj,
    innate,
    group_cells_by = "ct_cluster",
    ordering_type = "none",
    max.size = 3
) +
    theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
    )) +
    ggtitle("Innate_markers")
    
ggsave(file.path(plotsDir, "Innate_markers_Melanie.pdf"))

Bcells <-
    c("Igha",
      "Ighd",
      "Ighe",
      "Ighg1",
      "Ighg2a",
      "Ighg2b",
      "Ighg2c",
      "Ighg3",
      "Ighm")

plot_genes_by_group(
    monocle.obj,
    Bcells,
    group_cells_by = "ct_cluster",
    ordering_type = "none",
    max.size = 3
) +
    theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
    )) + 
    ggtitle("Bcell_markers")


ggsave(file.path(plotsDir, "Bcell_markers_Melanie.pdf"))

Other <- paste0("Mmp",seq(27))
plot_genes_by_group(
    monocle.obj,
    Other,
    group_cells_by = "ct_cluster",
    ordering_type = "none",
    max.size = 3
) +
    theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
    )) +
    ggtitle("other_markers")

ggsave(file.path(plotsDir, "Other_markers_Melanie.pdf"))

Dendritic <-
  c(
    "H2-Ab1",
    "Itgax",
    "CD86",
    "CD80",
    "Tlr3",
    "Tlr7",
    "Tlr9",
    "Flt3",
    "Ly75",
    "Irf8",
    "Siglech",
    "Ptprc",
    "Bst2"
  )

plot_genes_by_group(
  monocle.obj,
  Dendritic,
  group_cells_by = "ct_cluster",
  ordering_type = "none",
  max.size = 3
) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  )) +
  ggtitle("DC_markers")

ggsave(file.path(plotsDir, "Dendritic_markers_Melanie.pdf"))

# plot_genes_by_group(
#   monocle.obj,
#   "Nt5e",
#   group_cells_by = "ct_cluster",
#   ordering_type = "none",
#   max.size = 3)


# Fraction assigned ------------------------------------------------------
for (tx in c("treatment",
             "treatment.agg")) {
  monocle.obj@colData %>%
    as_tibble %>%
    filter(run != "run1") %>% 
    group_by(celltype,organ, run, treat = get(tx) == "Undefined") %>%
    summarize(n = n()) %>%
    mutate(sum = sum(n), fraction_assigned = n / sum) %>%
    filter(treat == F) %>%
    dplyr::select(-c(treat, n, sum, organ)) %>%
    ggplot() +
    geom_col(aes(celltype, fraction_assigned, fill = run),
             position = "dodge") +
    theme_my() +
    theme(axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1
    )) +
    facet_wrap(~organ, nrow = 2) +
    ggtitle(paste0("Fraction of assigned cells | ", tx))
  ggsave(file.path(plotsDir,paste0("fraction_assigned_", tx,".pdf")))
}
 



abdata <- monocle.obj@colData %>%
  as_tibble %>%
  pivot_longer(c(hashid, cd4cd8, cd45x, gdcd4),
               names_to = "abgroup",
               values_to = "abtype") %>%
  group_by(ct_cluster, abgroup, abtype) %>%
  summarize(n = n()) %>%
  mutate(fraction = n / sum(n))

abdata %>% 
  ggplot +
  geom_tile(aes(abtype, ct_cluster, fill = fraction),
            col = "white", 
            width = 0.75,
            show.legend = F) +
  scale_fill_gradient(low = "ivory2", high = "red") +
  facet_grid(~abgroup,
             scales = "free",
             space = "free"
  ) +
  theme_my() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      size = 8,
      hjust = 1,
      vjust = 0.5
    )
  )

ggsave(file.path(plotsDir, "heatmap_clusters_abtype_unfiltered.pdf"), height = 11)


bind_rows(hashs, cd45s, cd4s) %>% 
  ggplot +
  geom_point(aes(id, ct_cluster, fill = fraction, size = log10(cells)),
             shape = 21, col = "black") +
  scale_fill_gradient(low = "ivory2", high = "red") +
  facet_grid(~abtype,
             scales = "free",
             space = "free"
  ) +
  theme_my() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      size = 8,
      hjust = 1,
      vjust = 0.5
    )
  )


ggsave(file.path(plotsDir, "heatmap_clusters_abtype_cellcounts.pdf"), height = 14)

x <- monocle.obj@colData %>% 
  as_tibble() %>% 
  filter(treatment.agg != "Undefined") %>%  
  group_by(experiment,
           celltype,
           treatment.agg,
           HDAC = grepl("HDAC", treatment.agg)) %>%
  summarize(n = n())

monocle.obj@colData %>% 
  as_tibble() %>% 
  filter(treatment.agg != "Undefined") %>%  
  group_by(experiment,
           celltype,
           organ,
           treatment.agg) %>% 
  summarize(n = n()) %>% 
  spread(treatment.agg, n) %>% 
  mutate(WTvsCTRL = HDAC_WT/NoT,
         cKOvsCTRL = HDAC_cKO/NoT) %>% 
  pivot_longer(cols = c(WTvsCTRL, cKOvsCTRL), names_to = "comparison", values_to = "FC") %>% 
  filter(experiment == "HDAC1") %>% 
  ggplot() +
  geom_tile(aes(comparison, celltype, fill = log2(FC), width = 0.75)) +
  scale_fill_gradient2(low = "#00688b", mid = "ivory2", high = "#A11504", na.value = 'darkgreen') +
  facet_wrap(~organ) + 
  theme_my() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      size = 8,
      hjust = 1,
      vjust = 0.5
    )
  )

ggsave(file.path(plotsDir, "heatmap_ct_treatgroup.pdf"))

monocle.obj@colData %>% 
  as_tibble() %>% 
  filter(treatment.agg != "Undefined") %>%  
  group_by(experiment,
           celltype,
           organ,
           treatment.agg) %>% 
  summarize(n = n()) %>% 
  mutate(total = sum(n), fraction = n/total*100) %>% 
  filter(experiment == "HDAC1") %>% 
  ggplot() +
  geom_col(aes(celltype, fraction, fill = treatment.agg), position = "stack")+
  facet_wrap(~organ) + 
  ylab("fraction in %") +
  ggtitle("Fraction of assigned treatment per celltype in HDAC1 samples")+
  scale_y_continuous(expand = expansion(mult = c(0,0.05),
                                        add = c(1,0)))+
  theme_my() +
  theme(
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      size = 8,
      hjust = 1,
      vjust = 0.5
    ),
    panel.spacing.y = unit(0,"mm")
  )

ggsave(file.path(plotsDir,"assigned_treatment_ct.pdf"))


monocle.obj@colData %>% 
  as_tibble() %>% 
  group_by(sample, ct_cluster, treatment.agg) %>% 
  summarise(n = n()) %>% 
  mutate(total = sum(n)) %>% 
  filter(treatment.agg == "Undefined") %>% 
  mutate(percent.unassigned = n/total*100) %>% 
  ggplot() +
  geom_tile(aes(sample,ct_cluster, fill = percent.unassigned)) + 
  scale_fill_gradient2(low = "#00688b", mid = "ivory2", high = "#A11504", na.value = 'darkgreen') +
  theme_my() +
  theme(
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      size = 8,
      hjust = 1,
      vjust = 0.5
    ),
    panel.background=element_rect(fill="grey45", colour="black")
  ) +
  ggtitle("Percent Undefined Treatment Grp per sample and cluster")


ggsave(file.path(plotsDir, "percent_undefined.pdf"), height = 12)

monocle.obj@colData %>% 
  as_tibble() %>% 
  ggplot() +
  geom_boxplot(aes(ct_cluster, gd.counts+1), fill = "salmon") +
  scale_y_log10() +
  theme_my() +
  ggtitle("gamma_delta_tag_counts")

ggsave(file.path(plotsDir, "gamma_delta_counts.pdf"))

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

# sc vs bulk --------------------------------------------------------------

res_bulk <- readRDS("res_bulk.rds")
res_MA <- readRDS(file.path(dataDir,"DGE_res_MA.rds"))
res_MA_pb <- readRDS(file.path(dataDir, "DGE_res_MA_pb.rds"))

res_bulk$celltype <- case_when(res_bulk$coef %in% c("disease") ~ "Fibroblasts",
                               res_bulk$coef %in% c("disease_Keratinocytes") ~ "Keratinocytes")
res_bulk <- res_bulk %>% filter(!is.na(celltype))
res_bulk$P.Value <- 
  case_when(
    res_bulk$P.Value == 0 ~ 5e-324,
    TRUE ~ res_bulk$P.Value
  )

res_bulk$adj.P.Val <- 
  case_when(
    res_bulk$adj.P.Val == 0 ~ 5e-324,
    TRUE ~ res_bulk$adj.P.Val
  )
res_MA_pb <- res_MA_pb %>% filter(coef == "HDAC_WT")
res_bulk <- res_bulk %>% filter(rn %in% res_MA_pb$rn)
res_bulk <- res_bulk %>% filter(adj.P.Val < 0.05)

# res_MA <- res_MA %>% filter(rn %in% res_bulk$rn)
res_MA_pb <- res_MA_pb %>% filter(rn %in% res_bulk$rn)

res_MA <- res_MA_pb

cor_list <- list()
gene_list <- list()
for(ct in unique(res_MA$celltype)) {
  subset_MA <- res_MA %>% filter(celltype == ct,
                                 coef == "HDAC_WT")
  subset_bulk <- res_bulk %>% filter(celltype == ct,
                                     grepl("disease", coef),
                                     !grepl("\\:", coef))
  for (clusterx in unique(subset_MA$cluster)) {
    subset_sc <- 
      subset_MA %>% 
      filter(cluster == as.numeric(clusterx),
             rn %in% subset_bulk$rn
             ) %>% 
      mutate(signed.p = -log10(P.Value) * sign(logFC)) %>% 
      select(rn, logFC, signed.p)
    
    subset_b <-
      subset_bulk %>% 
      filter(rn %in% subset_sc$rn) %>% 
      mutate(signed.p = -log10(P.Value) * sign(logFC)) %>% 
      select(rn, logFC, signed.p)
    
    df <- inner_join(subset_sc, subset_b, by = "rn")
    
    colnames(df) <-
      case_when(
        grepl("\\.x", colnames(df)) ~ str_replace(colnames(df), "\\.x", "_sc"),
        grepl("\\.y", colnames(df)) ~ str_replace(colnames(df), "\\.y", "_b"),
        TRUE ~ colnames(df)
      )
  
    cor_list[[ct]][[clusterx]][["logFC"]] <- cor(df$logFC_sc, df$logFC_b, method = "spearman")
    cor_list[[ct]][[clusterx]][["-log10(P.Value) * sign(logFC)"]] <- cor(df$signed.p_sc, df$signed.p_b, method = "spearman")
    
    gene_list[[ct]][[clusterx]] <- df

  }
}

cors <- bind_rows(lapply(cor_list, function(l1) {
  bind_rows(l1, .id = "cluster")
}),
.id = "celltype")

genes <- bind_rows(lapply(gene_list, function(l1) {
  bind_rows(l1, .id = "cluster")
}),
.id = "celltype")

cors <- cors %>% pivot_longer(c(logFC, "-log10(P.Value) * sign(logFC)"), values_to = "Correlation") %>% filter(!is.na(Correlation))

cors %>% 
  ggplot() +
  geom_col(aes(cluster, Correlation),
           col = "black",
           fill = "ivory2") +
  facet_grid(name~celltype, scales = "free") +
  theme_my()+
  ggtitle("Significant bulk Genes")

ggsave("sig_b_cor_pb_bulk.pdf")

genes %>% 
  ggplot(aes(logFC_sc, logFC_b)) +
  geom_point() + 
  #geom_text_repel(aes(label = rn)) +
  facet_wrap(~ cluster, scales = "free") +
  geom_smooth(method = "lm", se = T) +
  theme_my()

ggsave(file.path(plotsDir, "logFC_cor_pb.pdf"))

genes %>% 
  ggplot(aes(signed.p_sc, signed.p_b)) +
  geom_point() + 
  #geom_text_repel(aes(label = rn)) +
  facet_wrap(~ cluster, scales = "free") +
  geom_smooth(method = "lm", se = T) +
  theme_my()

ggsave(file.path(plotsDir, "pval_cor_pb.pdf"))

enriched_terms <- enrichr(
  genes %>%
    filter(cluster %in% c(24,48)) %>% pull(rn),
  c(
    "KEGG_2019_Mouse",
    "WikiPathways_2019_Mouse",
    "GO_Biological_Process_2021",
    "MSigDB_Hallmark_2020",
    "Disease_Signatures_from_GEO_up_2014",
    "Jensen_DISEASES",
    "Rare_Diseases_GeneRIF_Gene_Lists"
  )
)

enriched_terms <- rbindlist(enriched_terms, idcol = "database")

# enriched_terms %>% group_by(database) %>% slice_max(n = 3, Combined.Score) %>% View()

plotEnrich(enriched_terms[which(enriched_terms$Genes %>% str_count(";")>2),],showTerms = 10, orderBy = "Combined.Score")
ggsave(file.path(plotsDir, "enrichment_bulk_vs_sc.pdf"))



