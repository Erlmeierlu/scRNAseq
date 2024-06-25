library(renv)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(monocle3)
library(ggrepel)
library(data.table)
library(limma)
library(Matrix.utils)
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
    strip.text.x = element_text(size = 10, margin = margin(b = 2, t = 2)),
    strip.background = element_rect(
      fill = "#9FD7D2",
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

monocle.obj <- readRDS(file.path(dataDir, "/scRNAseq_3_monocle_more_hash_cutoff.cds"))

t.cds <- readRDS(file.path(dataDir, "/scRNAseq_3_t_subset.cds"))

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

for(samplex in unique(monocle.obj@colData$sample)) {
  for (ct in unique(monocle.obj@colData$celltype)) {
    for (treat in unique(monocle.obj@colData$treatment55[monocle.obj@colData$treatment55 != "Undefined"])) {
      
      subsetx <- monocle.obj[, monocle.obj@colData %>% subset(celltype == ct &
                                                                treatment55 == treat &
                                                                sample == samplex) %>% rownames()]
      
      if(ncol(subsetx) == 0) next
      
      x <-  counts(subsetx) %>% rowSums() %>% as("sparseMatrix")
      
      
      countsx <- cbind(countsx, x)
      
      
      colnames(countsx)[ncol(countsx)] <-
        paste(samplex, str_replace_all(ct, " |\\/", "_"), treat, sep = "_")
      
    }
  }
}

countsx <- countsx[,-1]

designx <- monocle.obj@colData %>% as_tibble() %>% filter(treatment08 != "Undefined") %>%  group_by(celltype, treatment08, sample, sex, organ, experiment) %>% 
  summarise(
    n = n(),
    .groups = "drop") %>% mutate(rn = paste(sample, str_replace_all(celltype, " |\\/", "_"), treatment08, sep = "_"))
designx <- designx %>% tibble::column_to_rownames("rn")
designx <- designx[colnames(countsx),]
designx <- designx %>% rename(treatment08 = "treatment")
designx$treatment <- forcats::fct_relevel(designx$treatment, c("NoT", "HDAC_WT", "HDAC_cKO")) %>% droplevels()

# saveRDS(countsx, file = file.path(dataDir, "counts_pseudobulk55.rds"))
saveRDS(designx, file = file.path(dataDir, "design_pseudobulk_08.rds"))

### Limma -------------------------------------------------------------------

countsx <- readRDS(file = file.path(dataDir, "counts_pseudobulk_08.rds"))
designx <- readRDS(file = file.path(dataDir, "design_pseudobulk.rds"))

comp <- combn(unique(levels(
  monocle.obj[, colData(monocle.obj) %>%
                subset(treatment != "Undefined") %>% row.names]@colData$treatment %>%
    droplevels
)),
2,
FUN = paste,
collapse = ";")

results.DGE.pseudo <- list()
for (ct in unique(designx$celltype)) {
  for (organx in unique(designx$organ)) {
    for (compx in comp) {
      for (expx in unique(designx$experiment)) {
        
        compxx <- unlist(str_split(compx, ";"))
        
        subsetx <- designx %>% 
          subset(
            celltype == ct &
            organ == organx &
            treatment %in% compxx &
            experiment == expx
        )
        
        if(nrow(subsetx) < 2) next
        
        countsx.pseudo <- countsx[,rownames(subsetx)]
        
        # corx <- cor(as.matrix(countsx.pseudo))
        # 
        # for(i in compxx){
        #   
        # if(nrow(as.matrix(corx[grepl(i, rownames(corx)), grepl(i, colnames(corx))])) < 2) next
        #   
        # cor_heat <- pheatmap::pheatmap(corx[grepl(i, rownames(corx)),grepl(i, colnames(corx))],
        #                                main = paste0(
        #                                  i,
        #                                  "_",
        #                                  str_replace_all(ct, " |\\/", "_"),
        #                                  "_",
        #                                  organx,
        #                                  "_",
        #                                  expx
        #                                ))
        # dev.off()
        # 
        # pdf(file.path(
        #   plotsDir,
        #   paste0(
        #     "cor_heat_",
        #     i,
        #     "_",
        #     str_replace_all(ct, " |\\/", "_"),
        #     "_",
        #     organx,
        #     "_",
        #     expx,
        #     ".pdf"
        #   )))
        # 
        # print(cor_heat)
        # 
        # dev.off()
        # }
        # avg_cor <- colSums(corx)/ncol(corx)
        #countsx.pseudo <- countsx.pseudo[,avg_cor > 0.8]
        
        if (length(unique(subsetx$treatment)) >= 2 &
            length(unique(subsetx$sex)) >= 2) {
          frame <-
            model.frame(~ treatment + sex,
                        subsetx %>% droplevels())
          designmat <- model.matrix(~ treatment + sex, frame)
          
          if (ncol(designmat) >= nrow(designmat)){
            
            
            next
          }
          
        } else if (length(unique(subsetx$treatment)) >= 2 &
                   length(unique(subsetx$sex)) < 2) {
          frame <-
            model.frame(~ treatment, subsetx %>% droplevels())
          designmat <- model.matrix(~ treatment, frame)
          
          if (ncol(designmat) >= nrow(designmat)){
            
            next
          }
          
        } else {
          
          next
        }
        colnames(designmat) <-
          colnames(designmat) %>%
          str_remove_all("\\(|\\)|sample|treatment|\\/|[0-9]\\/[0-9]|sex") %>%
          str_replace_all("\\-", "_")
        
        
        # gfilter <-
        #   filterByExpr(
        #     countsx,
        #     design = designmat,
        #     min.prop = 0.1,
        #     min.count = 1
        #   )
        
        gfilter <- rowSums(countsx.pseudo[,rownames(subsetx)]) > 0
        
        voomx <-
          voom(
            countsx.pseudo[gfilter,],
            design = designmat,
            plot = F,
            save.plot = T
          )
        
        
        
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
        #     "pseudo_voom_",
        #     str_replace_all(paste0(unlist(
        #       str_split(compx, ";")
        #     ), collapse = "_"), "-|\\/", "_"),
        #     "_",
        #     str_replace_all(ct, " |\\/", "_"),
        #     "_",
        #     organx,
        #     "_",
        #     expx,
        #     ".pdf"
        #   )
        # ))
        
        fitx <- lmFit(voomx, design = designmat)
        fitx <- eBayes(fitx)
        
        for (coefx in colnames(designmat)) {
          topx <- topTable(fitx, number = Inf, coef = coefx)
          topx$rn <- rownames(topx)
          results.DGE.pseudo[[ct]][[organx]][[compx]][[expx]][[coefx]] <-
            data.table(topx)
        }
        
      }
      
    }
    
  }
  
}

### Res ---------------------------------------------------------------------

res.pb <-
  rbindlist(lapply(results.DGE.pseudo, function(l1) {
    rbindlist(lapply(l1, function(l2) {
      rbindlist(lapply(l2, function(l3) {
        rbindlist(lapply(l3, function(l4) {
          rbindlist(l4, idcol = "coef")}), 
          idcol = "experiment")
      }), idcol = "treatment")
    }), idcol = "organ")
  }), idcol = "celltype")

saveRDS(res.pb, file.path(dataDir, "DGE_raw_pb.rds"))

res.pb$treatment <-
  case_when(
    res.pb$treatment %in% "NoT;HDAC_WT" ~ "WT_vs_ctrl",
    res.pb$treatment %in% "NoT;HDAC_cKO" ~ "KO_vs_ctrl",
    res.pb$treatment %in% "HDAC_WT;HDAC_cKO" ~ "KO_vs_WT",
    TRUE ~ res.pb$treatment
  )

res.pb$P.Value <- 
  case_when(
    res.pb$P.Value == 0 ~ 5e-324,
    TRUE ~ res.pb$P.Value
  )

res.pb$adj.P.Val <- 
  case_when(
    res.pb$adj.P.Val == 0 ~ 5e-324,
    TRUE ~ res.pb$adj.P.Val
  )

res.pb$celltype <- res.pb$celltype %>% 
  str_replace_all(" |\\/|\\-", "_") %>% 
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



saveRDS(res.pb, file.path(dataDir, "DGE_res_pb_55.rds"))

shiny_list <- res.pb %>% filter(!grepl("M|Inter", coef))

saveRDS(shiny_list, file.path(dataDir, "list_for_shiny_55.rds"))



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
              "MA_pseudo_voom_",
              clusterx,
              str_replace_all(ct, " |\\/", "_"),
              "_",
              ".pdf"
            )
          ))
          
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

res_bulk <- readRDS(file.path(bulkDir,"res_bulk.rds"))
res_MA <- readRDS(file.path(dataDir,"DGE_res_MA.rds"))
res_MA_pb <- readRDS(file.path(dataDir, "DGE_res_MA_pb.rds"))

res_bulk$celltype <- case_when(res_bulk$coef %in% c("Intercept", "disease") ~ "Fibroblasts",
                               TRUE ~ "Keratinocytes")

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

res_bulk <- res_bulk %>% filter(rn %in% res_MA$rn)

res_bulk <- res_bulk %>% filter(adj.P.Val < 0.05)

res_MA <- res_MA %>% filter(rn %in% res_bulk$rn)
res_MA_pb <- res_MA_pb %>% filter(rn %in% res_bulk$rn)

# res_MA <- res_MA_pb

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

ggsave(file.path(plotsDir, "sig_b_cor_pb_bulk.pdf"))

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
    "TRANSFAC_and_JASPAR_PWMs",
    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
    "TRRUST_Transcription_Factors_2019",
    "MSigDB_Hallmark_2020",
    "Disease_Signatures_from_GEO_up_2014",
    "Jensen_DISEASES",
    "Mouse_Gene_Atlas",
    "OMIM_Disease",
    "Rare_Diseases_GeneRIF_Gene_Lists"
    
  )
)

enriched_terms <- rbindlist(enriched_terms, idcol = "database")

# enriched_terms %>% group_by(database) %>% slice_max(n = 3, Combined.Score) %>% View()

plotEnrich(enriched_terms[which(enriched_terms$Genes %>% str_count(";")>2),],showTerms = 10)
ggsave(file.path(plotsDir, "enrichment_bulk_vs_sc.pdf"))



# test --------------------------------------------------------------------

test <- readRDS(file.path(dataDir, "list_for_shiny.rds"))
test[, treatment := factor(treatment, levels = c("WT_vs_ctrl", "KO_vs_ctrl", "KO_vs_WT"))]

test[rn == "Foxp3", .(round(AveExpr), celltype, treatment, experiment, organ)]


# for(ct in unique(test$celltype)) {
  for (organx in unique(test$organ)) {
    for (expx in unique(test$experiment)) {
      subsetx <- test %>% filter(
        # celltype == ct,
                                   organ == organx,
                                   experiment == expx,
                                   !grepl("Intercept|M", coef)) 
      if (nrow(subsetx) == 0) next
      
      
      subsetx %>% 
        ggplot() +
        geom_histogram(aes(P.Value, fill = factor(round(AveExpr))),
                       bins = 30
        ) +
        facet_wrap(~ treatment) +
        theme_my() +
        theme(panel.grid.major = element_blank()) +
        ggtitle(paste("pseudo_histo", 
                      # ct, 
                      organx, expx, sep = "_"))
      
      ggsave(file.path(plotsDir,
                       paste0(
                         paste("pseudo_histo", 
                               # ct, 
                               organx, expx, sep = "_"),
                         ".pdf"
                       )))
      
    }
  }
 # }
