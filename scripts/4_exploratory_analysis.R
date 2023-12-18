library(dplyr)
library(tidyr)
library(stringr)
library(monocle3)
library(data.table)
library(ggplot2)
library(readr)

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
gfsDir <- '/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin'
plotsDir <- file.path(gfsDir, 'plots')
tablesDir <- file.path(gfsDir, 'tables')
oldDir <- "/vscratch/scRNAseq/data/old"
shinyDir <- 'dge-app'
dataDir <-"data"
resDir <- "results"


# Load Data ---------------------------------------------------------------
LUT <- read_rds(file.path(dataDir, 'doublet_LUT.rds'))
setDT(LUT)

monocle.obj <- read_rds(file.path(dataDir, "/3_annotated_monocle.cds"))
mon.old <- read_rds(file.path(oldDir, "scRNAseq_3_monocle_more_hash_cutoff.cds"))
# new.clusters <- read.csv(file.path(tablesDir, "NewClusterML.csv"))
# laia.clusters <- read.csv(file.path(tablesDir, "Ludwig4.csv"))

mon.old <- mon.old[,!grepl("fLN_40B3", colnames(mon.old))]
mon.old@colData <- mon.old@colData %>% droplevels()
mon.old@colData$doublet <- LUT[match(colnames(mon.old), rn)]$doublet_classification
# mon.old <- mon.old[,colnames(mon.old) %in% colnames(monocle.obj)]

# stopifnot(all(rownames(colData(monocle.obj)) == rownames(colData(mon.old))))
# monocle.obj@colData$Cluster_old_unfiltered <- mon.old@colData$Cluster

# monocle.obj@colData$same_ct <- monocle.obj@colData$celltype == monocle.obj@colData$celltype_old
# Exploratory Analysis ----------------------------------------------------

# dt <- as.data.table(monocle.obj@colData, keep.rownames = TRUE)

# dt <- dt[,.(celltype, celltype_old, Cluster, Cluster_old_unfiltered, cd45x, cd4cd8, rn)]
# 
# dt[celltype != celltype_old, .N, keyby = .(celltype, celltype_old)] %>% 
#   ggplot(aes(celltype, celltype_old)) + 
#   geom_point(aes(size = log10(N))) + 
#   theme_my() +
#   ggtitle("old vs new celltype") +
#   coord_flip()
# ggsave(file.path(dataDir, "ct_old_vs_new.pdf"))
# 
# pdf(file.path(dataDir, "ct_numbers_diff.pdf"))
# c("diff" = sum(dt[,celltype != celltype_old]), "same" = sum(dt[,celltype == celltype_old])) %>% 
#   barplot()
# dev.off()
# 

# Save Object -------------------------------------------------------------
# saveRDS(monocle.obj, file = file.path(dataDir, "scRNAseq_2a_monocle.cds"))

## Assignment ------------------------------------------------------


abdata <- monocle.obj@colData %>%
    as_tibble %>%
    pivot_longer(c(hashid, cd4cd8, cd45x, gdcd4),
                 names_to = "abgroup",
                 values_to = "abtype") %>%
    group_by(ct_cluster, abgroup, abtype) %>%
    summarize(n = n()) %>%
    mutate(fraction = n / sum(n))

abdata %>% 
  filter(abtype != "Undefined") %>% 
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

ggsave(file.path(plotsDir, "heatmap_clusters_abtype.pdf"), height = 11)

# UMAPS - CITE seq --------------------------------------------------------
coords <- reducedDims(monocle.obj)
coords <- data.frame(coords$UMAP)
colnames(coords) <- c("UMAP1", "UMAP2")
coords$experiment <- colData(monocle.obj)$experiment
coords$cd45x <- colData(monocle.obj)$cd45x
coords$cd4cd8 <- colData(monocle.obj)$cd4cd8
coords$hashid <- colData(monocle.obj)$hashid
coords$gdcd4 <- colData(monocle.obj)$gdcd4
coords$treatment <- colData(monocle.obj)$treatment
coords$hashid.agg <- colData(monocle.obj)$hashid.agg
coords$treatment.agg <- colData(monocle.obj)$treatment.agg


ggplot(coords) + 
    geom_hex(data = coords %>% dplyr::select(-cd45x), 
             aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
    geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
    scale_fill_gradient(low = "deepskyblue2", high = "red", limits = c(0,40), oob = scales::squish)+
    facet_wrap( ~cd45x, nrow = 2) +
    theme_my() 

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

plot_cells(monocle.obj, color_cells_by = "cluster")
ggsave(file.path(plotsDir, "UMAP_cluster.jpg"))


plot_cells(monocle.obj, group_cells_by="cluster", 
           color_cells_by="celltype",
           label_groups_by_cluster=F,
           group_label_size = 4) 

ggsave(file.path(plotsDir, "UMAP_annotated.jpg"))

plot_cells(monocle.obj, group_cells_by="cluster", 
           color_cells_by="celltype_raw",
           label_groups_by_cluster=F,
           label_cell_groups = F,
           group_label_size = 4) 

ggsave(file.path(plotsDir, "UMAP_annotated_raw.jpg"))


plot_cells(mon.old, group_cells_by = 'cluster',
           color_cells_by="celltype",
           label_groups_by_cluster=F,
           group_label_size = 4) 

ggsave(file.path(plotsDir, "old_UMAP_annotated.jpg"))

old_co <- reducedDims(mon.old)
old_co <- data.frame(old_co$UMAP)
colnames(old_co) <- c("UMAP1", "UMAP2")
old_co$experiment <- colData(mon.old)$experiment
old_co$cd45x <- colData(mon.old)$cd45x
old_co$cd4cd8 <- colData(mon.old)$cd4cd8
old_co$hashid <- colData(mon.old)$hashid
old_co$gdcd4 <- colData(mon.old)$gdcd4
old_co$treatment <- colData(mon.old)$treatment
old_co$hashid.agg <- colData(mon.old)$hashid.agg
old_co$treatment.agg <- colData(mon.old)$treatment.agg
old_co$doublet <- colData(mon.old)$doublet

ggplot(old_co) + 
  geom_hex(data = old_co %>% dplyr::select(-doublet), 
           aes(x = UMAP1, y = UMAP2), fill = "ivory3", bins  = 200) + 
  geom_hex(aes(x = UMAP1, y = UMAP2), bins = 200) + 
  scale_fill_gradient(low = "deepskyblue2", high = "red", limits = c(0,40), oob = scales::squish) +
  facet_wrap( ~doublet, nrow = 2) +
  theme_my()

ggsave(file.path(plotsDir, "old_UMAP_doublet.jpg"))
# plot_cells(monocle.obj, group_cells_by="cluster", 
#            color_cells_by="same_ct",
#            label_groups_by_cluster=F,
#            group_label_size = 4) 
# 
# ggsave(file.path(plotsDir, "UMAP_same_ct.jpg"))


# Markers -----------------------------------------------------------------

## Krt5 expression among ct_clusters

monocle.obj[, monocle.obj@colData %>% subset(organ == "LN" &
                                               experiment == "HDAC1") %>% droplevels %>% rownames()] %>%
  plot_genes_by_group(
    c("Krt5", "Cd4", "Foxp3"),
    group_cells_by = "ct_cluster",
    ordering_type = "none",
    max.size = 3
  ) + theme_my() +
  ggtitle("HDAC1 LN - Krt5 expression across cts")

ggsave(file.path(plotsDir, "krt_HDAC1_LN.pdf"))

## Top Marker Genes
marker_test_res <- top_markers(monocle.obj, group_cells_by="Cluster",
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(4, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(monocle.obj,
                    top_specific_marker_ids,
                    group_cells_by="ct_cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  )
  )

ggsave(file.path(plotsDir, "top_markers_diag_order.pdf"), height = 15, width = 10)


##Custom Markers

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



# Lists for cloupe --------------------------------------------------------

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

write.csv(UMAP_projection, file = file.path(tablesDir, "UMAP-Projection-Nov-2023.csv"), quote = F, row.names = F)

graph_based <- monocle.obj@colData %>% as_tibble(rownames = "Barcode") %>% dplyr::select(Barcode, ct_cluster)
graph_based$Nov2023 <- paste("Cluster", graph_based$ct_cluster)
graph_based <- graph_based %>% dplyr::select(-ct_cluster)

graph_based <- rename_barcodes(graph_based)

write.csv(graph_based, file = file.path(tablesDir, "Graph-Based-Nov-2023.csv"), quote = F, row.names = F)

treatments <- monocle.obj@colData %>% as_tibble(rownames = "Barcode") %>% dplyr::select(Barcode, treatment.agg)

treatments <- rename_barcodes(treatments)

write.csv(treatments, file = file.path(tablesDir, "Treatments-Nov-2023.csv"), row.names = F)

experiments <- monocle.obj@colData %>% as_tibble(rownames = "Barcode") %>% dplyr::select(Barcode, experiment)

experiments <- rename_barcodes(experiments)

write.csv(experiments, file = file.path(tablesDir, "experiments-Nov-2023.csv"), row.names = F)

exp_treat <- inner_join(experiments,treatments) %>%
    mutate(exp_treat = paste(experiment,
                             str_extract(treatment.agg, "[:alpha:]+$"),
                             sep = "_")) %>%
    dplyr::select(Barcode, exp_treat)

write.csv(exp_treat, file = file.path(tablesDir, "exp-treat-Nov-2023.csv"), row.names = F)


