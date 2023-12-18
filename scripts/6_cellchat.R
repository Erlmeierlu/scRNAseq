library(monocle3)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(CellChat)
library(Seurat)
library(readr)
library(data.table)
library(stringr)
library(foreach)

#personal theme
theme_my <- function(...) {
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
    ) + theme(...)
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
monocle.obj <- read_rds(file.path(dataDir, '3_annotated_monocle.cds'))

# Create CellChat Object --------------------------------------------------

doParallel::registerDoParallel(cores = 7)
foreach(organ = c('LN', 'Skin')) %dopar% {
    foreach(experiment = c('HDAC1', 'HDAC2')) %dopar% {
        foreach(treatment = c('NoT', 'HDAC_WT', 'HDAC_cKO')) %dopar% {
            
            subset <- monocle.obj[, monocle.obj@colData %>% subset(
                organ == organ & 
                    experiment == experiment & 
                    treatment.agg == treatment
            ) %>% droplevels %>% rownames()]
            
            subset@colData <- subset@colData %>% droplevels()
            
            counts <- counts(subset)
            counts <- normalizeData(counts)
            meta <- data.frame(subset@colData)
            
            subset <- createCellChat(
                counts,
                meta = meta,
                group.by = 'celltype'
            )
            
            #prep
            subset@DB <- CellChatDB.mouse
            subset <- subsetData(subset)
            
            subset <- subset %>% 
                identifyOverExpressedGenes() %>% 
                identifyOverExpressedInteractions()
            
            subset <- subset %>% 
                computeCommunProb() %>% 
                filterCommunication()
            
            subset <- computeCommunProbPathway(subset)
            subset <- aggregateNet(subset)
            
            subset <- netAnalysis_computeCentrality(subset)
            df <- subsetCommunication(subset)
            
            write_rds(subset, paste0('/vscratch/scRNAseq/data/', 
                                     'cellchat_',
                                     organ, '_', 
                                     experiment, '_', 
                                     treatment, '.rds')
            )
            
            fwrite(df, file.path(tablesDir, 
                                 paste0('cellchat_', 
                                        organ, '_',
                                        experiment, '_',
                                        treatment, '.csv'))
            )
        }
    }
}

    
files <- list.files(tablesDir, pattern = 'cellchat')

cc_dat <- sapply(files,
                 function(x)
                     fread(file.path(tablesDir, x)),
                 simplify = F)

cc_dat <- rbindlist(cc_dat, idcol = 'id')

cc_dat[, c('organ', 'experiment', 'treatment') := tstrsplit(
    str_replace_all(
        str_remove_all(id, '^[:alpha:]+_|.[:alpha:]+$'), 
        'HDAC_', 'HDAC-'),
    '_')]

cc_dat[, c('treatment', 'id') := .(factor(
    str_replace(treatment, '-', '_'), levels = c('NoT', 'HDAC_WT', 'HDAC_cKO')),
    NULL)]

cc_dat[, keyby = .(source, target), N := .N][pval < 0.05] %>%
    ggplot(aes(source, target, fill = N)) +
    geom_tile() +
    facet_grid(organ + experiment ~ treatment) +
    theme_my()


# ### HEATMAP: Difference in pathway counts of all treatments - Source Cells ###
# #Filter + calculate pathway interaction count difference between data frames 
# 
# df.net.noT <- df.net.noT %>%
#     filter(pval < 0.05) %>%
#     group_by(source, pathway_name) %>%
#     summarise(interaction_count_noT = n()) %>%
#     filter(interaction_count_noT > 5) 
# 
# df.net.ko <- df.net.ko %>%
#     filter(pval < 0.05) %>%
#     group_by(source, pathway_name) %>%
#     summarise(interaction_count_ko = n()) %>%
#     filter(interaction_count_ko > 5)
# 
# df.net.wt <- df.net.wt %>%
#     filter(pval < 0.05) %>%
#     group_by(source, pathway_name) %>%
#     summarise(interaction_count_wt = n()) %>%
#     filter(interaction_count_wt > 5)
# 
# merged_df <- inner_join(df.net.noT, df.net.ko, by = c("source", "pathway_name")) %>%
#     inner_join(., df.net.wt, by = c("source", "pathway_name")) %>%
#     mutate(
#         diff_noT_wt = interaction_count_wt - interaction_count_noT,
#         diff_noT_ko = interaction_count_ko - interaction_count_noT,
#         diff_wt_ko = interaction_count_ko - interaction_count_wt
#     )
# 
# min_value <- min(c(merged_df$diff_noT_wt, merged_df$diff_noT_ko, merged_df$diff_wt_ko))
# max_value <- max(c(merged_df$diff_noT_wt, merged_df$diff_noT_ko, merged_df$diff_wt_ko))
# 
# #Heatmap noT vs wt
# heatmap.noT.wt <- merged_df %>%
#     ggplot(aes(x = source, y = pathway_name, fill = diff_noT_wt)) +
#     geom_tile() +
#     ggtitle("Wildtype vs. No T cells") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
#     scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value))
# 
# #Heatmap noT vs ko
# heatmap.noT.ko <- merged_df %>%
#     ggplot(aes(x = source, y = pathway_name, fill = diff_noT_ko)) +
#     geom_tile() +
#     ggtitle("Knockout vs. No T cells") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
#     scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value))
# 
# #Heatmap wt vs ko
# heatmap.wt.ko <- merged_df %>%
#     ggplot(aes(x = source, y = pathway_name, fill = diff_wt_ko)) +
#     geom_tile() +
#     ggtitle("Knockout vs. Wildtype") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
#     scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value))
# 
# 
# grid.arrange(heatmap.noT.wt, heatmap.noT.ko, heatmap.wt.ko, ncol = 3)
# 
# ### HEATMAP: Difference in pathway counts of all treatments - Target Cells ###
# 
# df.net.noT <- df.net.noT %>%
#     filter(pval < 0.05) %>%
#     group_by(target, pathway_name) %>%
#     summarise(interaction_count_noT = n()) %>%
#     filter(interaction_count_noT > 5)
# 
# df.net.ko <- df.net.ko %>%
#     filter(pval < 0.05) %>%
#     group_by(target, pathway_name) %>%
#     summarise(interaction_count_ko = n()) %>%
#     filter(interaction_count_ko > 5)
# 
# df.net.wt <- df.net.wt %>%
#     filter(pval < 0.05) %>%
#     group_by(target, pathway_name) %>%
#     summarise(interaction_count_wt = n()) %>%
#     filter(interaction_count_wt > 5)
# 
# merged_df <- inner_join(df.net.noT, df.net.ko, by = c("target", "pathway_name")) %>%
#     inner_join(., df.net.wt, by = c("target", "pathway_name")) %>%
#     mutate(
#         diff_noT_wt = interaction_count_wt - interaction_count_noT,
#         diff_noT_ko = interaction_count_ko - interaction_count_noT,
#         diff_wt_ko = interaction_count_ko - interaction_count_wt
#     )
# 
# min_value <- min(c(merged_df$diff_noT_wt, merged_df$diff_noT_ko, merged_df$diff_wt_ko))
# max_value <- max(c(merged_df$diff_noT_wt, merged_df$diff_noT_ko, merged_df$diff_wt_ko))
# 
# #Heatmap noT vs wt
# heatmap.noT.wt <- merged_df %>%
#     ggplot(aes(x = target, y = pathway_name, fill = diff_noT_wt)) +
#     geom_tile() +
#     ggtitle("Wildtype vs. No T cells") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
#     scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value) )
# 
# #Heatmap noT vs ko
# heatmap.noT.ko <- merged_df %>%
#     ggplot(aes(x = target, y = pathway_name, fill = diff_noT_ko)) +
#     geom_tile() +
#     ggtitle("Knockout vs. No T cells") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
#     scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value))
# 
# #Heatmap wt vs ko
# heatmap.wt.ko <- merged_df %>%
#     ggplot(aes(x = target, y = pathway_name, fill = diff_wt_ko)) +
#     geom_tile() +
#     ggtitle("Knockout vs. Wildtype") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
#     scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value))
# 
# 
# grid.arrange(heatmap.noT.wt, heatmap.noT.ko, heatmap.wt.ko, ncol = 3)
# 
# 
# ######################################### ANALYSE WITH CELLCHAT ################################################
# 
# #### Set Pathways ####
# ######################
# pathways.show.noT <- c("LAMININ") 
# pathways.show.wt <- c("LAMININ") 
# pathways.show.ko <- c("LAMININ") 
# 
# ### Circle Plot: look at individual signalling pathways ###
# 
# cellchat_objects <- list(cellchat_noT, cellchat_wt, cellchat_ko)
# pathways <- list(pathways.show.noT, pathways.show.wt, pathways.show.ko)
# titles <- c("No T Cells", "Wildtype", "Knockout")
# 
# # Set up a 1x3 layout (1 row, 3 columns)
# par(mfrow = c(1, 3))
# 
# # Loop through the objects and create the plots and titles
# for (i in 1:3) {
#     plot <- netVisual_aggregate(cellchat_objects[[i]], signaling = pathways[[i]], layout = "circle")
#     title(main = titles[i])
# }
# 
# #Look at circle plot for each cell type 
# groupSize <- as.numeric(table(cellchat_wt@idents))
# mat <- cellchat_wt@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#     mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#     mat2[i, ] <- mat[i, ]
#     netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }
# 
# # Reset the layout to its default state 
# par(mfrow = c(1, 1))
# 
# ### Heatmaps of Pathway ###
# #noT
# par(mfrow=c(1,1))
# heatmap.noT <- netVisual_heatmap(cellchat_noT, signaling = pathways.show.noT, color.heatmap = "Reds")
# #Wildtype 
# par(mfrow=c(1,1))
# heatmap.wt <- netVisual_heatmap(cellchat_wt, signaling = pathways.show.noT, color.heatmap = "Reds")
# #Knockout 
# par(mfrow=c(1,1))
# heatmap.ko <- netVisual_heatmap(cellchat_ko, signaling = pathways.show.noT, color.heatmap = "Reds")
# 
# ### Contribution of L-R pairs to pathways ###
# #noT
# contribution_noT <- netAnalysis_contribution(cellchat_noT, signaling = pathways.show.noT, thresh = 0.05, return.data =  TRUE, title = "No T cell: LAMININ L-R Pairs")
# #Wildtype 
# contribution_wt <- netAnalysis_contribution(cellchat_wt, signaling = pathways.show.wt, thresh = 0.05, return.data =  TRUE, title = "Wildtype: LAMININ L-R Pairs")
# #Knockout 
# contribution_ko <- netAnalysis_contribution(cellchat_ko, signaling = pathways.show.ko, thresh = 0.05, return.data =  TRUE, title = "Wildtype: LAMININ L-R Pairs")
# 
# # Make data frame for graphing 
# nf <- merge(
#     contribution_noT$LR.contribution,
#     contribution_wt$LR.contribution,
#     suffixes=c(".noT", ".wt"),
#     by="name")
# 
# nf <- merge(nf, contribution_ko$LR.contribution, by="name") 
# 
# colnames(nf) <- c("names", "No T Cells", "Wildtype", "Knockout")
# 
# nf %>%
#     pivot_longer(!names, names_to = "condition", values_to = "contribution") %>%
#     group_by(condition) %>%
#     arrange(desc(contribution)) %>%
#     top_n(n= 15) %>%
#     ggplot(aes(x=contribution, y=names, fill=condition)) + geom_col(position="dodge") + 
#     theme_classic() +
#     scale_fill_manual(values=c("salmon3", 
#                                "navajowhite1", 
#                                "palegreen3"))  
# 
# 
# ### Violin Plot: Gene Expression ###
# pathways.show.noT <- c("LAMININ") 
# pathways.show.wt <- c("LAMININ") 
# pathways.show.ko <- c("LAMININ") 
# 
# plotGeneExpression(cellchat_noT, signaling = ("MHC-I"), enriched.only = TRUE, features = c("Klrd1")) + 
#     plotGeneExpression(cellchat_wt, signaling = ("MHC-I"), enriched.only = TRUE, features = c("Klrd1")) +
#     plotGeneExpression(cellchat_ko, signaling = ("MHC-I"), enriched.only = TRUE, features = c("Klrd1")) 
# 
