library(monocle3)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(CellChat)
library(Seurat)
library(readr)

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
vDir <- ("/vscratch/scRNAseq")
plotsDir <- ("/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin/plots")
tablesDir <- file.path(vDir, "tables")
oldDir <- file.path(vDir, "data/old")
shinyDir <- ('dge-app')
dataDir <-("data")
resDir <- ("results")

# Load Data ---------------------------------------------------------------
monocle.obj <- read_rds(file.path(dataDir, '3_annotated_monocle.cds'))

# Create CellChat Object --------------------------------------------------

ctrl_subset <- monocle.obj[, monocle.obj@colData %>% subset(
    organ == 'Skin' & experiment == 'HDAC2' & treatment.agg == 'NoT'
) %>% rownames]

countsx <- counts(ctrl_subset)
countsx <- normalizeData(countsx)
meta <- data.frame(ctrl_subset@colData)

ctrl_subset <- createCellChat(
    countsx,
    meta = meta,
    group.by = "celltype"
)

# Data Prep
ctrl_subset@DB <- CellChatDB.mouse

ctrl_subset <- subsetData(ctrl_subset)

future::plan("multisession", workers = 4) # do parallel
ctrl_subset <- identifyOverExpressedGenes(ctrl_subset)
ctrl_subset <- identifyOverExpressedInteractions(ctrl_subset)
ctrl_subset


cellchat_noT <- subsetData(cellchat_noT) 
cellchat_noT <- identifyOverExpressedGenes(cellchat_noT)
cellchat_noT <- identifyOverExpressedInteractions(cellchat_noT)
cellchat_noT <- computeCommunProb(cellchat_noT)
cellchat_noT  <- filterCommunication(cellchat_noT , min.cells = 10) 
cellchat_noT <- computeCommunProbPathway(cellchat_noT)
cellchat_noT <- aggregateNet(cellchat_noT)
cellchat_noT <- netAnalysis_computeCentrality(cellchat_noT)
df.net.noT <- subsetCommunication(cellchat_noT)

### WILDTYPE ###
# Data Prep
cellchat_wt@DB <- CellChatDB.mouse
cellchat_wt <- subsetData(cellchat_wt) 
cellchat_wt <- identifyOverExpressedGenes(cellchat_wt)
cellchat_wt <- identifyOverExpressedInteractions(cellchat_wt)
cellchat_wt <- computeCommunProb(cellchat_wt)
cellchat_wt  <- filterCommunication(cellchat_wt , min.cells = 10) 
cellchat_wt <- computeCommunProbPathway(cellchat_wt)
cellchat_wt <- aggregateNet(cellchat_wt)
cellchat_wt <- netAnalysis_computeCentrality(cellchat_wt)
df.net.wt <- subsetCommunication(cellchat_wt)


### KNOCKOUT ###
# Data Prep
cellchat_ko@DB <- CellChatDB.mouse
cellchat_ko <- subsetData(cellchat_ko) 
cellchat_ko <- identifyOverExpressedGenes(cellchat_ko)
cellchat_ko <- identifyOverExpressedInteractions(cellchat_ko)
cellchat_ko <- computeCommunProb(cellchat_ko)
cellchat_ko  <- filterCommunication(cellchat_ko , min.cells = 10) 
cellchat_ko <- computeCommunProbPathway(cellchat_ko)
cellchat_ko <- aggregateNet(cellchat_ko)
cellchat_ko <- netAnalysis_computeCentrality(cellchat_ko)
df.net.ko <- subsetCommunication(cellchat_ko)

#########################################  ANALYSE  ################################################

### CELL-CELL INTERACTIONS HEAT MAP ###

df.net.noT %>%
    filter(pval < 0.05) %>%
    group_by(source, target) %>%
    count() %>%
    ggplot(aes(x=source, y=target, fill=n)) + geom_tile() +
    ggtitle("No T Cells") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) +
    scale_fill_gradient(low = "navy",
                        high = "springgreen",
                        guide = "colorbar",
                        limits = c(0, 100)) +
    df.net.wt %>%
    filter(pval < 0.05) %>%
    group_by(source, target) %>%
    count() %>%
    ggplot(aes(x=source, y=target, fill=n)) + geom_tile() +
    ggtitle("Wildtype") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
    scale_fill_gradient(low = "navy",
                        high = "springgreen",
                        guide = "colorbar",
                        limits = c(0, 100)) +
    df.net.ko %>%
    filter(pval < 0.05) %>%
    group_by(source, target) %>%
    count() %>%
    ggplot(aes(x=source, y=target, fill=n)) + geom_tile() +
    ggtitle("Knockout") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
    scale_fill_gradient(low = "navy",
                        high = "springgreen",
                        guide = "colorbar",
                        limits = c(0, 100)) 


### HEATMAP: Difference in pathway counts of all treatments - Source Cells ###
#Filter + calculate pathway interaction count difference between data frames 

df.net.noT <- df.net.noT %>%
    filter(pval < 0.05) %>%
    group_by(source, pathway_name) %>%
    summarise(interaction_count_noT = n()) %>%
    filter(interaction_count_noT > 5) 

df.net.ko <- df.net.ko %>%
    filter(pval < 0.05) %>%
    group_by(source, pathway_name) %>%
    summarise(interaction_count_ko = n()) %>%
    filter(interaction_count_ko > 5)

df.net.wt <- df.net.wt %>%
    filter(pval < 0.05) %>%
    group_by(source, pathway_name) %>%
    summarise(interaction_count_wt = n()) %>%
    filter(interaction_count_wt > 5)

merged_df <- inner_join(df.net.noT, df.net.ko, by = c("source", "pathway_name")) %>%
    inner_join(., df.net.wt, by = c("source", "pathway_name")) %>%
    mutate(
        diff_noT_wt = interaction_count_wt - interaction_count_noT,
        diff_noT_ko = interaction_count_ko - interaction_count_noT,
        diff_wt_ko = interaction_count_ko - interaction_count_wt
    )

min_value <- min(c(merged_df$diff_noT_wt, merged_df$diff_noT_ko, merged_df$diff_wt_ko))
max_value <- max(c(merged_df$diff_noT_wt, merged_df$diff_noT_ko, merged_df$diff_wt_ko))

#Heatmap noT vs wt
heatmap.noT.wt <- merged_df %>%
    ggplot(aes(x = source, y = pathway_name, fill = diff_noT_wt)) +
    geom_tile() +
    ggtitle("Wildtype vs. No T cells") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
    scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value))

#Heatmap noT vs ko
heatmap.noT.ko <- merged_df %>%
    ggplot(aes(x = source, y = pathway_name, fill = diff_noT_ko)) +
    geom_tile() +
    ggtitle("Knockout vs. No T cells") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
    scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value))

#Heatmap wt vs ko
heatmap.wt.ko <- merged_df %>%
    ggplot(aes(x = source, y = pathway_name, fill = diff_wt_ko)) +
    geom_tile() +
    ggtitle("Knockout vs. Wildtype") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
    scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value))


grid.arrange(heatmap.noT.wt, heatmap.noT.ko, heatmap.wt.ko, ncol = 3)

### HEATMAP: Difference in pathway counts of all treatments - Target Cells ###

df.net.noT <- df.net.noT %>%
    filter(pval < 0.05) %>%
    group_by(target, pathway_name) %>%
    summarise(interaction_count_noT = n()) %>%
    filter(interaction_count_noT > 5)

df.net.ko <- df.net.ko %>%
    filter(pval < 0.05) %>%
    group_by(target, pathway_name) %>%
    summarise(interaction_count_ko = n()) %>%
    filter(interaction_count_ko > 5)

df.net.wt <- df.net.wt %>%
    filter(pval < 0.05) %>%
    group_by(target, pathway_name) %>%
    summarise(interaction_count_wt = n()) %>%
    filter(interaction_count_wt > 5)

merged_df <- inner_join(df.net.noT, df.net.ko, by = c("target", "pathway_name")) %>%
    inner_join(., df.net.wt, by = c("target", "pathway_name")) %>%
    mutate(
        diff_noT_wt = interaction_count_wt - interaction_count_noT,
        diff_noT_ko = interaction_count_ko - interaction_count_noT,
        diff_wt_ko = interaction_count_ko - interaction_count_wt
    )

min_value <- min(c(merged_df$diff_noT_wt, merged_df$diff_noT_ko, merged_df$diff_wt_ko))
max_value <- max(c(merged_df$diff_noT_wt, merged_df$diff_noT_ko, merged_df$diff_wt_ko))

#Heatmap noT vs wt
heatmap.noT.wt <- merged_df %>%
    ggplot(aes(x = target, y = pathway_name, fill = diff_noT_wt)) +
    geom_tile() +
    ggtitle("Wildtype vs. No T cells") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
    scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value) )

#Heatmap noT vs ko
heatmap.noT.ko <- merged_df %>%
    ggplot(aes(x = target, y = pathway_name, fill = diff_noT_ko)) +
    geom_tile() +
    ggtitle("Knockout vs. No T cells") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
    scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value))

#Heatmap wt vs ko
heatmap.wt.ko <- merged_df %>%
    ggplot(aes(x = target, y = pathway_name, fill = diff_wt_ko)) +
    geom_tile() +
    ggtitle("Knockout vs. Wildtype") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title = element_blank(), legend.title = element_blank()) + 
    scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = c(min_value, max_value))


grid.arrange(heatmap.noT.wt, heatmap.noT.ko, heatmap.wt.ko, ncol = 3)


######################################### ANALYSE WITH CELLCHAT ################################################

#### Set Pathways ####
######################
pathways.show.noT <- c("LAMININ") 
pathways.show.wt <- c("LAMININ") 
pathways.show.ko <- c("LAMININ") 

### Circle Plot: look at individual signalling pathways ###

cellchat_objects <- list(cellchat_noT, cellchat_wt, cellchat_ko)
pathways <- list(pathways.show.noT, pathways.show.wt, pathways.show.ko)
titles <- c("No T Cells", "Wildtype", "Knockout")

# Set up a 1x3 layout (1 row, 3 columns)
par(mfrow = c(1, 3))

# Loop through the objects and create the plots and titles
for (i in 1:3) {
    plot <- netVisual_aggregate(cellchat_objects[[i]], signaling = pathways[[i]], layout = "circle")
    title(main = titles[i])
}

#Look at circle plot for each cell type 
groupSize <- as.numeric(table(cellchat_wt@idents))
mat <- cellchat_wt@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Reset the layout to its default state 
par(mfrow = c(1, 1))

### Heatmaps of Pathway ###
#noT
par(mfrow=c(1,1))
heatmap.noT <- netVisual_heatmap(cellchat_noT, signaling = pathways.show.noT, color.heatmap = "Reds")
#Wildtype 
par(mfrow=c(1,1))
heatmap.wt <- netVisual_heatmap(cellchat_wt, signaling = pathways.show.noT, color.heatmap = "Reds")
#Knockout 
par(mfrow=c(1,1))
heatmap.ko <- netVisual_heatmap(cellchat_ko, signaling = pathways.show.noT, color.heatmap = "Reds")

### Contribution of L-R pairs to pathways ###
#noT
contribution_noT <- netAnalysis_contribution(cellchat_noT, signaling = pathways.show.noT, thresh = 0.05, return.data =  TRUE, title = "No T cell: LAMININ L-R Pairs")
#Wildtype 
contribution_wt <- netAnalysis_contribution(cellchat_wt, signaling = pathways.show.wt, thresh = 0.05, return.data =  TRUE, title = "Wildtype: LAMININ L-R Pairs")
#Knockout 
contribution_ko <- netAnalysis_contribution(cellchat_ko, signaling = pathways.show.ko, thresh = 0.05, return.data =  TRUE, title = "Wildtype: LAMININ L-R Pairs")

# Make data frame for graphing 
nf <- merge(
    contribution_noT$LR.contribution,
    contribution_wt$LR.contribution,
    suffixes=c(".noT", ".wt"),
    by="name")

nf <- merge(nf, contribution_ko$LR.contribution, by="name") 

colnames(nf) <- c("names", "No T Cells", "Wildtype", "Knockout")

nf %>%
    pivot_longer(!names, names_to = "condition", values_to = "contribution") %>%
    group_by(condition) %>%
    arrange(desc(contribution)) %>%
    top_n(n= 15) %>%
    ggplot(aes(x=contribution, y=names, fill=condition)) + geom_col(position="dodge") + 
    theme_classic() +
    scale_fill_manual(values=c("salmon3", 
                               "navajowhite1", 
                               "palegreen3"))  


### Violin Plot: Gene Expression ###
pathways.show.noT <- c("LAMININ") 
pathways.show.wt <- c("LAMININ") 
pathways.show.ko <- c("LAMININ") 

plotGeneExpression(cellchat_noT, signaling = ("MHC-I"), enriched.only = TRUE, features = c("Klrd1")) + 
    plotGeneExpression(cellchat_wt, signaling = ("MHC-I"), enriched.only = TRUE, features = c("Klrd1")) +
    plotGeneExpression(cellchat_ko, signaling = ("MHC-I"), enriched.only = TRUE, features = c("Klrd1")) 

