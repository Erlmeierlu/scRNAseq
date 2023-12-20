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
library(viridis)

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

doParallel::registerDoParallel(cores = 5)
foreach(organ = c('LN', 'Skin')) %dopar% {
    foreach(experiment = c('HDAC1', 'HDAC2')) %dopar% {
        foreach(treatment = c('NoT', 'HDAC_WT', 'HDAC_cKO')) %dopar% {
            
            subset <- monocle.obj[, monocle.obj$organ == organ & 
                                      monocle.obj$experiment == experiment &
                                      monocle.obj$treatment.agg == treatment]
            
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



# Exploratory Analsysis ---------------------------------------------------

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
cc_dat[, keyby = .(source, target, organ, experiment, treatment), 
       summarized_interactions := .N]


cc_dat %>%
    ggplot(aes(source, target, fill = summarized_interactions), na.value = 'lightgrey') +
    geom_tile() +
    facet_grid(organ + experiment ~ treatment, scales = 'free', space = 'free') +
    theme_my(panel.grid.major = element_blank(),
             panel.background = element_rect(fill = 'grey83')) +
    scale_fill_viridis(option = 'F')

ggsave(file.path(plotsDir, 'ligand_receptor_interactions_count.pdf'),
       width = 11, height = 9.5)


#Filter + calculate pathway interaction count difference between data frames 
cc_dat[, keyby = .(source, pathway_name, organ, experiment, treatment), source_pathway_count := .N]

cc_dat[, dcast(.SD, source + pathway_name + organ + experiment ~ treatment, 
               value.var = 'source_pathway_count',
               subset = .(source_pathway_count > 5))
       ][, .(source, 
            pathway_name, 
            organ, 
            experiment, 
            'WT_vs_NoT' = HDAC_WT - NoT,
            'cKO_vs_NoT' = HDAC_cKO - NoT,
            'cKO_vs_WT' = HDAC_cKO - HDAC_WT)
         ][,
           melt(.SD, measure.vars = patterns('vs'),
                variable.name = 'comparison',
                value.name = 'interaction_count_difference')
           ] %>% 
    ggplot() +
    geom_tile(aes(source, pathway_name, fill = interaction_count_difference)) +
    facet_grid(experiment + organ ~ comparison,
               scales = 'free', 
               space = 'free') +
    theme_my(panel.grid.major = element_blank(),
             panel.background = element_rect(fill = 'grey83')) +
    scale_fill_viridis(option = 'A')

ggsave(file.path(plotsDir, 'source_interactions_pathways_per_comparison.pdf'),
       width = 14, height = 17)


#Filter + calculate pathway interaction count difference between data frames 
cc_dat[, keyby = .(target, pathway_name, organ, experiment, treatment), target_pathway_count := .N]

cc_dat[, dcast(.SD, target + pathway_name + organ + experiment ~ treatment, 
               value.var = 'target_pathway_count',
               subset = .(target_pathway_count > 5))
      ][, .(target, 
            pathway_name, 
            organ, 
            experiment, 
            'WT_vs_NoT' = HDAC_WT - NoT,
            'cKO_vs_NoT' = HDAC_cKO - NoT,
            'cKO_vs_WT' = HDAC_cKO - HDAC_WT)
      ][,
        melt(.SD, measure.vars = patterns('vs'),
             variable.name = 'comparison',
             value.name = 'interaction_count_difference')
      ] %>% 
          ggplot() +
          geom_tile(aes(target, pathway_name, fill = interaction_count_difference)) +
          facet_grid(experiment + organ ~ comparison,
                     scales = 'free', 
                     space = 'free') +
          theme_my(panel.grid.major = element_blank(),
                   panel.background = element_rect(fill = 'grey83')) +
          scale_fill_viridis(option = 'A')

ggsave(file.path(plotsDir, 'target_interactions_pathways_per_comparison.pdf'),
       width = 14, height = 18)




# Cellchat Analysis -------------------------------------------------------


files <- list.files('/vscratch/scRNAseq/data/', pattern = 'cellchat')

cc_analysis <- function(object, pathways, thresh = 0.05){
  
    pathway_plots <- sapply(pathways,
                    function(x) {
                        
                        ###THIS PART IS FROM THE FUNCTION netVisual_aggregate
                        ###TO SKIP IF THERE IS NO SIG PATHWAY INSTEAD OF STOPPING
                        pairLR <- searchPair(signaling = x, 
                                             pairLR.use = object@LR$LRsig, 
                                             key = "pathway_name", 
                                             matching.exact = T, 
                                             pair.only = T)
                        net <- object@net 
                        pairLR.use.name <- dimnames(net$prob)[[3]]
                        pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
                        pairLR <- pairLR[pairLR.name, ]
                        prob <- net$prob
                        pval <- net$pval
                        prob[pval > thresh] <- 0
                        if (length(pairLR.name) > 1) {
                            pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
                        } else {
                            pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
                        }
                        ###UNTIL HERE. RETURN NULL INSTEAD OF STOPPING FUNCTION
                        
                        if (length(pairLR.name.use) == 0) return(NULL)
                        netVisual_aggregate(object, 
                                            signaling = x, 
                                            layout = "circle")
                    }, simplify = F)
    
    dev.off()
    
    groupSize <- as.numeric(table(object@idents))
    mat <- object@net$weight
    
    limit <- nrow(mat)
    
    dims <- c(ceiling(limit/6), 6)
    
    
    par(mfrow = dims, xpd=TRUE)
    for (i in 1:nrow(mat)) {
        mat2 <-
            matrix(
                0,
                nrow = nrow(mat),
                ncol = ncol(mat),
                dimnames = dimnames(mat)
            )
        mat2[i,] <- mat[i,]
        ct_plots <- netVisual_circle(
            mat2,
            vertex.weight = groupSize,
            weight.scale = T,
            edge.weight.max = max(mat),
            title.name = rownames(mat)[i]
        )
    }
    
    dev.off()
    
    par(mfrow=c(1,1))
    
    pathway_heatmaps <- sapply(pathways,
                               function(x) {
                                   if(!x %in% dimnames(slot(object, 'netP')$prob)[[3]]) return(NULL)
                                   netVisual_heatmap(object, 
                                                     signaling = x)
                               })
    dev.off()
    par(mfrow=c(1,1))
    contributions <- sapply(pathways,
                                 function(x) {
                                     ###THIS PART IS FROM THE FUNCTION netVisual_aggregate
                                     ###TO SKIP IF THERE IS NO SIG PATHWAY INSTEAD OF STOPPING
                                     pairLR <- searchPair(signaling = x, 
                                                          pairLR.use = object@LR$LRsig, 
                                                          key = "pathway_name", 
                                                          matching.exact = T, 
                                                          pair.only = T)
                                     net <- object@net 
                                     pairLR.use.name <- dimnames(net$prob)[[3]]
                                     pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
                                     pairLR <- pairLR[pairLR.name, ]
                                     prob <- net$prob
                                     pval <- net$pval
                                     prob[pval > thresh] <- 0
                                     if (length(pairLR.name) > 1) {
                                         pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
                                     } else {
                                         pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
                                     }
                                     ###UNTIL HERE. RETURN NULL INSTEAD OF STOPPING FUNCTION
                                     
                                     if (length(pairLR.name.use) == 0) return(NULL)
                                     netAnalysis_contribution(object, 
                                                              signaling = x, 
                                                              thresh = 0.05, 
                                                              return.data =  TRUE)
                                     
                                 }, simplify = F)
    
    dev.off()
    contribution_data <- lapply(contributions, function(x) x$LR.contribution)
    contribution_data <- rbindlist(contribution_data, idcol = 'pathway')
    
    contribution_plots <- lapply(contributions, function(x) x$gg.obj)
    
    list(pathway_plots = pathway_plots, 
         celltype_plots = ct_plots, 
         pathway_heatmaps = pathway_heatmaps,
         contribution_data = contribution_data,
         contribution_plots = contribution_plots)
}

cc <- sapply(files, function(x) read_rds(file.path('/vscratch/scRNAseq/data/', x)),
             simplify = F)


cc_res <- sapply(cc, function(x) {
    cc_analysis(x, 
                pathways = c('MHC-I', 'MHC-II', 'LAMININ', 'COLLAGEN')
                )},
    simplify = F)

contributions <- rbindlist(lapply(cc_res, function(x) x$contribution_data), idcol = 'object')

contributions[, c('organ', 'experiment', 'treatment') := tstrsplit(
  str_replace_all(
    str_remove_all(object, '^[:alpha:]+_|.[:alpha:]+$'), 
    'HDAC_', 'HDAC-'),
  '_')]

contributions[, c('treatment', 'object') := .(factor(
  str_replace(treatment, '-', '_'), levels = c('NoT', 'HDAC_WT', 'HDAC_cKO')),
  NULL)]

contributions[, by = .(organ, experiment, treatment), 
              .(rank = frank(-contribution, ties.method = 'dense'),
                pathway, name, contribution)
              ][
                rank < 10
                ] %>% 
  ggplot() +
  geom_col(aes(name, contribution, fill = treatment),
           position = 'dodge2') +
  facet_wrap(experiment ~ organ, scales = 'free') +
  theme_my() + 
  scale_fill_manual(values=c("salmon3",
                              "navajowhite1",
                              "palegreen3"))

ggsave(file.path(plotsDir, 'contributions.pdf'),
       width = 10, height = 7)

#save all the plots
lapply(seq_along(cc_res), function(x, a, b){
  name <- str_remove_all(b[x], '^cellchat_|.rds$')
  
  lapply(seq_along(a[[x]]$pathway_plots), function(y, n, m){
    
    if(is.null(n[[y]])) return(NULL)
    pdf(file.path(plotsDir, paste(m[y], name, 'pathway_plot.pdf', sep = '_')))
    print(n[[y]])
    dev.off()
    
  }, n = a[[x]]$pathway_plots, m = names(a[[x]]$pathway_plots))
  
  
  pdf(file.path(plotsDir, paste(name, 'celltype_plot.pdf', sep = '_')),
      width = 16, height = 14)
  print(a[[x]]$celltype_plots)
  dev.off()
  
  lapply(seq_along(a[[x]]$pathway_heatmaps), function(i, c, d){
    
    if(is.null(c[[i]])) return(NULL)
    pdf(file.path(plotsDir, paste(d[i], name, 'pathway_heatmap.pdf', sep = '_')))
    print(c[[i]])
    dev.off()
  }, c = a[[x]]$pathway_heatmaps, d = names(a[[x]]$pathway_heatmaps))
  
  lapply(seq_along(a[[x]]$contribution_plots), function(j, e, f){
    
    if(is.null(e[[j]])) return(NULL)
    pdf(file.path(plotsDir, paste(f[j], name, 'contribution_plots.pdf', sep = '_')))
    print(e[[j]])
    dev.off()
    
  }, e = a[[x]]$contribution_plots, f = names(a[[x]]$contribution_plots))
  
}, a = cc_res, b = names(cc_res))
