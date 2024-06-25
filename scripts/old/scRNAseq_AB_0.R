library(renv)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(hexbin)
library(patchwork)

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
        panel.spacing.x = unit(10,"mm"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.key = element_blank(),
        text = element_text(size = 15),
        strip.text.x = element_text(size = 10, margin = margin(b = 2, t = 2)),
        strip.background = element_rect(fill = "#9FD7D2", colour = "black", size = 1)
    )
}

#Functions
make_col_stats <- function(counts.df) {
    x <- counts.df
    hash <- x[grepl("-[0-9]+$", rownames(x)),]
    cite45 <- x[grepl("CD45.*$", rownames(x)),]
    cite451 <- x[grepl("CD45.1$", rownames(x)),]
    cite4 <- x[grepl("CD4$", rownames(x)),]
    cite48 <- x[grepl("CD[1-9]$", rownames(x)),]
    #gdt <- x[grepl("GD", rownames(x)),]
    
    hashtags <-
        data.frame(sum = colSums(hash),
                   max = colMaxs(hash),
                   ratio = colMaxs(hash) / colSums(hash))
    
    cd45.1 <- 
        data.frame(sum = colSums(cite45),
                   cd45.1 = cite451,
                   ratio = cite451 / colSums(cite45))
    
    cd4 <-
        data.frame(sum = colSums(cite48),
                   cd4 = cite4,
                   ratio = cite4 / colSums(cite48))
    
    #gd <- 
    #data.frame(sum = gdt)
    
    return(list(hashtags = hashtags, cd45x = cd45.1, cd4_cd8 = cd4
                #, gdt = gd
    ))
}

make_row_stats <- function(counts.df) {
    x <- counts.df
    hash <- x[grepl("-[0-9]+$", rownames(x)),]
    cite45 <- x[grepl("CD45.*$", rownames(x)),]
    cite48 <- x[grepl("CD[1-9]$", rownames(x)),]
    gdt <- x[grepl("GD", rownames(x)),]
    
    hashtags <-
        data.frame(sum = rowSums(hash),
                   max = rowMaxs(hash),
                   ratio = rowMaxs(hash) / rowSums(hash),
                   cells = ncol(x),
                   rpc = rowSums(hash)/ncol(x))
    
    cd45x <- 
        data.frame(sum = rowSums(cite45),
                   max = rowMaxs(cite45),
                   ratio = rowMaxs(cite45) / rowSums(cite45),
                   cells = ncol(x),
                   rpc = rowSums(cite45)/ncol(x))
    
    cd4_cd8 <-
        data.frame(sum = rowSums(cite48),
                   max = rowMaxs(cite48),
                   ratio = rowMaxs(cite48) / rowSums(cite48),
                   cells = ncol(x),
                   rpc = rowSums(cite48)/ncol(x))
    
    gd <- 
        data.frame(sum = sum(gdt, na.rm = T),
                   max = ifelse(length(gdt) != 0, max(gdt, na.rm = T), 0),
                   ratio = ifelse(length(gdt) != 0, max(gdt, na.rm = T) / sum(gdt, na.rm = T), 0),
                   cells = ncol(x),
                   rpc = sum(gdt, na.rm = T)/ncol(x))
    
    
    return(list(hashtags = hashtags, cd45x = cd45x, cd4_cd8 = cd4_cd8, gdt = gd))
}

make_CITE_stats <- function(counts.df) {
    x <- as.matrix(counts.df)
    
    ratios <- 
        data.frame(cd45.1_.2_ratio = (x[grepl("CD45.1$", rownames(x)),]+1)/(x[grepl("CD45.2$", rownames(x)),]+1),
                   cd4_8_ratio = (x[grepl("CD4$", rownames(x)),]+1)/(x[grepl("CD8$", rownames(x)),]+1))
    
    return(ratios)
}



# Load Data ---------------------------------------------------------------


baseDir <- getwd()
rawDir <- ("/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin/raw_data/scRNA_from_BSF/COUNT")
plotsDir <- file.path(baseDir, "plots/")
tablesDir <- file.path(baseDir, "tables/")
dataDir <- file.path(baseDir, "data/")


CountsAB <- readRDS(file.path(dataDir, "ab_counts.rds"))


# AB Distri ---------------------------------------------------------------

# Raw Quality Control ---------------------------------------------------------
ab_col_statlist <- list()
ab_row_statlist <- list()

for (experiment in names(CountsAB)){
    ab_col_statlist[[experiment]] <- lapply(CountsAB[[experiment]], make_col_stats)
    ab_row_statlist[[experiment]] <- lapply(CountsAB[[experiment]], make_row_stats)
}


ab_col_stats <-  bind_rows(lapply(ab_col_statlist, function(l1){
    bind_rows(
        lapply(l1, function(l2){
            bind_rows(l2, .id = 'abtype')
        }), .id = 'sample')
}), .id = 'experiment')

ab_row_stats <- bind_rows(lapply(ab_row_statlist, function(l1){
    bind_rows(
        lapply(l1, function(l2){
            bind_rows(l2, .id = 'abtype')
        }), .id = 'sample')
}), .id = 'experiment')

CITE_statlist <- list()
for (experiment in names(CountsAB)){
    CITE_statlist[[experiment]] <- lapply(CountsAB[[experiment]], make_CITE_stats)
}
CITE_stats <- bind_rows(lapply(CITE_statlist, bind_rows, .id = "sample"), .id = "experiment")

##All ABs ---------------------------------------------------------------

ab_row_stats %>% group_by(abtype, sample, experiment) %>% summarise(sum = sum(sum)) %>% 
    ggplot() +
    geom_col(aes(abtype, log10(sum+1), fill = sample), col = "black", position = "dodge") + 
    facet_wrap(~ experiment, scales = "free") + 
    theme_my() +
    labs(x = "Antibody Type", y = "log10(Reads+1)", title = "Antibody reads per sample")

ggsave(file.path(plotsDir, "AB_reads_per_sample.pdf"))

ab_row_stats %>% group_by(abtype, sample, experiment, cells) %>% summarise(sum = sum(sum)) %>% mutate(rpc = sum/cells) %>% 
    ggplot() +
    geom_col(aes(abtype, log10(rpc+1), fill = sample), col = "black", position = "dodge") + 
    facet_wrap(~ experiment, scales = "free") + 
    theme_my() +
    labs(x = "Antibody Type", y = "log10(rpc+1)", title = "Antibody reads per cell")

ggsave(file.path(plotsDir, "AB_reads_per_cell.pdf"))


##Hashtag ABs ------------------------------------------------------------
hashtag_plot_list <- list()
for (experimentx in unique(ab_col_stats$experiment)) {
    hashtag_plot_list[[experimentx]] <-
        ab_col_stats %>% filter(abtype == "hashtags", experiment == experimentx) %>%  group_by(sum, sample) %>% dplyr::count() %>%
        ggplot() +
        geom_col(aes(log10(sum + 1), n), col = "gray25", fill = "ivory2") +
        facet_wrap( ~ sample, scales = "free", ncol = 3) +
        theme_my() +
        labs(
            title = paste("Hashtag reads per sample", experimentx),
            x = "log10(Reads + 1)",
            y = "Number of Cells"
        )
}

for (plot in names(hashtag_plot_list)){
    ggsave(file.path(plotsDir, paste0("Hashtag_reads_per_sample_",plot,".pdf")), plot = hashtag_plot_list[[plot]])
}

ab_row_stats %>% filter(abtype == "hashtags") %>% 
    ggplot() +
    geom_col(aes(purrr::map_chr(
        rownames(ab_row_stats %>% filter(abtype == "hashtags")) %>% str_extract_all("^HTO|[0-9]+$"),
        ~ str_c(.x, collapse = "-")
    ), log10(sum),
    fill = sample),
    col = "gray25", position = "dodge") +
    facet_wrap(~ experiment, scales = "free") +
    theme_my() +
    theme(axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1,
        size = 10
    )) +
    labs(title = "Hashtag subtype distribution", x = "Hashtag subtype", y = "log10(Reads)")

ggsave(file.path(plotsDir, "Hashtag_subtype_distribution.pdf"))


hashtag_ratio_plot_list <- list()

for (experimentx in unique(ab_col_stats$experiment)) {
    
    hashtag_ratio_plot_list[[experimentx]] <- ab_col_stats %>% filter(abtype == "hashtags", sum != 0, experiment == experimentx) %>% 
        ggplot() +
        stat_binhex(aes(sum, ratio, fill = log(..count..)), bins = 30, col = "black") +
        facet_wrap( ~ sample, scales = "free_x", ncol = 3) +
        scale_fill_gradient(low = "ivory2", high = "red") +
        geom_hline(yintercept = 0.6, linetype = "dashed", colour = "blue")+
        scale_x_log10() +
        theme_my() +
        labs(x = "Reads", y = "Ratio [MaxReads/SumReads]", title = paste("Hashtag ratio per sample", experimentx))
    
}

for (plot in names(hashtag_ratio_plot_list)){
    ggsave(file.path(plotsDir, paste0("Hashtag_ratio_per_sample_",plot,".pdf")), plot = hashtag_ratio_plot_list[[plot]])
}


hashtag_distribution_plot_list <- list()
for (experimentx in unique(ab_row_stats$experiment)) {
    hashtag_distribution_plot_list[[experimentx]] <- ab_row_stats %>% filter(abtype == "hashtags", sum != 0, experiment == experimentx) %>% 
        ggplot() +
        geom_point(aes(sum, ratio), shape = 21, alpha = 0.5, col = "black", size = 2, fill = "ivory2") +
        facet_wrap( ~ sample, scales = "free_x") +
        theme_my() +
        labs(title = paste("Distribution of reads per hashtag", experimentx), x = "Reads per hashtag", y = "Ratio [MaxReads/SumReads]")
}

for(plot in names(hashtag_distribution_plot_list)){
    ggsave(file.path(plotsDir, paste0("Distribution_of_reads_per_hashtag",plot,".pdf")), plot = hashtag_distribution_plot_list[[plot]])
}

##CD45 ABs ---------------------------------------------------------------
cd45_plot_list <- list()
for (experimentx in unique(ab_col_stats$experiment)){
    cd45_plot_list[[experimentx]] <- ab_col_stats %>%
        filter(abtype == "cd45x", experiment == experimentx) %>%
        group_by(sample, sum == 0, ratio %in% c(0, 1), sum != 0 & !(ratio %in% c(0, 1))) %>%
        summarise(
            sum = sum(`sum == 0`),
            ratio = sum(`ratio %in% c(0, 1)`),
            both = sum(`sum != 0 & !(ratio %in% c(0, 1))`)
        ) %>%
        mutate(
            count = sum(c(sum, ratio, both), na.rm = T),
            ratio = replace_na(ratio, 0),
            Distinct_Tags = case_when(
                sum == 0 &
                    ratio == 0 ~ "2",
                sum == 0 & both == 0 ~ "1",
                ratio == 0 & both == 0 ~ "0"
            )
        ) %>%
        ggplot() +
        geom_col(aes(Distinct_Tags, log10(count))) +
        facet_wrap(~ sample, scales = "free_y") +
        labs(x = "Number of distinct tags", title = paste("Distinct CD45.1 / CD45.2 tags per sample", experimentx)) +
        theme_my()
}
for (plot in names(cd45_plot_list)){
    ggsave(file.path(plotsDir, paste0("Distinct_CD45_tags_per_sample_",plot,".pdf")), plot = cd45_plot_list[[plot]])
}


cd45_ratio_plot_list <- list()

for (experimentx in unique(ab_col_stats$experiment)) {
    cd45_ratio_plot_list[[experimentx]] <- ab_col_stats %>% filter(abtype == "cd45x", !(ratio %in% c(NA, NaN)), experiment == experimentx) %>% 
        ggplot() + 
        geom_violin(aes(sample, ratio, fill = sample)) +
        facet_wrap(~experiment, scales = "free_y") +
        theme_my() +
        labs(title = paste("Ratio of CD45.1 / CD45.x reads", experimentx), y = "Ratio [CD45.1 Reads/ Sum Reads)]")
}

for(plot in names(cd45_ratio_plot_list)){
    ggsave(file.path(plotsDir, paste0("Ratio_of_CD45_reads", plot,".pdf")), plot = cd45_ratio_plot_list[[plot]])
}


cd45_reads_plot_list <- list()

for (experimentx in unique(ab_col_stats$experiment)) {
    cd45_reads_plot_list[[experimentx]] <- ab_col_stats %>% filter(abtype == "cd45x", experiment == experimentx) %>%  group_by(sum, sample) %>% dplyr::count() %>% 
        ggplot() +
        geom_col(aes(log10(sum+1), n), col = "gray25", fill = "ivory2") +
        facet_wrap(~ sample, scales = "free", ncol = 3) +
        theme_my() +
        labs(title = paste("CD45.x reads per sample", experimentx), x = "log10(Reads + 1)", y = "Number of Cells")
}

for(plot in names(cd45_reads_plot_list)){
    ggsave(file.path(plotsDir, paste0("CD45_reads_per_sample_",plot,".pdf")), plot = cd45_reads_plot_list[[plot]])
}

ab_row_stats %>% filter(abtype == "cd45x") %>% 
    ggplot() +
    geom_col(aes(rownames(ab_row_stats %>% filter(abtype == "cd45x")) %>% str_extract("CD45.[1-2]"), log10(sum), fill = sample),
             col = "gray25",
             position = "dodge") +
    facet_wrap(~ experiment, scales = "free_y") +
    theme_my() +
    theme(axis.text.x = element_text(
        size = 10
    )) +
    labs(title = "CD45.1 / CD45.2 subtype distribution", x = "CITE subtype", y = "log10(Reads)")

ggsave(file.path(plotsDir, "CD45_subtype_distribution.pdf"))

cd45_hex_plot_list <- list()
for (experimentx in unique(ab_col_stats$experiment)){
    cd45_hex_plot_list[[experimentx]] <- ab_col_stats %>% filter(abtype == "cd45x", sum != 0, experiment == experimentx) %>% 
        ggplot() +
        stat_binhex(aes(log10(sum + 1), ratio, fill = log(..count..)), bins = 30, col = "black") +
        facet_wrap( ~ sample, scales = "free_x", ncol = 3) +
        scale_fill_gradient(low = "ivory2", high = "red") +
        theme_my() +
        labs(x = "log10(Reads + 1)", title = paste("CD45.1 / CD45.x ratio per sample", experimentx), y = "Ratio [CD45.1 Reads/ Sum Reads)]")
}

for(plot in names(cd45_hex_plot_list)){
    ggsave(file.path(plotsDir, paste0("CD45_ratio_per_sample",plot,".pdf")), plot= cd45_hex_plot_list[[plot]])
}


##CD8/CD4 ABs -------------------------------------------------------------
cd4_plot_list <- list()
for (experimentx in unique(ab_col_stats$experiment)){
    cd4_plot_list[[experimentx]] <- ab_col_stats %>%
        filter(abtype == "cd4_cd8", experiment == experimentx) %>%
        group_by(sample, sum == 0, ratio %in% c(0, 1), sum != 0 & !(ratio %in% c(0, 1))) %>%
        summarise(
            sum = sum(`sum == 0`),
            ratio = sum(`ratio %in% c(0, 1)`),
            both = sum(`sum != 0 & !(ratio %in% c(0, 1))`)
        ) %>%
        mutate(
            count = sum(c(sum, ratio, both), na.rm = T),
            ratio = replace_na(ratio, 0),
            Distinct_Tags = case_when(
                sum == 0 &
                    ratio == 0 ~ "2",
                sum == 0 & both == 0 ~ "1",
                ratio == 0 & both == 0 ~ "0"
            )
        ) %>%
        ggplot() +
        geom_col(aes(Distinct_Tags, log10(count))) +
        facet_wrap(~ sample, scales = "free_y") +
        labs(x = "Number of distinct tags", title = paste("Distinct CD4 / CD8 tags per sample", experimentx)) +
        theme_my()
}
for(plot in names(cd4_plot_list)){
    ggsave(file.path(plotsDir, paste0("Distinct_CD4_tags_per_sample_",plot,".pdf")), plot = cd4_plot_list[[plot]])
}

cd4_ratio_plot_list <- list()

for (experimentx in unique(ab_col_stats$experiment)) {
    cd4_ratio_plot_list[[experimentx]] <- ab_col_stats %>% filter(abtype == "cd4_cd8", !(ratio %in% c(NA, NaN)), experiment == experimentx) %>% 
        ggplot() + 
        geom_violin(aes(sample, ratio, fill = sample)) +
        facet_wrap(~experiment, scales = "free_y") +
        theme_my() +
        labs(title = paste("Ratio of CD4 / CDx reads", experimentx), y = "Ratio [CD4 Reads / Sum Reads)]")
}
for (plot in names(cd4_ratio_plot_list)){
    ggsave(file.path(plotsDir, paste0("Ratio_of_CD4_reads_",plot,".pdf")), plot = cd4_ratio_plot_list[[plot]])
}

cd4_reads_plot_list <- list()
for (experimentx in unique(ab_col_stats$experiment)){
    cd4_reads_plot_list[[experimentx]] <- ab_col_stats %>% filter(abtype == "cd4_cd8", experiment == experimentx) %>%  group_by(sum, sample) %>% dplyr::count() %>% 
        ggplot() +
        geom_col(aes(log10(sum + 1), n), col = "gray25", fill = "ivory2") +
        facet_wrap(~ sample, scales = "free") +
        theme_my() +
        labs(title = paste("CD4 / CD8 reads per sample", experimentx), x = "log10(Reads + 1)", y = "Number of Cells")
}
for(plot in names(cd4_reads_plot_list)){
    ggsave(file.path(plotsDir, paste0("CD4_reads_per_sample_",plot,".pdf")), plot = cd4_reads_plot_list[[plot]])
}


ab_row_stats %>% filter(abtype == "cd4_cd8") %>% 
    ggplot() +
    geom_col(aes(rownames(ab_row_stats %>% filter(abtype == "cd4_cd8")) %>% str_extract("CD[0-9]"), log10(sum), fill = sample),
             col = "gray25",
             position = "dodge") +
    facet_wrap(~ experiment, scales = "free_y") +
    theme_my() +
    theme(axis.text.x = element_text(
        size = 10
    )) +
    labs(title = "CD4 / CD8 subtype distribution", x = "CITE subtype", y = "log10(Reads)")

ggsave(file.path(plotsDir, "CD4_subtype_distribution.pdf"))

cd4_hex_plot_list <- list()
for (experimentx in unique(ab_col_stats$experiment)){
    cd4_hex_plot_list[[experimentx]] <- ab_col_stats %>% filter(abtype == "cd4_cd8", sum != 0, experiment == experimentx) %>% 
        ggplot() +
        stat_binhex(aes(log(sum + 1), ratio,fill = log(..count..)), bins = 30, col = "black") +
        facet_wrap( ~ sample, scales = "free_x", ncol = 3) +
        scale_fill_gradient(low = "ivory2", high = "red") +
        theme_my() +
        labs(x = "log10(Reads +1)", title = paste("CD4 / CDx ratio per sample", experimentx), y = "Ratio [CD4 Reads / Sum Reads)]")
}
for(plot in names(cd4_hex_plot_list)){
    ggsave(file.path(plotsDir, paste0("CD4_ratio_per_sample_",plot,".pdf")), plot = cd4_hex_plot_list[[plot]])
}


##CD45.1/.2 - CD4/8 combinations ------------------------------------------
CITE_plot_list <- list()
for (experimentx in unique(CITE_stats$experiment)){
    CITE_plot_list[[experimentx]] <- CITE_stats %>% filter(experiment == experimentx) %>% 
        ggplot() + 
        stat_bin_hex(aes(log2(cd45.1_.2_ratio + 0.01), log2(cd4_8_ratio + 0.01), fill = log(..count..)), bins = 30, col = "black") + 
        facet_wrap(~sample, scales = "fixed", ncol = 3) + 
        scale_fill_gradient(low = "ivory2", high = "red") +
        theme_my() +
        labs(x = "log2(Ratio [CD45.1:CD45.2] + 0.01)", y = "log10(Ratio [CD4:CD8] + 0.01)", title = paste("Combination of CITE-tags", experimentx))
}
for(plot in names(CITE_plot_list)){
    ggsave(file.path(plotsDir, paste0("CITE_tag_combinations_",plot,".pdf")), plot = CITE_plot_list[[plot]])
}
