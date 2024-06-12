#setting up directories
gfsDir <- '/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin'
rawDir <- file.path(gfsDir, 'raw_data/scRNA_from_BSF/COUNT')
plotsDir <- file.path(gfsDir, 'plots')
tablesDir <- file.path(gfsDir, 'tables')
shinyDir <- 'dge-app'
dataDir <-'data'
resDir <- 'results'
vDir <- '/vscratch/scRNAseq/data'
oldDir <- file.path(vDir, 'old')

#ggplot theme
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