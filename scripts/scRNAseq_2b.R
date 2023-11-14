library(tidyr)
library(dplyr)
library(stringr)
library(monocle3)
library(data.table)

#setting up directories
baseDir <- getwd()
rawDir <- ("/media/AGFORTELNY/PROJECTS/Gratz_InflammedSkin/raw_data/scRNA_from_BSF/COUNT")
plotsDir <- file.path(baseDir, "plots/")
tablesDir <- file.path(baseDir, "tables/")
dataDir <- file.path(baseDir, "data/")

# Load Data ---------------------------------------------------------------
monocle.obj <- readRDS(file.path(dataDir, "/scRNAseq_2_monocle.cds"))
mon.old <- readRDS(file.path(dataDir, "old/scRNAseq_2_monocle_old.cds"))
new.clusters <- read.csv(file.path(tablesDir, "NewClusterML.csv"))
laia.clusters <- read.csv(file.path(tablesDir, "Ludwig4.csv"))

mon.old <- mon.old[,!grepl("fLN_40B3", colnames(mon.old))]
mon.old@colData <- mon.old@colData %>% droplevels()

# Assign New Clusters -----------------------------------------------------
new.clusters <- new.clusters %>% filter(NewClusterML != "")

new.clusters <- new.clusters %>% mutate(cluster = as.numeric(str_extract(new.clusters$NewClusterML, "\\d+")))

new.clusters$Barcode <- case_when(grepl("\\-1$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-1", "-1_fLN_40B2"),
                                  grepl("\\-2$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-2", "-1_fLN_40B3"),
                                  grepl("\\-3$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-3", "-1_fLN_41B1"),
                                  grepl("\\-4$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-4", "-1_fLN_41B2"),
                                  grepl("\\-5$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-5", "-1_fLN_41B3"),
                                  grepl("\\-6$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-6", "-1_fLN_B1"),
                                  grepl("\\-7$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-7", "-1_fSkin_40B2"),
                                  grepl("\\-8$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-8", "-1_fSkin_40B3"),
                                  grepl("\\-9$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-9", "-1_fSkin_41B1"),
                                  grepl("\\-10$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-10", "-1_fSkin_41B2"),
                                  grepl("\\-11$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-11", "-1_fSkin_41B3"),
                                  grepl("\\-12$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-12", "-1_fSkin_B1"),
                                  grepl("\\-13$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-13", "-1_mLN_40B2"),
                                  grepl("\\-14$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-14", "-1_mLN_40B3"),
                                  grepl("\\-15$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-15", "-1_mLN_41B1"),
                                  grepl("\\-16$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-16", "-1_mLN_41B3"),
                                  grepl("\\-17$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-17", "-1_mLN_41B4"),
                                  grepl("\\-18$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-18", "-1_mLN_B1"),
                                  grepl("\\-19$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-19", "-1_mSkin_40B2"),
                                  grepl("\\-20$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-20", "-1_mSkin_40B3"),
                                  grepl("\\-21$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-21", "-1_mSkin_41B1"),
                                  grepl("\\-22$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-22", "-1_mSkin_41B3"),
                                  grepl("\\-23$",new.clusters$Barcode) ~ str_replace(new.clusters$Barcode, "\\-23", "-1_mSkin_B1"))




clusters <- setNames(new.clusters$cluster, new.clusters$Barcode)
clusters <- clusters[colnames(monocle.obj)]
clusters <- as.factor(clusters)

clusters <- clusters %>% recode("57" = "34",
                                "42" = "14",
                                "51" = "39",
                                "9" = "7",
                                "62" = "47",
                                "40" = "12",
                                "43" = "19",
                                "27" = "26",
                                "18" = "13",
                                "30" = "2",
                                "3" = "1",
                                "5" = "1", 
                                "6" = "1",
                                "16" = "1",
                                "22" = "1"
                                )

#Just to check which clusters got removed and need replacement in order to have the clusters in a sequence
ol_c <- clusters[which(!(clusters %in% seq(length(levels(clusters)))))] %>% droplevels %>% levels
new_c <- seq(length(levels(clusters)))[which(!(seq(length(levels(clusters))) %in% clusters))]

ol_c
new_c

clusters <- clusters %>% recode_factor(
    "53" = "3",
    "54" = "5",
    "55" = "6",
    "56" = "9",
    "58" = "16",
    "59" = "18",
    "60" = "22",
    "61" = "27",
    "63" = "30",
    "64" = "40",
    "65" = "42",
    "66" = "43",
    "67" = "51"
)

#having levels in right order
for (i in as.character(length(levels(clusters)):1)){
    clusters <- relevel(clusters, i)
}




laia.clusters <- laia.clusters %>% filter(Ludwig3 != "")
laia.clusters <- laia.clusters %>% mutate(cluster = as.numeric(str_extract(laia.clusters$Ludwig3, "\\d+")))

laia.clusters$Barcode <- case_when(grepl("\\-1$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-1", "-1_fLN_40B2"),
                                   grepl("\\-2$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-2", "-1_fLN_40B3"),
                                   grepl("\\-3$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-3", "-1_fLN_41B1"),
                                   grepl("\\-4$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-4", "-1_fLN_41B2"),
                                   grepl("\\-5$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-5", "-1_fLN_41B3"),
                                   grepl("\\-6$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-6", "-1_fLN_B1"),
                                   grepl("\\-7$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-7", "-1_fSkin_40B2"),
                                   grepl("\\-8$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-8", "-1_fSkin_40B3"),
                                   grepl("\\-9$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-9", "-1_fSkin_41B1"),
                                   grepl("\\-10$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-10", "-1_fSkin_41B2"),
                                   grepl("\\-11$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-11", "-1_fSkin_41B3"),
                                   grepl("\\-12$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-12", "-1_fSkin_B1"),
                                   grepl("\\-13$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-13", "-1_mLN_40B2"),
                                   grepl("\\-14$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-14", "-1_mLN_40B3"),
                                   grepl("\\-15$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-15", "-1_mLN_41B1"),
                                   grepl("\\-16$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-16", "-1_mLN_41B3"),
                                   grepl("\\-17$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-17", "-1_mLN_41B4"),
                                   grepl("\\-18$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-18", "-1_mLN_B1"),
                                   grepl("\\-19$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-19", "-1_mSkin_40B2"),
                                   grepl("\\-20$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-20", "-1_mSkin_40B3"),
                                   grepl("\\-21$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-21", "-1_mSkin_41B1"),
                                   grepl("\\-22$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-22", "-1_mSkin_41B3"),
                                   grepl("\\-23$",laia.clusters$Barcode) ~ str_replace(laia.clusters$Barcode, "\\-23", "-1_mSkin_B1"))



# setdiff(
#     clusters %>% unname %>% unique %>% sort,
#     laia.clusters %>% filter(!is.na(cluster)) %>% pull(cluster) %>% unique %>% sort
# )

#laia.clusters[which(laia.clusters$Barcode %in% names((clusters[clusters == 3]))),] %>% pull(Ludwig3) %>% unique

laia.clusters$cluster <- case_when(laia.clusters$Ludwig3 == "GC B cells" ~ 50,
                                   laia.clusters$Ludwig3 == "Plasma cells" ~ 3,
                                   laia.clusters$Ludwig3 == "Pre-plasmablasts" ~ 53,
                                   TRUE ~ laia.clusters$cluster)


clusters <- setNames(laia.clusters$cluster, laia.clusters$Barcode)
clusters <- clusters[colnames(monocle.obj)]
clusters <- as.factor(clusters)

# Assigning to monocle.obj ------------------------------------------------
monocle.obj@clusters$UMAP$clusters <- clusters

monocle.obj@colData$Cluster <- unname(clusters(monocle.obj[,rownames(colData(monocle.obj))]))



monocle.obj@colData$celltype <- case_when(monocle.obj@colData$Cluster == 34 ~ "CD45.1+ T cells",
                                          monocle.obj@colData$Cluster == 44 ~ "Unassigned1",
                                          monocle.obj@colData$Cluster == 11 ~ "Unassigned2",
                                          monocle.obj@colData$Cluster == 5  ~ "Unassigned3",
                                          monocle.obj@colData$Cluster == 18 ~ "Myofibroblasts",
                                          monocle.obj@colData$Cluster == 50 ~ "GC B cells",
                                          monocle.obj@colData$Cluster == 3 ~ "Plasma cells",
                                          monocle.obj@colData$Cluster == 53 ~ "Pre-plasmablasts",
                                          TRUE ~ monocle.obj@colData$celltype)



monocle.obj@colData <- transform(
    monocle.obj@colData,
    ct_cluster = paste0(
        as.character(monocle.obj@colData$Cluster),
        " (",
        monocle.obj@colData$celltype
    ,
    ")"
    )
)

ordered_factor <-
    monocle.obj@colData %>% as_tibble() %>% arrange(celltype, Cluster) %>% select(celltype, Cluster)  %>%
    mutate(ct_cl = forcats::fct_inorder(paste0(as.character(Cluster),
                                               " (",
                                               celltype
                                               ,
                                               ")")))

monocle.obj@colData$ct_cluster <- factor(monocle.obj@colData$ct_cluster, levels = levels(ordered_factor$ct_cl))

# Save Object -------------------------------------------------------------

saveRDS(monocle.obj, file.path(dataDir, "scRNAseq_2a_monocle_with_fLN_40B3.cds"))

monocle.obj <- monocle.obj[,colnames(monocle.obj)[!grepl("fLN_40B3", colnames(monocle.obj))]]
monocle.obj@colData$sample <- monocle.obj@colData$sample %>% droplevels()

additional.meta <- readRDS(file.path(dataDir, "additional_meta.RDS"))
monocle.obj@colData <- cbind(monocle.obj@colData,additional.meta)

saveRDS(monocle.obj, file.path(dataDir, "scRNAseq_3_monocle_more_hash_cutoff.cds"))
