library(dplyr)
library(stringr)
library(monocle3)
library(Seurat)
library(data.table)
library(readr)

# Directories and Used Functions -----------------------------------------------------------
source('functions/directories_and_theme.R')
source('functions/load_raw_data.R')

#Unique Functions
#This is a huge function. The main things happening here are:
#(i) As the name suggest, assigning the meta data to the object.
#(ii) Integrating transcriptome data with CITE-seq data.
#(iii) Aggregating ab counts for treatment assignment
AssignMetadata <- function(sx, 
                           abcounts,
                           aggregated_treatment_threshold = 0.6,
                           treatment_threshold = 0.6,
                           cd4_threshold = 0.75,
                           cd45_threshold = 0.6,
                           gd_threshold = 0.75) {
    
    #Here we make sure we only use CITE-seq counts for barcodes withing
    #The data, as we removed a sample from our data. Additionally, we make 
    #sure that they are in the same order.
    barcodes <- row.names(sx@meta.data) %>% str_extract("^[:alpha:]+\\-\\d")
    abcounts <- abcounts[, colnames(abcounts) %in% barcodes]
    stopifnot(all(barcodes == colnames(abcounts)))
    
    #Here the aggregation for hashtags happens
    ab.hash <- abcounts[grepl("-[0-9]+$", row.names(abcounts)), ]
    ab.aggr <-
        list(
            "HTO-1" = 1:3,
            "HTO-2" = 4:6,
            "HTO-3" = 7:9
        ) %>%
        sapply(function(x)
            colSums(ab.hash[x, ])) %>%
        t() %>%
        as("sparseMatrix")
    
    stopifnot(colSums(ab.aggr) == colSums(ab.hash))
    
    #Seperating different CITE-seq counts
    ab.cd4 <- abcounts[grepl("CD[1-9]$", row.names(abcounts)), ]
    ab.cd45 <- abcounts[grepl("CD45.*$", row.names(abcounts)), ]
    
    #Setting to NULL in case the variable is assigned. Code might fail
    #otherwise for samples without gamma-delta CITE-seq data.
    ab.gd <- NULL
    ab.gdcd4 <- NULL
    ab.gdcd8 <- NULL
    #assigning gd in case they exist
    if (any(grepl("GD", rownames(abcounts)))) {
        ab.gd <- abcounts[grepl("GD", rownames(abcounts)), ]
        ab.gdcd4 <- abcounts[grepl("GD|CD4$", rownames(abcounts)), ]
        ab.gdcd8 <- abcounts[grepl("GD|CD8$", rownames(abcounts)), ]
    }
    
    #Function to assign CITE-seq to each cell
    assign_ab <- function(ab, threshold) {
        ratio <- t(t(as.matrix(ab)) / colSums(ab))
        #which ab fulfills the threshold criteria in each cell?
        apply(ratio, 2, function(col) {
            id <- which(col > threshold)
            if (length(id) == 0)
                return("Undefined")
            names(id)
        })
    }
    
    #The above function allows us to set thresholds for assigning a cell
    #You can see/set the in the function call. Standard values can be seen
    #on top.
    aggr <- assign_ab(ab.aggr, aggregated_treatment_threshold)
    hash <- assign_ab(ab.hash, treatment_threshold)
    cd4cd8 <- assign_ab(ab.cd4, cd4_threshold)
    cd45x <- assign_ab(ab.cd45, cd45_threshold)
    
    #same procedure for gamma delta, just with fail safe
    gdcd4 <- NULL
    gdcd8 <- NULL
    if (!is.null(ab.gd)) {
        gdcd4 <- assign_ab(ab.gdcd4, gd_threshold)
        gdcd8 <- assign_ab(ab.gdcd8, gd_threshold)
    }
    
    #Function to calculate the ratio of the hashtag with the highest bound
    #fraction (i.e. the assigned hashtag). 
    max_ratio <- function(x) {
        colMaxs(as.matrix(x)) / colSums(x)
    }
    
    #Here we simply create the metadata for the monocle object and assign the 
    #antibody data.
    sx@meta.data <- sx@meta.data %>%
        mutate(
            hash,
            aggr,
            cd4cd8,
            cd45x,
            gdcd4 = if (!is.null(ab.gd)) {
                gdcd4
            } else
                NA,
            gd.ratio = if (!is.null(ab.gd)) {
                ab.gdcd4["HTO-GD.TCR",] / colSums(ab.gdcd4)
            } else
                NA,
            gdcd8 = if (!is.null(ab.gd)) {
                gdcd8
            } else
                NA,
            gd.cd8ratio = if (!is.null(ab.gd)) {
                ab.gdcd4["HTO-GD.TCR",] / colSums(ab.gdcd8)
            } else
                NA,
            hash.sum = colSums(ab.hash),
            hash.ratio = max_ratio(ab.hash),
            hash.agg.ratio = max_ratio(ab.aggr),
            CD45 = ab.cd45["HTO-CD45.1",] / colSums(ab.cd45),
            CD45.sum = colSums(ab.cd45),
            CD4 = ab.cd4['HTO-CD4',] / colSums(ab.cd4),
            CD4.sum = colSums(ab.cd4)
        )
    #return the obejct at the end
    sx
}

# Create Data ---------------------------------------------------------------

#list all SoupX corrected files to load
fileDirs <- list.files(file.path(vDir, "countsSoupX"))

#load the data
data_list <- sapply(fileDirs, 
                    load_data,
                    path = file.path(vDir, "countsSoupX"),
                    USE.NAMES = FALSE)

#load CITE-seq data
CountsAB <- read_rds(file.path(dataDir, "0_ab_counts.rds"))

#LookupTable for empty droplets
LUT <- read_rds(file.path(dataDir, "0_doublet_LUT.rds"))
setDT(LUT)

# Assign Metadata ---------------------------------------------------------

#Assigning metadata with function defined above
data_list <- mapply(AssignMetadata, data_list, CountsAB)

#Determine %mito genes (important for QC/Identification of viable cells)
data_list <- lapply(data_list, function(x) {
    x$percent.mt <- PercentageFeatureSet(x, pattern = "^mt-")
    x
})

#cell cycle scoring. - self explanatory
data_list <- sapply(seq_along(data_list), function(x, y, i) {
    obj <- x[[i]]
    obj$sample <- y[[i]]
    obj <- NormalizeData(obj, verbose = T)
    obj <- CellCycleScoring(
        obj,
        s.features = str_to_title(cc.genes$s.genes),
        g2m.features = str_to_title(cc.genes$g2m.genes),
        set.ident = TRUE
    )
    setNames(list(obj), y[[i]])
}, x = data_list, y = names(data_list), USE.NAMES = F)


#filter dead cells; doublets;...
data_subset <-
    lapply(data_list,
           function(x) {
               obj <- x[,!colnames(x) %in% LUT[doublet_classification == "Doublet", rn]]
               subset(obj,
                      subset = nFeature_RNA > 200 &
                          percent.mt < 10 &
                          hash.sum > 10)
           })

#Create Monocle Obj
monocle.list <- lapply(data_subset, function(x){
    mat <- x@assays$RNA@counts
    cds <- new_cell_data_set(expression_data = mat,
                             cell_metadata = x@meta.data)
    cds
})

#Integrating all samples together in one large data set
monocle.obj <- combine_cds(cds_list = monocle.list, cell_names_unique = TRUE)
rowData(monocle.obj)$gene_short_name <- row.names(rowData(monocle.obj))

# Metadata refinement
#This part is a little bit messy. But the metadata gets refined for
#better later use (e.g. factors instead of character columns). Additional
#information such as organ and experiment group (HDAC1/HDAC2) gets annotated
#aswell.... and more!
monocle.obj@colData$hash <- as.factor(monocle.obj@colData$hash)
monocle.obj@colData$aggr <- as.factor(monocle.obj@colData$aggr)
monocle.obj@colData$gdcd4 <- as.factor(monocle.obj@colData$gdcd4)
monocle.obj@colData$gdcd8 <- as.factor(monocle.obj@colData$gdcd8)
monocle.obj@colData$cd4cd8 <- as.factor(monocle.obj@colData$cd4cd8)
monocle.obj@colData$cd45x <- as.factor(monocle.obj@colData$cd45x)
monocle.obj@colData$sample <- as.factor(monocle.obj@colData$sample)
monocle.obj@colData$organ <- as.factor(str_extract(monocle.obj@colData$sample, "[A-Z][:alpha:]+"))
levels1 <- monocle.obj@colData %>% 
    as_tibble %>% 
    filter(organ == "Skin") %>% 
    pull(sample) %>% 
    droplevels() %>% 
    levels
levels2 <- monocle.obj@colData %>% 
    as_tibble %>% 
    filter(organ == "LN") %>% 
    pull(sample) %>% 
    droplevels() %>% 
    levels
levelsx <- paste(c(levels1, levels2))
monocle.obj@colData$sample <- monocle.obj@colData$sample %>% forcats::fct_relevel(levelsx)
monocle.obj@colData$Phase <- as.factor(monocle.obj@colData$Phase)
monocle.obj@colData$batch <-
    as.factor(str_extract(monocle.obj@colData$sample, "[A-Z][0-9]$"))
monocle.obj@colData$treatment <-
    as.factor(case_when(
        grepl("30[1-3]$", monocle.obj@colData$hash) ~ "HDAC_WT",
        grepl("30[4-6]$", monocle.obj@colData$hash) ~ "HDAC_cKO",
        grepl("30[7-9]$", monocle.obj@colData$hash) ~ "NoT",
        TRUE ~ "Undefined"
    ))
monocle.obj@colData$treatment.agg <-
    as.factor(case_when(
        grepl("\\-3$", monocle.obj@colData$aggr) ~ "NoT",
        grepl("\\-1$", monocle.obj@colData$aggr) ~ "HDAC_WT",
        grepl("\\-2$", monocle.obj@colData$aggr) ~ "HDAC_cKO",
        TRUE ~ "Undefined"
    ))
monocle.obj@colData$hashid <-
    as.factor(ifelse(
        monocle.obj@colData$hash == "Undefined",
        "Undefined",
        str_extract(monocle.obj@colData$hash, "[0-9]+$")
    ))
monocle.obj@colData$hashid.agg <-
    as.factor(ifelse(
        monocle.obj@colData$aggr == "Undefined",
        "Undefined",
        str_extract(monocle.obj@colData$aggr, "[0-9]+$")
    ))
monocle.obj@colData$sex <- as.factor(str_to_upper(str_extract(monocle.obj@colData$sample,"^[A-z]")))
monocle.obj@colData$treatment <-
    monocle.obj@colData$treatment %>% forcats::fct_relevel("NoT", "HDAC_WT", "HDAC_cKO", "Undefined")
monocle.obj@colData$treatment.agg <-
    monocle.obj@colData$treatment.agg %>% forcats::fct_relevel("NoT", "HDAC_WT", "HDAC_cKO", "Undefined")
monocle.obj@colData$sample_treat <- as.factor(paste(monocle.obj@colData$sample, monocle.obj@colData$treatment, sep = "___"))
monocle.obj@colData$sample_treat.agg <-
    as.factor(paste(
        monocle.obj@colData$sample,
        monocle.obj@colData$treatment.agg,
        sep = "___"
    ))
monocle.obj@colData$experiment <-
    factor(case_when(
        !grepl("[0-9]+B[0-9]", monocle.obj@colData$sample) ~ "HDAC1",
        grepl("40", monocle.obj@colData$sample) ~ "HDAC1",
        TRUE ~ "HDAC2"
    ))


# QC Table ----------------------------------------------------------------

#coming soon...

# Export Data -------------------------------------------------------------

#save data for later use in different script
write_rds(monocle.obj, file.path(vDir, "1_qc.cds"))
