library(dplyr)
library(stringr)
library(monocle3)
library(Seurat)
library(data.table)
library(readr)

#setting up directories
vDir <- ("/vscratch/scRNAseq")
plotsDir <- file.path(vDir, "plots")
tablesDir <- file.path(vDir, "tables")
oldDir <- file.path(vDir, "data/old")
dataDir <-("data")
resDir <- ("results")

# Create Data ---------------------------------------------------------------

#directories to loop over
fileDirs <- list.files(file.path(vDir, "data/countsSoupX"))

load_data <- function(dat, 
                      path, 
                      slot = c("none", "gene expression", "antibody capture")){
    data <- Read10X(file.path(path, dat))
    
    slot <- match.arg(slot)
    if (slot == "gene expression")
        data <- data[[1]]
    else if (slot == "antibody capture")
        data <- data[[2]]
    
    obj <- CreateSeuratObject(
        counts = data,
        project = "scRNAseq",
        min.cells = 0,
        min.features = 0
    )
    x <- dat %>% str_extract("^[:alpha:]+_[:alnum:]+")
    
    setNames(list(obj), x)
}

data_list <- sapply(fileDirs, 
                    load_data,
                    path = file.path(vDir, "data/countsSoupX"),
                    USE.NAMES = FALSE)

CountsAB <- read_rds(file.path(dataDir, "ab_counts.rds"))


# Assign Metadata ---------------------------------------------------------

AssignMetadata <- function(sx, abcounts) {
    
    barcodes <- row.names(sx@meta.data) %>% str_extract("^[:alpha:]+\\-\\d")
    abcounts <- abcounts[, colnames(abcounts) %in% barcodes]
    stopifnot(all(barcodes == colnames(abcounts)))
    
    ab.hash <- abcounts[grepl("-[0-9]+$", row.names(abcounts)), ]
    ab.aggr <-
        list(
            "HTO-fLN-40B2-1" = 1:3,
            "HTO-fLN-40B2-2" = 4:6,
            "HTO-fLN-40B2-3" = 7:9
        ) %>%
        sapply(function(x)
            colSums(ab.hash[x, ])) %>%
        t() %>%
        as("sparseMatrix")
    
    stopifnot(colSums(ab.aggr) == colSums(ab.hash))
    
    ab.cd4 <- abcounts[grepl("CD[1-9]$", row.names(abcounts)), ]
    ab.cd45 <- abcounts[grepl("CD45.*$", row.names(abcounts)), ]
    
    if (any(grepl("GD", rownames(abcounts)))) {
        ab.gd <- abcounts[grepl("GD", rownames(abcounts)), ]
        ab.gdcd4 <- abcounts[grepl("GD|CD4$", rownames(abcounts)), ]
    } else {
        ab.gd <- NULL
        ab.gdcd4 <- NULL
    }
    
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
    
    aggr <- assign_ab(ab.aggr, 0.6)
    hash <- assign_ab(ab.hash, 0.6)
    cd4cd8 <- assign_ab(ab.cd4, 0.75)
    cd45x <- assign_ab(ab.cd45, 0.6)
    if (!is.null(ab.gd)) {
        gdcd4 <- assign_ab(ab.gdcd4, 0.75)
    } else {
        gdcd4 <- NULL
    }
    
    max_ratio <- function(x) {
        colMaxs(as.matrix(x)) / colSums(x)
    }
    
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
            hash.sum = colSums(ab.hash),
            hash.ratio = max_ratio(ab.hash),
            hash.agg.ratio = max_ratio(ab.aggr),
            CD45 = ab.cd45["HTO-CD45.1",] / colSums(ab.cd45),
            CD45.sum = colSums(ab.cd45),
            CD4 = ab.cd4['HTO-CD4',] / colSums(ab.cd4),
            CD4.sum = colSums(ab.cd4)
        )
    sx
}

#Assign meta
data_list <- mapply(AssignMetadata, data_list, CountsAB)

#Determine %mito genes
data_list <- lapply(data_list, function(x) {
    x$percent.mt <- PercentageFeatureSet(x, pattern = "^mt-")
    x
})

#cell cycle scoring
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

#LookupTable for empty droplets
LUT <- readr::read_rds(file.path(dataDir, "doublet_LUT.rds"))

setDT(LUT)
#filter dead cells; droplets;...
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

monocle.obj <- combine_cds(cds_list = monocle.list, cell_names_unique = TRUE)
rowData(monocle.obj)$gene_short_name <- row.names(rowData(monocle.obj))

# Metadata refinement
monocle.obj@colData$hash <- as.factor(monocle.obj@colData$hash)
monocle.obj@colData$aggr <- as.factor(monocle.obj@colData$aggr)
monocle.obj@colData$gdcd4 <- as.factor(monocle.obj@colData$gdcd4)
monocle.obj@colData$cd4cd8 <- as.factor(monocle.obj@colData$cd4cd8)
monocle.obj@colData$cd45x <- as.factor(monocle.obj@colData$cd45x)
monocle.obj@colData$sample <- as.factor(monocle.obj@colData$sample)
monocle.obj@colData$organ <- as.factor(str_extract(monocle.obj@colData$sample, "[A-Z][:alpha:]+"))
levels1 <- monocle.obj@colData %>% as_tibble %>% filter(organ == "Skin") %>% pull(sample) %>% droplevels() %>%  levels
levels2 <- monocle.obj@colData %>% as_tibble %>% filter(organ == "LN") %>% pull(sample) %>% droplevels() %>%  levels
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


#save data
write_rds(monocle.obj, file.path(dataDir, "1_qc.cds"))

