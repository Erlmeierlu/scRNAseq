library(dplyr)
library(stringr)
library(monocle3)
library(Seurat)
library(data.table)
library(readr)
library(gt)

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

#same for raw data for the large QC plot at the end.
#for soupx the data got already filtered. Using them as the raw data
#reference would remove some information (it changes nothing in the
#analysis. Every otehr step is still performed with the SoupX corrected
#data. It's just for the QC table at the end)
raw_files <- paste0(
    grep(list.files(rawDir), pattern = "AGG", invert = T, value = T),
    "/filtered_feature_bc_matrix/"
)
raw_files <- grep("fLN_40B3", raw_files, value = T, invert = T)

raw_data <- sapply(raw_files, 
                   load_data,
                   path = rawDir, 
                   slot = "gene expression", 
                   USE.NAMES = FALSE)

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
#preprocessing for raw data
#Assigning metadata with function defined above
raw_data <- mapply(AssignMetadata, raw_data, CountsAB)

#Determine %mito genes (important for QC/Identification of viable cells)
raw_data <- lapply(raw_data, function(x) {
    x$percent.mt <- PercentageFeatureSet(x, pattern = "^mt-")
    x
})

#cell cycle scoring. - self explanatory
raw_data <- sapply(seq_along(raw_data), function(x, y, i) {
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
}, x = raw_data, y = names(raw_data), USE.NAMES = F)
#For later use, we extract the meta data of each sample
#and combine it into one large data table. This is the
#base of the large QC table produced at the end..
gt_data_raw <- lapply(raw_data,
                      function(x) {
                          getElement(x, name = 'meta.data') %>%
                              as.data.table(keep.rownames = 'barcode')
                      }) %>% rbindlist()

gt_data_raw[, experiment := ifelse(grepl('41', sample), 'HDAC2', 'HDAC1')]
gt_data_raw[, barcode := do.call(paste, args = list(barcode, sample, sep = '_'))]
gt_data_raw[, singlet := barcode %in% LUT[doublet_classification == 'Singlet', rn]]

gt_data_raw <- gt_data_raw[,
                           .(experiment,
                             sample,
                             nFeature_RNA,
                             percent.mt,
                             hash.sum,
                             hash.ratio,
                             hash.agg.ratio,
                             singlet)]

gt_data_raw <- gt_data_raw[,
                           by = sample,
                           .(
                               'Raw' = .N,
                               'Singlet' = sum(singlet),
                               'Unique Features' = sum(nFeature_RNA > 200),
                               'Mitochondrial Genes' = sum(percent.mt < 10),
                               'Hashtags Bound' = sum(hash.sum > 10),
                               'Hashtag Ratio' = sum(hash.ratio > 0.6, na.rm = T),
                               'Combined Hashtag Ratio' = sum(hash.agg.ratio > 0.6, na.rm = T),
                               'Filtered' = sum(nFeature_RNA > 200 &
                                                    percent.mt < 10 &
                                                    hash.sum > 10 &
                                                    singlet),
                               'Combined Hashtags' = sum(
                                   nFeature_RNA > 200 &
                                       percent.mt < 10 &
                                       hash.sum > 10 &
                                       singlet &
                                       hash.agg.ratio > 0.6
                               ),
                               'Individual Hashtags' = sum(
                                   nFeature_RNA > 200 &
                                       percent.mt < 10 &
                                       hash.sum > 10 &
                                       singlet &
                                       hash.ratio > 0.6
                               ),
                               Experiment = unique(experiment)
                           )]

gt_table <- gt(gt_data_raw, rowname_col = "sample", groupname_col = "Experiment")

gt_table <-
    gt_table %>%  tab_header(title = "Filtered Cells") %>%
    tab_spanner(
        label = "Cell Counts",
        columns = c(Raw, Filtered),
        gather = TRUE
    ) %>% tab_spanner(
        label = "Quality Control Criteria",
        columns = c(
            "Unique Features",
            "Mitochondrial Genes",
            "Hashtags Bound",
            'Singlet'
            ),
        gather = T) %>%
    tab_spanner(
        label = "Hashtag Ratios",
        columns = c(
            "Hashtag Ratio",
            "Combined Hashtag Ratio"),
        gather = T
    ) %>%
    tab_spanner(
        label = 'Assigned & Viable',
        columns = c(
            'Combined Hashtags', 
            'Individual Hashtags'
        )
    ) %>% 
    fmt_number(
        where(~ is.numeric(.x)),
        use_seps = T,
        decimals = 0) %>% 
    tab_style(style = list(cell_text(weight = "bold")),
              locations = cells_body(columns = Filtered)) %>%
    tab_style(style = list(cell_text(weight = "bold", style = "italic")),
              locations = cells_stub()) %>%
    tab_style(style = cell_fill(color = "gray90"),
              locations = list(cells_body(rows = contains(c(
                  "fSkin", "mSkin"
              ))),
              cells_stub(rows = contains(c(
                  "fSkin", "mSkin"
              ))))) %>%
    tab_footnote(footnote = md("Filtered with 'Quality Control Criteria'"),
                 location = cells_column_labels(columns = Filtered)) %>%
    tab_footnote(footnote = md("Filtered with 'Quality Control Criteria' & 'Hashtag Ratio'"),
                 location = cells_column_labels(columns = 'Individual Hashtags')) %>%
    tab_footnote(footnote = md("Filtered with 'Quality Control Criteria' & 'Combined Hashtag Ratio'"),
                 location = cells_column_labels(columns = 'Combined Hashtags')) %>%
    tab_footnote(footnote = md("Cells expressing more than **200 unique genes**"),
                 location = cells_column_labels(columns = 'Unique Features')) %>%
    tab_footnote(footnote = md("Cells with less than **10% mitochondrial reads**"),
                 location = cells_column_labels(columns = 'Mitochondrial Genes')) %>%
    tab_footnote(footnote = md("Cells with at least **10 hashtag reads**"),
                 location = cells_column_labels(columns = 'Hashtags Bound')) %>%
    tab_footnote(footnote = md("Cells with more than **60% reads** of **one hashtag**"),
                 location = cells_column_labels(columns = 'Hashtag Ratio')) %>%
    tab_footnote(footnote = md("Cells with more than **60% reads** of hashtags corresponding to a **treatment**"),
                 location = cells_column_labels(columns = 'Combined Hashtag Ratio')) %>% 
    tab_options(row_group.as_column = TRUE)

#To save the table we use the gtsave function. This creates a html file in this case.
#You can open the file and create screenshots, or save it as pdf. You will
#figure it out. It is also possible to save the table as png or pdf immediately.
#However, this is not working on the server and needs the 'webshot2' package which
#is not installed. 
gt_table %>% gtsave(file.path(plotsDir, "1_qc_table.html"))

# Export Data -------------------------------------------------------------

#save data for later use in different script
write_rds(monocle.obj, file.path(vDir, "1_qc.cds"))
