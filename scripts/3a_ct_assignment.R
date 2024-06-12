library(SingleR)
library(celldex)
library(tidyverse)
library(monocle3)
library(Seurat)
library(foreach)
library(data.table)
library(biomaRt)
library(scuttle)

# Directories and Used Functions -----------------------------------------------------------
source('functions/directories_and_theme.R')
source('functions/convert_genes.R')
source('functions/assign_ct_cluster.R')

#Unique Functions
#This function takes a subset of the monocle object as input. 
#It repeats the preprocessing steps in order to recluster the subsets.
recluster_subset <- function(cds, res = 0.0005) {
  cds <-
    preprocess_cds(cds, verbose = TRUE) %>%
    reduce_dimension(preprocess_method = "PCA", verbose = TRUE)
  
  # batch correction
  cds <- align_cds(cds,
                   alignment_group = "sample", 
                   residual_model_formula_str = "~Phase",
                   verbose = TRUE)
  
  cds <- reduce_dimension(cds,
                          reduction_method = "UMAP",
                          preprocess_method = "Aligned",
                          verbose = TRUE)
  
  cds <- cluster_cells(cds, resolution = res)
  cds
}

#This function relabels the clusters. Usually, in a monocle object the 
#clusters are ordered by cluster size. Meaning larger clusters have 
#smaller numbers. 
#In the end I want this to be true after reclustering the subsets and
#after reassigning the full dataset. 
relabel_clusters <- function(subset, cds){
  
  #This function increases the cluster levels by the 'increase' value
  #for example: c(1,2,3) with increase = 3 will result in c(4,5,6).
  increase_fac <- function(factor, increase){
    names <- names(factor)
    values <- unname(factor) %>% as.numeric
    values <- values + increase
    factor(setNames(values, names))
  }
  
  #This functions combines two sets of factors. Looking at the description
  #of ?forcats::fct_c() might make this function more clear:
  #"This is a useful way of patching together factors from multiple
  #sources that really should have the same levels but don't."
  comb_factors <- function(factor1, factor2){
    names <- c(names(factor1), names(factor2))
    fact <- forcats::fct_c(factor1, factor2) %>% droplevels
    
    names(fact) <- names
    fact
  }
  #This determines the value that we want to increase our factor levels by.
  #This is the number of clusters in our full dataset. First, we want the clusters
  #of the subset to start after the highest value in the full dataset.
  by_cl <- length(levels(cds@clusters@listData$UMAP$clusters))
  
  #First we increase the factor levels of the subset by the number of clusters
  #in our full dataset.
  subset@clusters@listData$UMAP$clusters <- increase_fac(clusters(subset), by_cl)
  
  #Then we combine the factors of the clusters of the full dataset with the
  #ones from the subset. To do so we seperate the cells that are in the subset
  #from those that are only in the full data set.
  cds@clusters$UMAP$clusters <- 
    comb_factors(cds@clusters$UMAP$clusters[!names(cds@clusters$UMAP$clusters) %in% names(clusters(subset))],
                 clusters(subset))
  
  #Make sure the values are in the right order.
  cds@clusters$UMAP$clusters <- 
    cds@clusters$UMAP$clusters[match(colnames(cds), names(cds@clusters$UMAP$clusters))]
  
  #Lastly, the factor levels are adjusted and ordered. First, the factors
  #are counted, and then ranked.
  tab <- fct_count(clusters(cds))
  r <- rank(-tab$n, ties.method = 'first')
  
  #Then we revalue the levels in the full dataset based on this ranking. 
  cds@clusters$UMAP$clusters <- clusters(cds) %>% 
    lvls_revalue(r %>% as.character) %>% 
    fct_inseq()
  
  #reset levels of subset clusters
  subset@clusters@listData$UMAP$clusters <- increase_fac(subset@clusters@listData$UMAP$clusters, 0)
  #assign Cluster column in metadata of the subset
  subset$Cluster <- clusters(subset) %>% unname()
  #in the end we return both the subset and the full dataset with new clusters
  #sorted by size.
  list(subset, cds)
}

#jsut a simple function to reshape singleR results for further use to wide 
#format.
cast_res <- function(dt) {
  dt[, dcast(.SD,
             cell ~ ref + label,
             value.var = 'labels')]
}
#function to set rownames. Used after merging meta data from an object
#with singleR results.
set_rownames <- function(cds){
  rownames(colData(cds)) <- colData(cds)$Row.names
  colData(cds)$Row.names <- NULL
  cds
}

#Based on the singleR results, we calculate the frequency of each
#annotated label within the clusters. The cells within each cluster
#are annotated seperately. This helps us decide visually how to assign the
#clusters at the end.
calc_freq <- function(results){
  #Calculate the frequency
  results <- results[,
                     keyby = .(ref, label, cluster, labels),
                     .N][,
                         by = .(ref, label, cluster),
                         freq := N / sum(N) * 100][freq > 5,-c('N')]
  #this just pastes two columns together into one. Important for later
  #use.
  results[,
          label := do.call(paste, c(.SD, sep = '_')),
          .SDcols = 1:2][]
}

#This function allows us to determine the majority of annotated celltype
#within each Cluster. This is done seperately for each reference database
#used. It is important that we supply the already annotated object, and the
#columns that should be used for the majority voting.
determine_majority <- function(cds, columns){
  #first we subset the meta data of the object to only have cluster and
  #annotation information
  out <- as.data.table(cds@colData)[, c('Cluster', ..columns)]
  #then we put it in long format 
  out <- out[, melt(.SD, 'Cluster', var = 'ref', val = 'label')]
  #counting of occurrence of each annotation label in each cluster/database 
  out <- out[,by = .(Cluster, ref, label), .N]
  #filter for the celltype that is mostly annotated in each cluster
  out <- out[out[, .I[which.max(N)], by = .(Cluster, ref)]$V1]
  
  #put it back into wide format (as we started with)
  out[, dcast(.SD, Cluster ~ ref, value.var = 'label')]
}

#This function makes the celltypes names more beautiful and 
#suitable for future analysis. Note that you might want to change
#what's inside this function. Here, I change two cell type names
#in a way that the main cell type (i.e. T cells) comes first, and 
#the subtype comes second: "NKT" becomes "T_NK". 
refine_strings <- function(results){
  #First we remove unwanted symbols from the names
  results[, celltype_raw := str_replace_all(celltype_raw, " |\\-", "_")]
  #Then we change specific names to fit the pattern. You might add more
  #here depending on the assigned names. Best is to follow the pattern
  #"MAINTYPE_SUBTYPE".
  results[, celltype_raw := fcase(grepl('NKT', celltype_raw), 'T_NK',
                                  grepl('Tgd', celltype_raw), 'T_GD',
                                  celltype_raw == 'CD40L_IL4_blasts', 'B_CD40L_IL4_blasts',
                                  rep(TRUE, .N), celltype_raw)]
  #Paste together cluster and celltype information to ct_cluster column
  results[,
          ct_cluster := do.call(paste, c(.SD, sep = '_')), 
          .SDcols = c('cluster', 'celltype_raw')]
  
  #assign factor levels
  fc_levels <- unique(results[order(celltype_raw, cluster), ct_cluster])
  results[,
          ct_cluster := factor(ct_cluster,
                               levels = fc_levels)]
  results[]
}

#This is commented out, because the plots generated with this function
#are redundent, and will be better displayed later on.
#The code that utilizes this function is also commented out in line 904++.

#remove comments if you want to use it.
# #Just a simple plotting function that will be used throughout the
# #script. It creates a heatmap indicating the frequency of assigned
# #celltypes withing each cluster.
# plot_ct_hm <- function(results){
#   results %>% 
#     ggplot(aes(x = labels, y = ct_cluster, fill = freq)) +
#     geom_tile(colour = "white", show.legend = FALSE) +
#     scale_fill_gradient(low = "ivory2", high = "red") +
#     theme_my() +
#     theme(
#       axis.title.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.text.x = element_text(
#         angle = 90,
#         size = 8,
#         hjust = 1,
#         vjust = 0.5
#       )
#     ) +
#     facet_grid(~ label,
#                scales = "free",
#                space = "free"
#     )
# }

# Load Data ---------------------------------------------------------------
monocle.obj <- read_rds(file.path(vDir, "2_unsupervised.cds"))

#The following reference object comes from this publication:
#https://doi.org/10.1016/j.jaci.2022.03.002
#Cite it!
reference.obj <- read_rds(file.path(vDir, '3a_annotation_reference_data.rds'))

# old.mon <- readRDS(file.path(oldDir, "scRNAseq_3_monocle_more_hash_cutoff.cds"))
# old_meta <- as.data.frame(colData(old.mon))
# old_counts <- assay(old.mon)
# old.mon <- CreateSeuratObject(counts = old_counts,
#                               project = "ol",
#                               assay = "counts",
#                               meta.data = old_meta)
# old.mon <-  as.SingleCellExperiment(old.mon)
# 
# colnames(colData(old.mon))[colnames(colData(old.mon)) == "celltype"] <- "label.fine"
# 
# old.mon@metadata <- list(type = 'sc')
# SingleR annotation ------------------------------------------------------

## Reference data ----------------------------------------------------------

### making a list of reference data sets  ----------------------------------------------

#Create a list of reference datasets. I commented most of them out, since
#They don't provide usable results in my opinion. If you want you can use
#all. This will take some time though.

ref_data <- list(
  #Human Data Sets: 
  ##HPCA generally / Keratinocytes
  # hpca = HumanPrimaryCellAtlasData(),
  # 
  ##DB ImmuneCells, comprehensive CD4+ subsets; only one B cell subset, no dendritic cells
  # dice = DatabaseImmuneCellExpressionData(),
  # 
  # ##Monaco Immunc cells
  # monaco = MonacoImmuneData(),
  
  #MouseData:
  ##MouseRNAseqData
  # mrsd = MouseRNAseqData(),
  
  ##ImmGen 
  immgen = ImmGenData(),
  
  ##reference Dataset
  struc = reference.obj
  
  ##old dataset
  # old = old.mon
)

#In case you used more datasets above you need to change gene names
#from human symbols to mouse symbols for the following. Not necessary if
#you choose the same reference datasets as me. 

#Human Gene names as Mouse Gene names
# rownames(ref_data$hpca) <- convert_gene_list(rownames(ref_data$hpca))
# rownames(ref_data$dice) <- convert_gene_list(rownames(ref_data$dice))
# rownames(ref_data$monaco) <- convert_gene_list(rownames(ref_data$monaco))

#Our Object where we wanna predict the celltypes.
sce <- monocle.obj
#needs logcounts assay. However, it is enough to use the normal counts for
#the analysis. Therefore, we assign our counts as logcounts.
sce@assays@data$logcounts <- sce@assays@data$counts

### run singleR with additional saving of original results data -------------------------------
# parallel computation
res <- list()
doParallel::registerDoParallel(cores = 7) #parallel package
#list where results will be stored in
sc_res <- foreach(ref = names(ref_data)) %dopar% { #for each reference database
  
  #e.g immgen has broad (label.main) and more detailed (label.fine) annotations
  #we want results for both first.
  for(labelx in c("label.main", "label.fine")){ 
    
    #This is done to not encounter an error if a dataset doesn't have a
    #label.main column. There might be better solutions but it works.
    if(!labelx %in% colnames(colData(ref_data[[ref]]))) next
    
    #For reference datasets that come from single-cell data, the wilcox
    #method in SingleR is more accurate. Datasets have a value for 
    #ref_data[[ref]]@metadata$type that is not NULL, and therefore the
    #wilcox method will be used. 
    use_method <- fifelse(is.null(ref_data[[ref]]@metadata$type),
                          'classic',
                          'wilcox')
    
    #If we have a single-cell reference dataset, this value will be TRUE,
    #and therefore, the datasets will be aggregated to pseudobulk samples.
    #This increases the speed by quite a bit.
    use_aggr <- fifelse(use_method == 'classic',
                        FALSE,
                        TRUE)
    
    #Here, the automated annotation is performed
    results <- SingleR(
      test = sce,
      ref = ref_data[[ref]],
      labels = colData(ref_data[[ref]])[, labelx],   
      aggr.ref = use_aggr,
      de.method = use_method
    )
    
    #store results in the list
    res[[labelx]] <- data.table(
      as_tibble(results, rownames = "cell"),
      ref = ref,
      labels = labelx
    )
    
    #make column names more beautiful
    colnames(res[[labelx]]) <- gsub("\\.", "_", colnames(res[[labelx]]))
    
    #extract relevent information
    res[[labelx]] <- res[[labelx]][, 
                                   .(cell, 
                                     labels, 
                                     tuning_scores_first, 
                                     tuning_scores_second)]
    #find numeric columns
    cols <- res[[labelx]][,
                          .SD,
                          .SDcols = is.numeric] %>%
      colnames()
    
    #round numeric columns with 2 digits to decrease size
    res[[labelx]][,
                  (cols) := round(.SD, 2),
                  .SDcols = cols]
    
  }
  #return list at the end
  res
}

sc_res <- setNames(sc_res, names(ref_data))
#combine the list into a dataframe
sc_res <- rbindlist(lapply(sc_res, 
                           rbindlist,
                           idcol = 'label'),
                    idcol = 'ref')

# Save Raw assignment data ------------------------------------------------
#We will use it in a later script aswell. Save it!
fwrite(sc_res, file.path(tablesDir, '3a_SingleR_res_full_dataset.csv'))

### assign labels to colData ----------------------------------------------------
#load data if needed
# sc_res <- fread(file.path(tablesDir, '3a_SingleR_res_full_dataset.csv'))

#Reshaping of the results in a longer format, so it fits the metadata
#format of the monocle object. We get a column for each reference database.
sc_res <- cast_res(sc_res)


#Adding singleR results to meta data. Here we merge two dataframes. 
#to make sure we match the right cells we use the row names of the metadata
#and the cell id of the results for merging.
colData(monocle.obj) <- merge(colData(monocle.obj), 
                              sc_res, 
                              by.x = 0, #using rownames
                              by.y = 'cell', #using cell column
                              sort = F) #dont sort output. 

#In the step above we lost the rownames. lets reassign them.
monocle.obj <- set_rownames(monocle.obj)

#For the first raw annotation I assign cells to KCs/FBs if they are annotated
#as such in the reference dataset from the publication above. The rest is 
#annotated as suggested by the immgen main reference database.
monocle.obj@colData$celltype_raw <- 
  fifelse(monocle.obj@colData$struc_label.fine %in% c('Keratinocytes', 'Fibroblasts'),
          monocle.obj@colData$struc_label.fine,
          monocle.obj@colData$immgen_label.main)


#Here we perform the majority voting for the assigned cell type per cluster.
ct_majority <- determine_majority(monocle.obj, 'celltype_raw')

#Lets assign the major celltype to each cluster.
#with the match we make sure we assign them in the right order.
monocle.obj@colData$celltype <- ct_majority[match(monocle.obj@colData$Cluster, 
                                                  ct_majority$Cluster),
                                            ]$celltype_raw

#Some renaming of the cells in order to preserve the pattern of e.g. T_cell,
#T_NK, etc.. Also, removing spaces and other symbols from the names..
monocle.obj@colData$celltype <-
  str_replace_all(monocle.obj@colData$celltype,
                  " |\\-",
                  "_") %>%
  case_when(. == 'NKT' ~ 'T_NK',
            . == 'Tgd' ~ 'T_GD',
            .default = .) %>% 
  factor() #Store it as a factor

#Let's take a look at the first annotation. We can already identify
#distinct clusters, such as b cells, t cells, and our structural cells.
plot_cells(monocle.obj, 
           color_cells_by = 'celltype', 
           label_groups_by_cluster = F,
           group_label_size = 3.5)
#We can save our first annotation.. UMAPs are best saved as JPGs, because
#PDFs would result in HUGE files due to every cell being plotted individually
ggsave(file.path(plotsDir, '3a_first_annotation_umap_full_dataset.jpg'))

#Create a column in the metadata that combines Cluster and celltype
#information. This function orders the factors by cluster and celltype name.
monocle.obj@colData$ct_cluster <- assign_ct_cluster(monocle.obj)

#If needed you can save the object now. However, it will not be used in a
#later script.
# write_rds(monocle.obj, file.path(vDir, "3a_pre_annotation_full_dataset.cds"))

# Subclustering -----------------------------------------------------------

# monocle.obj <- read_rds(file.path(vDir, "3a_pre_annotation_full_dataset.cds"))

#Explanation for barcodes that we load in the following:
#In the original analysis I used the function monocle3::choose_cells(monocle.obj).
#This function opens a window, that plots the UMAP of the object, and allows
#the user to select cells individually (or select areas). I looked at the
#assigned cell types in the full data set and then selected the most striking
#Clusters for each subset. I then exported the cell IDs of each subset in order
#to be able to reproduce these selections (without having to achieve the
# same selection again)

#T Subset
t_used <- read_rds(file.path(vDir, '3a_t_barcodes.rds'))
t_subset <- monocle.obj[, colnames(monocle.obj) %in% t_used]
t_subset@colData <- t_subset@colData %>% droplevels

t_subset <- recluster_subset(t_subset)

# write_rds(t_subset, file.path(vDir, '3a_pre_annotation_t_subset.cds'))

#B Subset
b_used <- read_rds(file.path(vDir, '3a_b_barcodes.rds'))
b_subset <- monocle.obj[, colnames(monocle.obj) %in% b_used]
b_subset@colData <- b_subset@colData %>% droplevels

b_subset <- recluster_subset(b_subset)

# write_rds(b_subset, file.path(vDir, '3a_pre_annotation_b_subset.cds'))

#M Subset
m_used <- read_rds(file.path(vDir, '3a_m_barcodes.rds'))
m_subset <- monocle.obj[, colnames(monocle.obj) %in% m_used]
m_subset@colData <- m_subset@colData %>% droplevels

m_subset <- recluster_subset(m_subset)

# write_rds(m_subset, file.path(vDir, '3a_pre_annotation_m_subset.cds'))

# Integrating New Clusters back into data ---------------------------------

#Most of the following is done to order the factor levels in the final monocle
#object by factor size 
# monocle.obj <- read_rds(file.path(vDir, "3a_pre_annotation_full_dataset.cds"))

#The object of the full dataset gets progressively changed. Please make
#sure to follow every step one after the other. 
## T_subset ----------------------------------------------------------------

# t_subset <- read_rds(file.path(vDir, '3a_pre_annotation_t_subset.cds'))
#these functions are explained on top of the script. 
relabel_res <- relabel_clusters(t_subset, monocle.obj)
t_subset <- relabel_res[[1]]
monocle.obj <- relabel_res[[2]]

## M_subset ----------------------------------------------------------------
# m_subset <- read_rds(file.path(vDir, '3a_pre_annotation_m_subset.cds'))
relabel_res <- relabel_clusters(m_subset, monocle.obj)
m_subset <- relabel_res[[1]]
monocle.obj <- relabel_res[[2]]

## B_subset ----------------------------------------------------------------

# b_subset <- read_rds(file.path(vDir, '3a_pre_annotation_b_subset.cds'))
relabel_res <- relabel_clusters(b_subset, monocle.obj)
b_subset <- relabel_res[[1]]
monocle.obj <- relabel_res[[2]]

# Reassigning colData -----------------------------------------------------
monocle.obj$Cluster <- clusters(monocle.obj) %>% unname()
monocle.obj$ct_cluster <- assign_ct_cluster(monocle.obj)

# Ref Data for T subsets --------------------------------------------------

### TH Raw data
#This reference dataset for t cells comes from the following publication
#https://doi.org/10.1186/s13062-015-0045-x
#Please Cite this paper, I guess.

#download data
th_counts <- fread('https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/582/E-MTAB-2582/Files/Teichmann-ThExpress_rawCounts.txt')

#Use the right BioMart dataset with the correct gene IDs that are used in
#our dataset. 
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                   dataset="mmusculus_gene_ensembl", 
                   host = 'http://oct2014.archive.ensembl.org')

#here we extract the gene IDs for each ensemble ID in the ref data
genes <- getBM(filters = 'ensembl_gene_id', 
               attributes= c("ensembl_gene_id",
                             'external_gene_name'),
               values = th_counts$ENSEMBL_ID,
               mart = ensembl)

#then we assign the IDs to the data
th_counts <- th_counts[ENSEMBL_ID %in% genes$ensembl_gene_id, #for each row which ENSEMBL_ID is in the ensemble database...
                       gene_id := genes$external_gene_name #...assign the corresponding gene ID
                        ][, 
                          na.omit(.SD) #then, remove rows that were not in the ensemble database
                        ][,
                          lapply(.SD, mean), #Some gene IDs are duplicated
                          by = gene_id,      #I simply take the average
                          .SDcols = -1.      #of the counts and combine them.
                        ][gene_id != '']    #at the end I remove genes with an empty name

#Set the counts as a dataframe and assign gene IDs as rownames
th_counts <- as.data.frame(th_counts)
rownames(th_counts) <- th_counts$gene_id
th_counts$gene_id <- NULL

#Convert to matrix first, and then to a sparseMatrix
th_counts <- as(as.matrix(th_counts), 'sparseMatrix')

#Now we can load the meta data
th_meta <- fread('https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/582/E-MTAB-2582/Files/E-MTAB-2582.sdrf.txt')

#getting the metadata in a shape that allows us to combine it to a 
#singlecellexperiment object
th_meta <- th_meta[, .(`Extract Name`, label.fine = str_extract(`Extract Name`, '^[:alnum:]+'))] %>% unique()
th_meta <- data.frame(label.fine = th_meta$label.fine, row.names = colnames(th_counts))

th_raw <- SingleCellExperiment(list(counts = th_counts), colData = DataFrame(th_meta))
th_raw <- logNormCounts(th_raw)

### Rest of T data
#This reference dataset comes from following publication:
#https://doi.org/10.1016/j.immuni.2019.05.014
#Maybe also cite.

#Chris Tibbit sent me this file personally. Please find it in the tables
#directory
t_sort_info <- readxl::read_xlsx(file.path(tablesDir, 'Sort format for Single cells HDM.xlsx'))
#remove unwanted columns and rows
t_sort_info <- setDT(t_sort_info[-1, -2])
#pivot data into long format so we end up with labels assigned to each sample
t_sort_info <- t_sort_info[, melt(.SD, #long format
                                  id.vars = 1,
                                  value.name = 'label.fine')][, 
                                                              sample := do.call(paste, c(.SD, sep = '')), #make sample column
                                                              .SDcols = 1:2][, -(1:2)]
#if we have NA as label change it to unknown
t_sort_info[is.na(label.fine), label.fine := 'unknown'] 

# set as data frame and assign sample column as row name
t_sort_info <- as.data.frame(t_sort_info)
rownames(t_sort_info) <- t_sort_info$sample
t_sort_info$sample <- NULL

#The counts can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131935
#this is done automatically here
t_counts <- fread('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131935&format=file&file=GSE131935%5FSS2%5F17%5F449%5Frpkms%2Etab%2Egz')
#again take the mean expression for duplicate genes
t_counts <- t_counts[, lapply(.SD, mean), by = gene, .SDcols = -1]
#set rownames..
t_counts <- as.data.frame(t_counts)
rownames(t_counts) <- t_counts$gene
t_counts$gene <- NULL
#convert to sparse Matrix
t_counts <- as(as.matrix(t_counts), 'sparseMatrix')

#Create a SCE object out of counts and meta data.
th_rest <- SingleCellExperiment(list(counts = t_counts),
                                 colData = DataFrame(t_sort_info))
th_rest <- logNormCounts(th_rest)

#we put the reference datasets into a list again
ref_data <- list(
  th_raw = th_raw,
  th_rest = th_rest
)

sce <- monocle.obj
#needs logcounts assay. However, it is enough to use the normal counts for
#the analysis. Therefore, we assign our counts as logcounts.
sce@assays@data$logcounts <- sce@assays@data$counts

#basically same approach for automated annotation than before
#explanations are not repeated here.
res <- list()
doParallel::registerDoParallel(cores = 7)
sc_res <- foreach(ref = names(ref_data)) %dopar% {
  for(labelx in c("label.main", "label.fine")){
    
    if(!labelx %in% colnames(colData(ref_data[[ref]]))) next
    
    use_method <- fifelse(is.null(ref_data[[ref]]@metadata$type),
                          'classic',
                          'wilcox')
    
    use_aggr <- fifelse(use_method == 'classic',
                        FALSE,
                        TRUE)
    
    results <- SingleR(
      test = sce,
      ref = ref_data[[ref]],
      labels = colData(ref_data[[ref]])[, labelx],   
      aggr.ref = use_aggr,
      de.method = use_method
    )
    
    res[[labelx]] <- data.table(
      as_tibble(results, rownames = "cell"),
      ref = ref,
      labels = labelx
    )
    
    colnames(res[[labelx]]) <- gsub("\\.", "_", colnames(res[[labelx]]))
    
    res[[labelx]] <- res[[labelx]][, 
                                   .(cell, 
                                     labels, 
                                     tuning_scores_first, 
                                     tuning_scores_second)]
    
    cols <- res[[labelx]][,
                          .SD,
                          .SDcols = is.numeric] %>%
      colnames()
    
    res[[labelx]][,
                  (cols) := round(.SD, 2),
                  .SDcols = cols]
    
  }
  res
}

sc_res <- setNames(sc_res, names(ref_data))

sc_res <- rbindlist(lapply(sc_res, rbindlist,
                           idcol = 'label'),
                    idcol = 'ref')

## Save Raw assignment data ------------------------------------------------
#we will use it later! save
fwrite(sc_res, file.path(tablesDir, '3a_SingleR_res_t_subset.csv'))

# B Subsets Ref Data ------------------------------------------------------

## Bregs --------------------------------------------------------------------

#I didn't use this data because it was not really usable, because I 
#could not really determine the cell type names..
#Maybe you wanna try it again and figure out what the cell types are

# #Data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8476928/
# #Cite?
# 
# #count data
# url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE174739&format=file&file=GSE174739%5Ffiltered%5F10X%5FB%5Fcounts%2Etxt%2Egz'
# breg_counts <- as.matrix(fread(url), rownames = 1)
# breg_counts <- as(breg_counts, 'sparseMatrix')
# 
# #meta data
# url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE174739&format=file&file=GSE174739%5Ffiltered%5F10X%5FB%5Fmetadata%2Etxt%2Egz'
# breg_meta <-  as.data.frame(fread(url))
# rownames(breg_meta) <- breg_meta$V1
# breg_meta$V1 <- NULL
# 
# breg <- SingleCellExperiment(list(counts = breg_counts), 
#                              colData = DataFrame(breg_meta))
# breg <- scuttle::logNormCounts(breg)
# breg@metadata <- list(type = 'sc')

# breg$label.fine <- ???????? #assign cell type here

## B rest---------------------------------------------------------------------

#Data from https://www.nature.com/articles/ni.3154
url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60927&format=file&file=GSE60927%5FRaw%5Fcounts%2Etxt%2Egz'
b <- fread(url)
#only use genes that are also in our subset
b <- b[Symbol %in% rownames(b_subset)]
#note some symbols are duplicated in data
#find out which Symbols are 
dupe_id <- b$Symbol %>% anyDuplicated()
#exclude from data
b <- b[!dupe_id]

#assign symbol as rowname
b <- as.data.frame(b)
row.names(b) <- b$Symbol
b$EntrezID <- b$Symbol <- b$GeneLength <- NULL
#convert to sparse matrix
b <- as(as.matrix(b), 'sparseMatrix')

#create meta data
b_meta <- data.frame(label.fine = str_remove(colnames(b), '_rep.?'),
                     row.names = colnames(b))

#fix cell type names to fit the pattern of "MAINTYPE_SUBTYPE"
b_meta$label.fine <- case_when(b_meta$label.fine == 'FoB' ~ 'B_Follicular',
                               b_meta$label.fine == 'MZB' ~ 'B_Marginal-zone',
                               b_meta$label.fine == 'GCB' ~ 'B_Germinal-center',
                               b_meta$label.fine == 'SplPB' ~ 'Plasmablasts-splenic',
                               b_meta$label.fine == 'SplPC' ~ 'Plasmacells-splenic',
                               b_meta$label.fine == 'BMPC' ~ 'Plasmacells-bonemarrow',
                               b_meta$label.fine == 'A20' ~ 'B_Germinal-center-A20',
                               b_meta$label.fine == 'MPC11' ~ 'Plasmacytoma',
                               .default = b_meta$label.fine)

b_rest <- SingleCellExperiment(list(counts = b),
                               colData = b_meta)
b_rest <- scuttle::logNormCounts(b_rest)

## Run singleR -------------------------------------------------------------

sce <- monocle.obj
#needs logcounts assay. However, it is enough to use the normal counts for
#the analysis. Therefore, we assign our counts as logcounts.
sce@assays@data$logcounts <- sce@assays@data$counts

ref_data <- list(b_data = b_rest
                 # ,
                 # b_reg = breg
                 )

#same approach as before..
res <- list()
doParallel::registerDoParallel(cores = 7)
b_res <- foreach(ref = names(ref_data)) %dopar% {
  for(labelx in c("label.main", "label.fine")){
    
    if(!labelx %in% colnames(colData(ref_data[[ref]]))) next
    
    use_method <- fifelse(is.null(ref_data[[ref]]@metadata$type),
                          'classic',
                          'wilcox')
    
    use_aggr <- fifelse(use_method == 'classic',
                        FALSE,
                        TRUE)
    
    results <- SingleR(
      test = sce,
      ref = ref_data[[ref]],
      labels = colData(ref_data[[ref]])[, labelx],   
      aggr.ref = use_aggr,
      de.method = use_method
    )
    
    res[[labelx]] <- data.table(
      as_tibble(results, rownames = "cell"),
      ref = ref,
      labels = labelx
    )
    
    colnames(res[[labelx]]) <- gsub("\\.", "_", colnames(res[[labelx]]))
    
    res[[labelx]] <- res[[labelx]][, 
                                   .(cell, 
                                     labels, 
                                     tuning_scores_first, 
                                     tuning_scores_second)]
    
    cols <- res[[labelx]][,
                          .SD,
                          .SDcols = is.numeric] %>%
      colnames()
    
    res[[labelx]][,
                  (cols) := round(.SD, 2),
                  .SDcols = cols]
    
  }
  res
}

b_res <- setNames(b_res, names(ref_data))

b_res <- rbindlist(lapply(b_res, rbindlist,
                           idcol = 'label'),
                    idcol = 'ref')

fwrite(b_res, file.path(tablesDir, '3a_SingleR_res_b_subset.csv'))


# Assign colData Labels ---------------------------------------------------

sc_res <- fread(file.path(tablesDir, '3a_SingleR_res_t_subset.csv'))
#for large set
full_res <- fread(file.path(tablesDir, '3a_SingleR_res_full_dataset.csv'))
#same for t subset
t_res <- sc_res[cell %chin% colnames(t_subset)]
#for m subset
m_res <- full_res[cell %chin% colnames(m_subset)]
#load b res data
full_b_res <- fread(file.path(tablesDir, '3a_SingleR_res_b_subset.csv'))
#only for b cells
b_res <- full_b_res[cell %chin% colnames(b_subset)]

#apply same data transformations to result tables
sc_res <- cast_res(sc_res)
full_res <- cast_res(full_res)
m_res <- cast_res(m_res)
t_res <- cast_res(t_res)
full_b_res <- cast_res(full_b_res)
b_res <- cast_res(b_res)

#merge the datasets with the singleR results
colData(monocle.obj) <- merge(colData(monocle.obj), 
                              sc_res, 
                              by.x = 0, 
                              by.y = 'cell',
                              sort = F)

colData(monocle.obj) <- merge(colData(monocle.obj), 
                              full_b_res, 
                              by.x = 'Row.names', 
                              by.y = 'cell',
                              sort = F)

colData(t_subset) <- merge(colData(t_subset), 
                              t_res, 
                              by.x = 0, 
                              by.y = 'cell',
                              sort = F)

colData(b_subset) <- merge(colData(b_subset), 
                           b_res, 
                           by.x = 0, 
                           by.y = 'cell',
                           sort = F)

monocle.obj <- set_rownames(monocle.obj)
t_subset <- set_rownames(t_subset)
b_subset <- set_rownames(b_subset)

res <- fread(file.path(tablesDir, '3a_SingleR_res_t_subset.csv'))
res2 <- fread(file.path(tablesDir, '3a_SingleR_res_full_dataset.csv'))
res3 <- fread(file.path(tablesDir, '3a_SingleR_res_b_subset.csv'))
res <- rbind(res,res2, res3)

res_t <- res[cell %chin% colnames(t_subset)]
res_b <- res[cell %chin% colnames(b_subset)]
res_m <- res[cell %chin% colnames(m_subset)]

res[, cluster := clusters(monocle.obj)[res$cell]]
res_t[, cluster := clusters(t_subset)[res_t$cell]]
res_b[, cluster := clusters(b_subset)[res_b$cell]]
res_m[, cluster := clusters(m_subset)[res_m$cell]]

res <- calc_freq(res)
res_t <- calc_freq(res_t)
res_b <- calc_freq(res_b)
res_m <- calc_freq(res_m)

#Maybe newer version of ct_majority with majority of all ref data bases
cols <- colnames(monocle.obj@colData) %>% grep(pattern = 'label', value = T)
cols_t <- colnames(t_subset@colData) %>% grep(pattern = 'label', value = T)
cols_b <- colnames(b_subset@colData) %>% grep(pattern = 'label', value = T)
cols_m <- colnames(m_subset@colData) %>% grep(pattern = 'label', value = T)

ct_majority <- determine_majority(monocle.obj, cols)
#combination of old data and one t subset ref data
ct_majority[, th_raw_label.fine := fcase(th_raw_label.fine == 'Naive', 'T_naive',
                                         th_raw_label.fine == 'iTreg', 'Treg_induced',
                                         rep(TRUE, .N), th_raw_label.fine)]

ct_majority[, celltype_raw := fcase(struc_label.fine == 'T cells', 
                                    th_raw_label.fine,
                                    struc_label.fine %in% c('Keratinocytes', 'Fibroblasts'),
                                    struc_label.fine,
                                    rep(TRUE, .N), immgen_label.main)]
ct_majority <- ct_majority[, .(Cluster, celltype_raw)]

#for t subset
ct_majority_t <- determine_majority(t_subset, cols_t)
#combination of old data and one t subset ref data
ct_majority_t[, th_raw_label.fine := fcase(th_raw_label.fine == 'Naive', 'T_naive',
                                         th_raw_label.fine == 'iTreg', 'Treg_induced',
                                         rep(TRUE, .N), th_raw_label.fine)]
ct_majority_t[, celltype_raw := th_raw_label.fine]
ct_majority_t <- ct_majority_t[, .(Cluster, celltype_raw)]

#for b subset
ct_majority_b <- determine_majority(b_subset,cols_b)
ct_majority_b <- ct_majority_b[, celltype_raw := b_data_label.fine]
ct_majority_b <- ct_majority_b[, .(Cluster, celltype_raw)]

#for m subset
ct_majority_m <- determine_majority(m_subset, cols_m)
ct_majority_m <- ct_majority_m[, celltype_raw := immgen_label.main]
ct_majority_m <- ct_majority_m[, .(Cluster, celltype_raw)]

res_t <- merge(res_t, ct_majority_t, by.x = 'cluster', by.y = 'Cluster')
res_b <- merge(res_b, ct_majority_b, by.x = 'cluster', by.y = 'Cluster')
res <- merge(res, ct_majority, by.x = 'cluster', by.y = 'Cluster')
res_m <- merge(res_m, ct_majority_m, by.x = 'cluster', by.y = 'Cluster')

res <- refine_strings(res)
res_t <- refine_strings(res_t)
res_b <- refine_strings(res_b)
res_m <- refine_strings(res_m)

# # plot a heatmap
# res %>% 
#   plot_ct_hm()
# ggsave(file.path(plotsDir, "all_clusters_all_ct_labels_HM.pdf"), height = 15, width = 28)
# 
# # plot a heatmap
# res_t %>% 
#   plot_ct_hm()
# ggsave(file.path(plotsDir, "t_all_clusters_all_ct_labels_HM.pdf"), height = 15, width = 28)
# 
# # plot a heatmap
# res_b %>% 
#   plot_ct_hm()
# ggsave(file.path(plotsDir, "b_all_clusters_all_ct_labels_HM.pdf"), height = 15, width = 28)
# 
# # plot a heatmap
# res_m %>% 
#   plot_ct_hm()
# ggsave(file.path(plotsDir, "m_all_clusters_all_ct_labels_HM.pdf"), height = 9, width = 28)

# First assignment of broad cell types ------------------------------------------------

#Full object
monocle.obj$celltype <- res$celltype_raw[match(monocle.obj$Cluster, res$cluster)]
monocle.obj$ct_cluster <- assign_ct_cluster(monocle.obj)

#T Subset
t_subset$celltype <- res_t$celltype_raw[match(t_subset$Cluster, res_t$cluster)]
t_subset$celltype <- case_when(t_subset$cd45x == 'HTO-CD45.1' ~ 
                                 'T_CD45.1',
                               t_subset$gdcd4 == 'HTO-GD.TCR' ~ 
                                 'T_GD',
                               grepl('^T', t_subset$celltype) ~ 
                                 'T_cells',
                               TRUE ~ 
                                 t_subset$celltype)

t_subset$ct_cluster <- assign_ct_cluster(t_subset)

#B Subset
b_subset$celltype <- res_b$celltype_raw[match(b_subset$Cluster, res_b$cluster)]

b_subset$ct_cluster <- assign_ct_cluster(b_subset)

#M Subset
m_subset$celltype <- res_m$celltype_raw[match(m_subset$Cluster, res_m$cluster)]

m_subset$ct_cluster <- assign_ct_cluster(m_subset)

# Save Objects ------------------------------------------------------------
write_rds(monocle.obj, file.path(vDir, '3a_first_annotation_full_dataset.cds'))
write_rds(b_subset, file.path(vDir, '3a_first_annotation_b_subset.cds'))
write_rds(t_subset, file.path(vDir, '3a_first_annotation_t_subset.cds'))
write_rds(m_subset, file.path(vDir, '3a_first_annotation_m_subset.cds'))
