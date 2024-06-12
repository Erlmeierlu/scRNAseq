assign_ct_cluster <- function(cds, 
                              cluster_col = 'Cluster', 
                              celltype_col = 'celltype'){
    
    ct_cluster <- paste(cds[[cluster_col]], cds[[celltype_col]], sep = '_')
    
    #sorting by character AND numeric in same string 
    fct_order <-
        data.table(stringi::stri_split_fixed(unique(ct_cluster),
                                             "_",
                                             2,
                                             simplify = T))[
                                                 order(V2, as.numeric(V1))
                                             ][,
                                               do.call(paste,
                                                       c(.SD,
                                                         sep = '_'))]
    
    factor(ct_cluster, levels = fct_order)
}