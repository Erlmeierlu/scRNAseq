load_celltypist <- function(filename){
    dt <- fread(filename)
    dt[, let(
        cell = V1,
        ref = 'CellTypist',
        label = 'label.fine',
        labels = predicted_labels
    )]
    
    dt[,
       .(cell,
         ref,
         label,
         labels)]
}

load_singler <- function(filename, subset){
    
    if(length(filename > 1)){
        dt <- lapply(filename, fread) %>% rbindlist()
    }else{dt <- fread(filename)}
    
    dt <- dt[cell %in% colnames(subset)]
    dt[,
       .(cell,
         ref,
         label,
         labels)]
}

prepare_res <- function(ds = NULL, dc = NULL, subset){
    dt <- rbind(ds, dc)
    
    dt[, let(cluster = clusters(subset)[dt$cell],
             celltype = subset$celltype[match(dt$cell, colnames(subset))])]
    
    dt <- dt[,
             keyby = .(ref, label, cluster, celltype, labels),
             .N][,
                 by = .(ref, label, cluster, celltype),
                 freq := N / sum(N) * 100][freq > 5, -c('N')]
    
    dt[,
       label := do.call(paste, c(.SD, sep = '_')),
       .SDcols = 1:2]
    
    dt[,
       ct_cluster := do.call(paste, c(.SD, sep = '_')), 
       .SDcols = c('cluster', 'celltype')][]
}

assign_factor_levels <- function(dt, order_by){
    ord <- dt[label == order_by,
              sum(freq),
              by = labels][order(V1, decreasing = TRUE), labels] 
    
    mat <- dt[,
              dcast(.SD,
                    label + cluster + celltype + ct_cluster ~ labels,
                    value.var = 'freq',
                    fill = 0,
                    subset = .(label == order_by))
    ]
    setcolorder(mat, c('label', 'cluster', 'celltype', 'ct_cluster', ord))
    
    #ct_cluster levels
    dist <- dist(mat[,-(1:4)])
    cl <- hclust(dist)
    fc_levels <- fct_inorder(mat$ct_cluster[cl$order])
    
    #assignment levels
    rest <- dt[!labels %in% ord][order(labels), labels] %>% unique()
    col_dist <- dist(t(mat[,-(1:4)]))
    cl_col <- hclust(col_dist)
    x_levels <- fct_inorder(c(cl_col$labels[cl_col$order], rest))
    
    dt <- copy(dt)
    dt[, ct_cluster := factor(ct_cluster,
                              levels = fc_levels)]
    dt[, labels := factor(labels, levels = x_levels)][]
}

plot_automated_anno <- function(dt, 
                                group = c('ct_cluster', 'celltype'),
                                use_ref = NULL, 
                                ignore_ref = NULL){
    
    group <- match.arg(group, c('ct_cluster', 'celltype'))
    
    if(!is.null(use_ref) & !is.null(ignore_ref)){
        ignore_ref <- NULL
        warning('"use_ref" and "ignore_ref" were provided. Setting "ignore_ref" to NULL.')
    }
    if(!is.null(use_ref)) dt <- dt[label %in% use_ref]
    if(!is.null(ignore_ref)) dt <- dt[label %notin% ignore_ref]
    
    ggplot(dt,
           aes(x = labels, y = !!ensym(group), fill = freq)) +
        geom_tile(colour = "white", show.legend = FALSE) +
        scale_fill_gradient(low = "ivory2", high = "red") +
        theme_my() +
        theme(
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(
                angle = 90,
                size = 8,
                hjust = 1,
                vjust = 0.5
            )
        ) +
        facet_grid( ~ label,
                    scales = "free",
                    space = "free") + 
        ggtitle('Automated cell type annotation')
}

generate_abdata <- function(subset, 
                            levels, 
                            abgroups = c('hashid', 'cd4cd8', 'cd45x'),
                            group = c('ct_cluster', 'celltype')){
    
    group = match.arg(group, 
                      c('ct_cluster', 'celltype'))
    
    abgroups <- match.arg(abgroups, 
                          c('hashid', 'cd4cd8', 'cd45x'), 
                          several.ok = TRUE)
    
    dt <- subset@colData %>% as.data.table()
    
    dt[,
       melt(.SD, 
            measure = abgroups,
            var = 'abgroup',
            val = 'abtype')
    ][,
      keyby = .(col, abgroup, abtype),
      .N,
      env = list(col = group)
    ][,
      keyby = .(col, abgroup),
      .(fraction = N/sum(N),
        abtype = abtype %>% 
            str_replace('Undefined', 'Equally Distributed') %>% 
            factor()),
      env = list(col = group)
    ][,
      let(group = factor(group, levels = levels),
          abtype = abtype %>% 
              fct_relevel('Equally Distributed', after = Inf)),
      env = list(col = group)
    ][]
}

plot_abs <- function(dt, group = c('ct_cluster', 'celltype')){
    group <- match.arg(group, c('ct_cluster', 'celltype'))
    ggplot(dt) +
        geom_tile(aes(abtype, !!ensym(group), fill = fraction),
                  col = "white", 
                  width = 0.75,
                  show.legend = F) +
        scale_fill_gradient(low = "ivory2", high = "red") +
        facet_wrap(~abgroup,
                   scales = "free"
        ) +
        theme_my() +
        theme(
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(
                angle = 90,
                size = 8,
                hjust = 1,
                vjust = 0.5,
            ),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
        ) + 
        ggtitle('CITE-seq Antibodies & Hashtags')
}
plot_genes <- function (cds, 
                        markers = list(), 
                        group_cells_by = c('ct_cluster', 'celltype'), 
                        lower_threshold = 0, 
                        max.size = 3, 
                        ordering_type = c("cluster_row_col",
                                          'maximal_on_diag',
                                          "none"),
                        cell_factor_levels,
                        pseudocount = 1,
                        scale_max = 3, 
                        scale_min = -3) {
    
    group_cells_by = match.arg(group_cells_by, c('ct_cluster', 'celltype'))
    
    ordering_type <- match.arg(ordering_type,
                               c("cluster_row_col",
                                 'maximal_on_diag',
                                 "none"))
    
    markers <- lapply(markers, unique)
    
    gene_ids <- lapply(markers, function(x) {
        as.data.frame(fData(cds)) %>% tibble::rownames_to_column() %>%
            dplyr::filter(rowname %in% x | gene_short_name %in%
                              x) %>% dplyr::pull(rowname)
    })
    
    cell_group <- colData(cds)[, group_cells_by]
    names(cell_group) = colnames(cds)
    
    exprs_mat <- lapply(gene_ids,
                        function(x){
                            m <- t(as.matrix(normalized_counts(cds)[x,]))
                            m <- reshape2::melt(m)
                            colnames(m) <- c("Cell", "Gene", "Expression")
                            m$Gene <- as.character(m$Gene)
                            m$Group <- cell_group[m$Cell]
                            m <- m %>% filter(!is.na(Group))
                        })
    
    ExpVal <- lapply(exprs_mat,
                     function(x){
                         v <- x %>% group_by(Group, Gene) %>% 
                             summarise(mean = mean(log(Expression + pseudocount)), 
                                       percentage = sum(Expression > lower_threshold)/length(Expression))
                         v$mean <- ifelse(v$mean < scale_min, scale_min, 
                                          v$mean)
                         v$mean <- ifelse(v$mean > scale_max, scale_max, 
                                          v$mean)
                         v$Gene <- fData(cds)[v$Gene, "gene_short_name"]
                         v$Group <- factor(v$Group, levels = cell_factor_levels)
                         v
                     })
    
    res <- lapply(ExpVal,
                  function(x){
                      r <- reshape2::dcast(x[, 1:4], Group ~ Gene, value.var = 'mean')
                      group_id <- r[, 1]
                      r <- r[, -1]
                      row.names(r) <- group_id
                      r
                  })
    
    if (ordering_type == "cluster_row_col") {
        ph <- lapply(res,
                     function(x){
                         row_dist <- stats::as.dist((1 - stats::cor(t(x)))/2)
                         row_dist[is.na(row_dist)] <- 1
                         
                         col_dist <- stats::as.dist((1 - stats::cor(x))/2)
                         col_dist[is.na(col_dist)] <- 1
                         
                         p <- pheatmap::pheatmap(x, useRaster = TRUE, cluster_cols = TRUE, 
                                                 cluster_rows = TRUE, show_rownames = FALSE, show_colnames = FALSE, 
                                                 clustering_distance_cols = col_dist, clustering_distance_rows = row_dist, 
                                                 clustering_method = "ward.D2", silent = TRUE, filename = NA)
                         p$tree_col$order
                     })
        
        ExpVal <- mapply(function(x, y, z){
            x$Gene <- factor(x$Gene, levels = colnames(y)[z])
            x
        }, ExpVal, res, ph, 
        SIMPLIFY = FALSE)
        
    }
    else if (ordering_type == "maximal_on_diag") {
        ExpVal <- lapply(ExpVal, function(x){
            x %>% 
                dplyr::group_by(Gene) %>% 
                dplyr::mutate(max_value = max(mean),
                              is_max = max_value == mean) %>% 
                dplyr::ungroup()
        })
        
        group_ordering_df <- lapply(ExpVal, function(x){
            x %>% 
                dplyr::filter(is_max) %>%
                dplyr::group_by(Group) %>% 
                dplyr::summarize(num_genes = n())
        })
        
        ExpVal <- mapply(function(x, y) {
            dplyr::left_join(x, y, by = c("Group")) %>% 
                dplyr::arrange(dplyr::desc(num_genes), 
                               dplyr::desc(max_value))
        }, ExpVal, group_ordering_df, 
        SIMPLIFY = F)
        
        ##Commented out to preserve factor levels for large plot. gene ordering
        ##Still makes it possible to get some diagonal line.
        # ExpVal <- lapply(ExpVal,
        #                  function(x){
        #                    x$Group <- factor(x$Group,
        #                                      levels = x %>% 
        #                                        dplyr::select(Group) %>% 
        #                                        dplyr::distinct() %>% 
        #                                        dplyr::pull(Group))
        #                    x
        #                  })
        
        gene_ordering <- lapply(ExpVal,
                                function(x){
                                    x %>% 
                                        dplyr::filter(is_max) %>% 
                                        dplyr::group_by(Gene) %>% 
                                        dplyr::slice_head(n = 1) %>% 
                                        dplyr::arrange(Group, dplyr::desc(max_value)) %>% 
                                        dplyr::pull(Gene)
                                })
        
        ExpVal <- mapply(function(x, y){
            x$Gene <- factor(x$Gene, levels = y)
            x
        }, ExpVal, gene_ordering,
        SIMPLIFY = F)
    }
    else if (ordering_type == "none") {
        ExpVal <- mapply(function(x, y) {
            x$Gene <- factor(x$Gene, levels = y)
            x
        }, ExpVal, markers,
        SIMPLIFY = FALSE)
    }
    
    ExpVal <- bind_rows(ExpVal, .id = 'marker_group')
    
    g <- ggplot(ExpVal, aes(y = Group, x = Gene)) +
        geom_point(aes(colour = mean, size = percentage)) +
        viridis::scale_color_viridis(name = "log(mean + 0.1)") +
        scale_size(name = "percentage", range = c(0,
                                                  max.size)) +
        facet_grid(~marker_group, scales = 'free', space = 'free') +
        monocle3:::monocle_theme_opts() + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
        theme_my() + 
        ylab(group_cells_by)
    
    g
}

get_top_markers <- function(subset, 
                            group_by = c('ct_cluster', 'celltype'), 
                            ref = 1000,
                            n = 4,
                            rank_by = c('pseudo_R2', 'specificity')){
    
    group_by <- match.arg(group_by, c('ct_cluster', 'celltype'))
    
    rank_by <- match.arg(rank_by, c('pseudo_R2', 'specificity'))
    
    markers <- top_markers(subset, group_cells_by = group_by, 
                           reference_cells = ref)
    
    markers <- markers %>%
        filter(fraction_expressing >= 0.10) %>%
        group_by(cell_group) %>%
        top_n(n, pseudo_R2)
    
    unique(markers %>% pull(gene_id))
}

prepare_organ_data <- function(subset, levels, group = c('ct_cluster',
                                                         'celltype')){
    group <- match.arg(group, c('ct_cluster',
                                'celltype'))
    dt <- subset@colData %>% as.data.table()
    dt[,
       by = .(ct,
              organ),
       .N,
       env = list(ct = group)
    ][, by = ct,
      freq := N / sum(N) * 100,
      env = list(ct = group)
    ][,
      let(ct = factor(ct, levels = levels),
          tmp = 'Organ'),
      env = list(ct = group)][]
}

plot_org <- function(dt, group = c('ct_cluster',
                                   'celltype')){
    group = match.arg(group, c('ct_cluster', 'celltype'))
    
    ggplot(dt) +
        geom_tile(aes(organ, !!ensym(group), fill = freq),
                  col = "white", 
                  width = 0.75) +
        scale_fill_gradient(low = "ivory2", high = "red") +
        facet_wrap(~tmp) +
        theme_my() + 
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = 'none'
        ) +
        ggtitle('% LN - Skin')
}
