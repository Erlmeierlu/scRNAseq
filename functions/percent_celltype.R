count_groups <- function(cds, group_col = c('raw_celltype',
                                            'celltype')){
    group_col <- match.arg(group_col, c('raw_celltype',
                                        'celltype'))
    
    dt <- as.data.table(cds@colData)[treatment.agg != 'Undefined']
    dt[,
       by = .(experiment,
              organ,
              treatment.agg,
              sex,
              batch,
              ct),
       .N,
       env = list(ct = group_col)]
}

expand_groups <- function(dt, group_col = c('raw_celltype',
                                            'celltype')){
    group_col <- match.arg(group_col, c('raw_celltype',
                                        'celltype'))
    #Expand DataTable 
    exp <- dt[,
              CJ(experiment,
                 organ,
                 treatment.agg,
                 sex,
                 batch,
                 ct,
                 unique = T),
              env = list(ct = group_col)]
    
    #Join with expanded DT on cols 
    dt[exp,
       on = c('experiment',
              'organ',
              'treatment.agg',
              'sex',
              'batch',
              group_col)] %>% 
        droplevels()
}

calculate_percentages <- function(dt, 
                                  group_col = c('raw_celltype', 
                                                'celltype')){
    group_col <- match.arg(group_col, 
                           c('raw_celltype', 
                             'celltype'))
    
    setnafill(dt, fill = 0, cols = 'N')
    
    dt[,
       let(cells = sum(N),
           pc.cells = N / sum(N) * 100),
       keyby = .(experiment,
                 organ,
                 treatment.agg,
                 sex,
                 batch)]
    dt[,
       grp := .GRP,
       by = .(experiment, organ, treatment.agg, ct),
       env = list(ct = group_col)]
    
    lut <- dt[,sum(N) == 0, by = grp][V1 == FALSE, grp]
    
    dt[,
       pc.box := fcase(grp %in% lut, 
                       pc.cells,
                       default = 100000)][]
}

calculate_statistics <- function(dt, group_by){
    
    if(any(group_by %notin% names(dt))) stop('All arguments in "group_by" must be a column in "dt"',
                                             call. = FALSE)
    
    apply_lm <- function(x){
        if(uniqueN(x[,treatment.agg]) < 2) return(data.table())
        
        catch_mod <- function(){
            tryCatch(
                mod <- lm(pc.cells ~ treatment.agg + sex, x),
                error = function(e){
                    mod <<- lm(pc.cells ~ treatment.agg, x)
                })
            
            out <- as.data.table(summary(mod)$coefficient, keep.rownames = 'end')
            out <- out[grepl("WT|cKO|NoT", end)]
            out[,end:=str_remove(end, 'treatment.agg')]
            
            new_names <- names(out) %>% 
                str_replace_all(" ", "_") %>% 
                str_remove("\\.")
            
            setnames(out, c(new_names[1:4], 'p_value'))
            start_val <- levels(x[,treatment.agg])[1]
            out[,start:=start_val][end != 'NoT']
        }
        
        stats <- catch_mod()
        x$treatment.agg <- relevel(x$treatment.agg, 'HDAC_WT')
        stats2 <- catch_mod()
        
        rbind(stats, stats2)
    }
    
    dt <- dt[pc.cells != 0]
    
    dt <- dt[, keyby = group_by, apply_lm(.SD)]
    dt[,
       # keyby = group_by, 
       adj_p := p.adjust(p_value, method = 'BH')]
    dt[,label := fcase(adj_p %between% c(0.05, 0.1) & adj_p != 0.1, '*',
                       adj_p %between% c(0.01, 0.05) & adj_p != 0.05, '**',
                       adj_p < 0.01, '***',
                       default = 'NS')][]
}

plot_percent <- function(dt,
                         group_by,
                         remove_junk = TRUE){
    
    if(remove_junk == TRUE) dt <- dt[!grepl('Unass|remo', ct),
                                     env = list(ct = group_by)]
    
    point_dat <- copy(dt)
    point_dat[pc.box == 100000, pc.cells := 100000]
    point_dat <- point_dat[pc.cells != 0]
    
    dt %>% 
        ggplot(aes(!!ensym(group_by), pc.cells)) +
        geom_boxplot(
            aes(fill = treatment.agg,
                y = pc.box), 
            outlier.shape = NA, 
            width = 0.75,
            position = position_dodge2(preserve = 'single')) +
        geom_point(data = point_dat,
                   aes(group = treatment.agg, col = batch, shape = sex),
                   size = 2, 
                   alpha = 0.8,
                   position = position_jitterdodge(jitter.width = 0.3, 
                                                   dodge.width = 0.75)) + 
        # geom_signif(data = stats,
        #             aes(xmin = start,
        #                 xmax = end,
        #                 annotations = label,
        #                 y_position = position,
        #                 ),
        #             show.legend = TRUE,
        #             manual = T) +
        scale_fill_grey(start = 0.4, name = "Treatment Group") +
        scale_shape_discrete(name = "Sex") +
        scale_color_discrete(name = "Batch") +
        scale_y_log10() +
        facet_wrap(organ ~ experiment, ncol = 2) +
        theme_my(strip.background = element_blank(),
                 strip.text.x = element_text(hjust = 0,
                                             size = 20),
                 axis.text.x = element_text(size = 20),
                 axis.text.y = element_text(size = 20),
                 legend.title = element_text(size = 20),
                 legend.text = element_text(size = 20)
        ) +
        xlab("Celltype") +
        ylab("% of Cells") +
        coord_cartesian(ylim = range(dt$pc.cells, na.rm = T) + c(-.25, .25))
}