library(S4Vectors)
library(shiny)
library(stringr)
library(DT)
library(enrichR)
library(ggplot2)
library(ggrepel)
library(shinyWidgets)
library(data.table)
library(memoise)
library(fgsea)
library(viridis)
library(ggsignif)
library(fst)


# Functions ---------------------------------------------------------------
#personal theme
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
        ),
        ... 
    )
}

enrichrGetGenesets <- function(databases) {
    lapply(databases, function(dbx) {
        cat(
            paste(
                "Fetching",
                dbx,
                "dataset from: maayanlab.cloud/Enrichr ...\n"
            )
        )
        fpath <-
            paste0(
                "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=",
                dbx
            )
        fhandle <- file(fpath)
        dblines <- tryCatch({
            readLines(con = fhandle)
        }, error = function(e) {
            message(e, "\nFailed reading database: ", dbx)
            NULL
        })
        close(fhandle)
        if (is.null(dblines)) {
            return(list())
        } else {
            res <- strsplit(dblines, "\t")
            names(res) <- sapply(res, function(x)
                x[1])
            res <-
                lapply(res, function(x)
                    x[3:length(x)])
            return(res)
        }
    })
}

cacheGenesets <- memoise(enrichrGetGenesets,
                         cache = cachem::cache_disk("./cached_data"))

my_fgsea <- function(stats, pathways) {
    y <- lapply(stats,
                function(x) {
                    lapply(
                        pathways,
                        fgsea,
                        stats = x,
                        minSize = 1,
                        nproc = 1
                    )
                })
    rbindlist(lapply(y,
                     rbindlist,
                     idcol = "Database"),
              idcol = "cell")
}

plot_volcano <- function(object,
                         intercept = 1,
                         highlight = NULL) {
    object %>% 
        ggplot(aes(logFC, p_transformed)) +
        geom_point(aes(fill = direction2),
                   shape = 21,
                   alpha = 0.5) +
        ylab("-log10(P.Value)") +
        facet_wrap( ~ treatment, scales = "free") +
        theme_my() +
        theme(
            panel.background = element_rect(fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                size = 1
            ),
            axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.position = "bottom",
            legend.key = element_blank(),
            axis.text.x = element_text(angle = 0, hjust = 0.5)
        ) +
        scale_fill_manual(
            values = c("#C8466D", "#3F60AE", "#E6E6E6"),
            name = "",
            limits = c("up", "down", "NS"),
            labels = c("Up", "Down", "Unaffected")
        ) +
        scale_colour_manual(
            values = c("white", "yellow", "black"),
            name = "",
            limits = c("up", "down", "NS"),
            labels = c("Up", "Down", "Unaffected")
        ) +
        geom_vline(
            xintercept = intercept,
            color = "black",
            linetype = "dotted",
            size = 0.5
        ) +
        geom_vline(
            xintercept = -intercept,
            color = "black",
            linetype = "dotted",
            size = 0.5
        ) +
        geom_point(
            data = object[rn %in% highlight],
            fill = "yellow",
            shape = 21,
            alpha = 1,
            size = 2.5
        ) +
        geom_label_repel(
            data = object[rn %in% highlight],
            aes(
                label = rn,
                fill = direction2,
                col = direction2
            ),
            size = 3,
            label.padding = 0.15,
            box.padding = 0.15,
            max.time = 5,
            show.legend = F
        )
    
}

plot_enr <- function(df) {
    
    ggplot(df, aes(treatment, Term)) +
        geom_point(aes(
            col = log(Odds.Ratio),
            size = pmin(5,-log10(Adjusted.P.value))
        )) +
        theme_my() +
        facet_grid(Database ~ celltype + organ + experiment,
                   scales = "free",
                   space = "free") +
        scale_x_discrete(drop = F) +
        scale_size_continuous(range = c(0, 5)) + guides(x = guide_axis(angle = 90)) +
        scale_color_gradientn(
            colors = c("#EEEEE0", "#DC3220"),
            values = scales::rescale(c(min(
                log(df$Odds.Ratio - 0.0001)
            ), max(
                log(df$Odds.Ratio + 0.0001)
            ))),
            oob = scales::squish
        )
}

plot_gene_plot <- function(data, stats) {
    ggplot(data, aes(treatment, NormExpr)) +
        geom_boxplot(
            size = 0.5,
            outlier.shape = NA,
            linewidth = 0.8,
            fill = "ivory",
            width = 0.5
        ) +
        geom_jitter(
            aes(fill = log2(n), shape = sex),
            width = 0.1,
            size = 2.5,
            height = 0,
            col = "black",
            show.legend = c(shape = T)
        ) +
        scale_shape_manual(values = c(21, 24)) +
        geom_signif(
            data = stats[label != "NS"],
            aes(
                xmin = start,
                xmax = end,
                annotations = label,
                y_position = position
            ),
            manual = T,
            tip_length = 0
        ) +
        theme_my() +
        theme(
            panel.background = element_rect(fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                size = 1
            ),
            axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.key = element_blank(),
            axis.text.x = element_text(angle = 0, hjust = 0.5)
        ) +
        ggtitle(data$rn) +
        scale_fill_viridis(
            name = "log2(Number Of Cells)",
            discrete = F,
            option = "H",
            begin = 0.2,
            end = 1
        ) 
}

plot_fgsea <- function(df) {
    
    ggplot(df, aes(treatment, pathway)) +
        geom_point(aes(col = NES, size = pmin(5,-log10(padj)))) +
        facet_grid(Database ~ celltype + organ + experiment,
                   scales = "free_y",
                   space = "free") +
        theme_my() +
        scale_color_gradient2(
            low = "#005AB5",
            mid = "#EEEEE0",
            high = "#DC3220",
            midpoint = 0
        )
}


# Load Data ---------------------------------------------------------------
shinyDir <- ('dge-app')
dataDir <- ('data')

#load data
design.res <- read_rds(file.path(shinyDir, 'design.rds'))

# UI Section --------------------------------------------------------------

ui <- fluidPage(
    tags$head(tags$style(
        HTML(
            "
                .custom-well {
                background-color: rgba(248, 248, 248, 0.85);
                position: absolute;
                z-index: 100;
                display: flex;
                justify-content: center;
                align-items: center;
                padding-top: 4px;
                padding-bottom: 4px;
                padding-left: 2.5px;
                padding-right: 2.5px;
                border-radius: 5px;
                border: 1px solid
                                }"
        )
    )),
    
    titlePanel("Volcano Plots"),
    fluidRow(
        sidebarPanel(width = 2,
                     h3("Plot Settings"),
                     selectizeInput(
                         "celltype",
                         "Celltype:",
                         choices = design.res$celltype,
                         multiple = FALSE
                     ),
                     selectInput(
                         "organ",
                         "Organ:",
                         choices = c("Lymphnode" = "LN",
                                     "Skin"),
                         multiple = FALSE,
                         selectize = FALSE,
                         size = 2
                     ),
                     selectInput(
                         "experiment",
                         "Experiment:",
                         choices = design.res$experiment,
                         multiple = FALSE,
                         selectize = FALSE,
                         size = 2
                     ),
                     numericInput(
                         "xintercept",
                         "logFC threshold:",
                         value = 1,
                         min = 0.3,
                         max = 3,
                         step = 0.1
                     ),
                     textInput(
                         "rn",
                         "Highlight genes:",
                         placeholder = "Gzma, Trac, Prf1"
                     ),
                     actionBttn(
                         "plotter",
                         "Generate",
                         color = "royal",
                         style = "fill",
                         block = T
                     )
        ),
        column(10,
               tabsetPanel(
                   tabPanel(
                       "Plot",
                       column(
                           7,
                           div(
                               style = "position:relative",
                               plotOutput(
                                   "Volcano_plot",
                                   brush = "plot_brush",
                                   hover = hoverOpts("plot_hover", 
                                                     delay = 10, 
                                                     delayType = "debounce"),
                                   click = "plot_click"
                               ),
                               uiOutput("hover_info")
                           ),
                           verbatimTextOutput("brush_info",
                                              placeholder = TRUE),
                           radioButtons(
                               "filetype",
                               "Select filetype:",
                               choices = c("pdf", "png"),
                               selected = "pdf"
                           ),
                           downloadLink("save_plot",
                                        "Save Plot")
                       ),
                       column(
                           5,
                           plotOutput("gene_plot"),
                           column(5,
                                  textInput("gene_input",
                                            label = NULL,
                                            placeholder = "Enter gene or click volcano plot")
                           ), 
                           column(2,
                                  actionButton(
                                      "plot_gene_input",
                                      label = "plot gene")
                           ),
                           column(12,
                                  downloadLink("save_gene_plot",
                                               "Save Plot"))
                       )
                   ),
                   tabPanel("Table", dataTableOutput("table"),
                            fluidRow(column(
                                12,
                                markdown("**Enter Filename**"),
                                fluidRow(column(
                                    4,
                                    textInput("filename",
                                              NULL,
                                              placeholder = "my cool genes")
                                ),
                                column(
                                    8,
                                    downloadButton('dl_table', "Export table")
                                ))
                            )))
               ))
    ),
    titlePanel("Enrichment Analysis"),
    fluidRow(
        sidebarPanel(width = 2,
                     h3("Enrichment Settings"),
                     checkboxInput("advanced",
                                   "Advanced Mode"),
                     selectizeInput(
                         "databases",
                         "Choose Databases:",
                         selected = c(
                             "KEGG_2019_Mouse",
                             "WikiPathways_2019_Mouse",
                             "GO_Molecular_Function_2023",
                             "GO_Biological_Process_2023",
                             "TRANSFAC_and_JASPAR_PWMs",
                             "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
                             "TRRUST_Transcription_Factors_2019",
                             "MSigDB_Hallmark_2020"
                         ),
                         choices = c("Type to search" = "",
                                     listEnrichrDbs()$libraryName),
                         multiple = TRUE
                     ),
                     conditionalPanel(
                         condition = "input.advanced==1",
                         radioButtons(
                             "enr_method",
                             "Enrichment Method",
                             choices = c("enrichR", "fgsea"),
                             selected = "enrichR"
                         ),
                         selectizeInput(
                             "ct_enr",
                             "Celltype(s):",
                             choices = design.res$celltype,
                             selected = design.res$celltype,
                             multiple = TRUE
                         ),
                         selectizeInput(
                             "organ_enr",
                             "Organ(s):",
                             choices = c("Lymphnode" = "LN",
                                         "Skin"),
                             multiple = TRUE,
                             selected = c("Lymphnode" = "LN",
                                          "Skin")
                         ),
                         selectizeInput(
                             "exp_enr",
                             "Experiment(s):",
                             choices = design.res$experiment,
                             multiple = TRUE,
                             selected = design.res$experiment
                         )
                         
                     ),
                     actionBttn(
                         "enrich",
                         "Enrich Genes",
                         color = "royal",
                         style = "fill",
                         block = T
                     )
                     
        ),
        mainPanel(width = 10, tabsetPanel(
            tabPanel(
                "Plot",
                plotOutput("enrichplot"),
                downloadLink("save_plot_enr",
                             "Save Plot")
            ),
            tabPanel("Table",
                     dataTableOutput("enrichres"),
                     fluidRow(column(
                         12,
                         markdown("**Enter Filename**"),
                         fluidRow(column(
                             4,
                             textInput("filename_enr",
                                       NULL,
                                       placeholder = "my enrichr results")
                         ),
                         column(
                             8,
                             downloadButton('dl_table_enr', "Export table")
                         ))
                     )))
        ))
    )
)

# Server Section ----------------------------------------------------------

server <- function(input, output, session) {
    
    object <- eventReactive(input$plotter, {
        res.pb <- read.fst(file.path(shinyDir,
                                     "data",
                                     paste0(
                                         input$celltype,
                                         "_",
                                         input$organ,
                                         "_",
                                         input$experiment,
                                         ".fst"
                                         
                                     )),
                           as.data.table = TRUE
        )
        
        res.pb[, ":="(
            direction2 =
                fcase(
                    logFC > threshold() & adj.P.Val < 0.05,
                    "up",
                    logFC < -threshold() &
                        adj.P.Val < 0.05,
                    "down",
                    default = "NS"
                ),
            p_transformed = -log10(P.Value)
        )]
        
    })
    
    designx <- reactive({
        setDT(readRDS(file.path(dataDir, "design_pseudobulk.rds")),
              keep.rownames = "rn")[celltype == input$celltype &
                                        organ == input$organ &
                                        experiment == input$experiment]
        
    }) %>%  
        bindEvent(input$plotter)
    
    countsx <- reactive({
        read_fst(file.path(
            "data",
            paste0(
                "voom_counts_",
                input$celltype,
                "_",
                input$organ,
                "_",
                input$experiment,
                ".fst"
            )
        ),
        as.data.table = TRUE)
    }) %>%
        bindEvent(designx())
    
    gene_id <- reactiveVal()
    
    observeEvent(input$plot_click, {
        click <- input$plot_click
        point <-
            nearPoints(object(),
                       click,
                       threshold = 5,
                       maxpoints = 1)
        if(nrow(point) == 0) return(gene_id(NULL))
        gene_id(point$rn)
    })
    
    observeEvent(input$plot_gene_input, {
        point <- str_to_sentence(unlist(str_split(
            str_remove_all(str_replace_all(input$gene_input,
                                           ";",
                                           ","),
                           " "),
            ","
        )))
        if (str_count(point) == 0) return(gene_id(NULL))
        gene_id(point)
    })
    
    observeEvent(input$plotter, {
        gene_id(NULL)
    })
    
    gene_data <- reactive({
        gene <- gene_id()
        if(is.null(gene)) return(NULL)
        data <- countsx()[rn == gene]
        data <- melt(
            data,
            id.vars = "rn",
            variable.name = "sample",
            value.name = "NormExpr"
        )
        
        data[, ":=" (
            treatment = factor(
                str_extract(sample, "[:alpha:]+$"),
                levels = c("NoT", "WT", "cKO")
            ),
            sex = str_extract(sample, "\\w"),
            n = designx()[match(rn, data[, sample]), n]
        )]
        
        stats <- object()[rn == gene]
        tmp <- data[,.SD[which.max(NormExpr)], by = treatment]
        
        stats[, ":=" (
            position = fcase(
                treatment == "WT_vs_ctrl",
                1.075 * max(tmp[treatment %in% c("NoT", "WT"), NormExpr]),
                treatment == "KO_vs_WT",
                1.15 * max(tmp[treatment %in% c("WT", "cKO"), NormExpr]),
                treatment == "KO_vs_ctrl",
                1.24 * max(tmp[, NormExpr])
            ),
            label = fcase(
                adj.P.Val %between% c(0.01, 0.05),
                "*",
                adj.P.Val %between% c(0.001, 0.01),
                "**",
                adj.P.Val < 0.001,
                "***",
                default = "NS"
            )
        )]
        return(list(data = data, stats = stats))
    }) %>%
        bindEvent(gene_id())
    
    goi <- eventReactive(input$plotter, {
        str_to_sentence(unlist(str_split(
            str_remove_all(str_replace_all(input$rn,
                                           ";",
                                           ","),
                           " "),
            ","
        )))
    })
    
    threshold <- eventReactive(input$plotter, {
        input$xintercept
    })
    
    my_plot <- eventReactive(input$plotter, {
        plot_volcano(object(),
                     intercept = threshold(),
                     highlight = goi()) + 
            ggtitle(paste(
                input$celltype,
                input$organ,
                input$experiment
            ))
        
    })
    
    output$Volcano_plot <- renderPlot({
        my_plot()
    })
    
    output$brush_info <- renderPrint({
        if (!is.null(input$plot_brush)) {
            brush = brushedPoints(object(), input$plot_brush)
            
            cat("Selected:\n",
                paste(brush$rn, collapse = ", "))
        }
        
    })
    
    output$hover_info <- renderUI({
        hover <- input$plot_hover
        point <-
            nearPoints(
                object(),
                hover,
                threshold = 5,
                maxpoints = 1,
                addDist = TRUE
            )
        
        if (nrow(point) == 0)
            return(NULL)
        
        # calculate point position INSIDE the image as percent of total dimensions
        # from left (horizontal) and from top (vertical)
        left_pct <-
            (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
        top_pct <-
            (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
        
        # calculate distance from left and bottom side of the picture in pixels
        left_px <-
            hover$range$left + left_pct * (hover$range$right - hover$range$left)
        top_px <-
            hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
        
        # create style property for tooltip
        
        style <-
            paste0("left:", left_px + 2, "px; top:", top_px + 2, "px;")
        
        div(
            class = "custom-well",
            style = style,
            div(style = "height: 11px; line-height: 11px; text-align: center; font-size: 11px",
                p(HTML(
                    paste0("<b>", point$rn, "</b>")
                )))
        )
    })
    
    
    output$save_plot <- downloadHandler(
        filename = function() {
            paste0(
                "volcano",
                "_",
                input$celltype,
                "_",
                input$organ,
                "_",
                input$experiment,
                ".",
                input$filetype
            )
        },
        content = function(file) {
            if (input$filetype == "png") {
                png(file, width = 720)
            } else{
                pdf(file, width = 10)
            }
            plot(my_plot())
            dev.off()
        }
    )
    
    my_gene_plot <- reactive({
        
        plot_gene_plot(gene_data()$data, gene_data()$stats)
        
    }) %>% 
        bindEvent(req(gene_data()))
    
    output$gene_plot <- renderPlot(my_gene_plot())
    
    output$save_gene_plot <- downloadHandler(
        filename = function() {
            paste0(
                "cpm",
                "_",
                gene_data()$data$rn,
                "_",
                input$celltype,
                "_",
                input$organ,
                "_",
                input$experiment,
                ".pdf"
            )
        },
        content = function(file) {
            pdf(file)
            plot(my_gene_plot())
            dev.off()
        }
    )
    
    my_table <- reactive({
        x <- brushedPoints(object(), input$plot_brush)
        x[, .(gene_name = rn,
              logFC = formatC(logFC, digits = 3),
              p_value = P.Value,
              adjusted_pval = adj.P.Val,
              direction2)] 
    })
    
    output$table <- renderDataTable({
        formatSignif(
            datatable(my_table()[, direction2 := NULL]),
            columns = c("p_value",
                        "adjusted_pval"),
            digits = 2
        )
    })
    
    output$dl_table <- downloadHandler(
        filename = function() {
            fifelse(
                input$filename == "",
                paste0("Genes_",
                       Sys.Date(),
                       ".csv"),
                paste0(
                    str_replace_all(input$filename,
                                    " |\\.|\\\\|\\/|\\,|\\-", "_"),
                    ".csv"
                )
            )
        },
        content = function(fname) {
            write.csv(my_table(),
                      fname,
                      row.names = F)
        }
    )
    
    res_enrich <- reactive({
        
        if (!input$advanced) {
            
            res <- enrichr(my_table()$gene_name,
                           input$databases)
            res <- rbindlist(res, idcol = "Database")
            cols <- grep("Old", names(res))
            res <- res[, !..cols][order(Adjusted.P.value)]
        } else if (input$enr_method == "enrichR") {
            showModal(
                modalDialog(
                    "Gene Enrichment... This might take several minutes",
                    footer = NULL
                )
            )
            
            res <- read.fst(file.path("data",
                                      "list_for_shiny.fst"),
                            as.data.table = TRUE)[direction2 == "up" &
                                                      celltype %in% input$ct_enr &
                                                      organ %in% input$organ_enr &
                                                      experiment %in% input$exp_enr][, .(rn = paste(rn, collapse = ",")),
                                                                                     by = .(celltype,
                                                                                            organ,
                                                                                            experiment,
                                                                                            treatment)]
            
            res[, ":=" (tmp = lapply(str_split(res$rn, ","),
                                     function(x) {
                                         enrichr(x, input$databases) %>%
                                             rbindlist(idcol = "Database")
                                     }),
                        rn = NULL)]
            
            res <-
                res[, rbindlist(tmp), by = .(celltype,
                                             organ,
                                             experiment,
                                             treatment)]
            removeModal()
            
            cols <- grep("Old", names(res))
            res <- res[, !..cols][order(Adjusted.P.value)]
            
            dat <- res[Adjusted.P.value < 0.05 &
                           str_count(Genes, ";") >= 1][
                               res[, .I[frankv(-Odds.Ratio) < 3], #CHANGED TIES METHOD CHECK HERE IF ERROR
                                   by = .(celltype, organ, experiment, Database)]
                               $V1, Term
                           ]
            
            
            res <- res[Term %in% dat]
            res[, experiment := factor(res$experiment,
                                       levels = c("HDAC1",
                                                  "HDAC2"))]
            gc()
            return(res)
            
        } else {
            showModal(
                modalDialog(
                    "Gene Enrichment... This might take several minutes",
                    footer = NULL
                )
            )
            
            paths <-
                sapply(input$databases, cacheGenesets)
            
            res <- read.fst(file.path("data",
                                      "list_for_shiny.fst"),
                            as.data.table = TRUE)[celltype %in% input$ct_enr &
                                                      organ %in% input$organ_enr &
                                                      experiment %in% input$exp_enr]
            
            res <- res[, .(data = .(.SD)),
                       by = .(celltype,
                              organ,
                              experiment,
                              treatment),
                       .SDcols = !"direction2"]
            
            res[, data := lapply(res$data, \(x) setNames(x[,t], x[,rn]))]
            
            
            res <- setNames(res$data,
                            paste(res$celltype,
                                  res$organ,
                                  res$experiment,
                                  res$treatment,
                                  sep = "-"))
            
            
            res <-
                my_fgsea(res, paths)
            
            res[, c("celltype",
                    "organ",
                    "experiment",
                    "treatment") := tstrsplit(cell, "-")]
            
            removeModal()
            
            res.filtered <- res[padj < 0.05][res[, .I[frankv(padj, ties.method = "dense") < 6], 
                                                 by = .(celltype, organ, experiment, treatment, Database)]$V1, pathway]
            
            res <- res[pathway %in% res.filtered]
            
            res[, ":=" (treatment = factor(treatment,
                                           levels = c("WT_vs_ctrl", "KO_vs_ctrl", "KO_vs_WT")),
                        leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ","))]
            return(res)
        }
        
    }) %>%
        bindEvent(input$enrich)
    
    
    
    plot_width <-
        eventReactive(input$enrich, {
            6 + 1.5 * length(unique(res_enrich()$celltype)) *
                length(unique(res_enrich()$organ)) *
                length(unique(res_enrich()$experiment))
        })
    
    
    plot_height <- eventReactive(input$enrich, {
        2.25 + 4.75 * length(unique(res_enrich()$Database))
    })
    
    
    fgsea_width <-
        eventReactive(input$enrich, {
            6 + 1.5 * length(unique(res_enrich()$celltype)) *
                length(unique(res_enrich()$organ)) *
                length(unique(res_enrich()$experiment))
        })
    
    
    fgsea_height <- eventReactive(input$enrich, {
        2.5 + 5.5 * length(unique(res_enrich()$Database))
    })
    
    my_enrich <- reactive({
        if (!input$advanced) {
            plotEnrich(res_enrich()) + theme_my()
        } else if (input$enr_method == "enrichR") {
            plot_enr(res_enrich())
        } else {
            plot_fgsea(res_enrich())
        }
    }) %>%
        bindEvent(res_enrich())
    
    
    
    output$enrichplot <- renderPlot({
        my_enrich()
    })
    
    output$save_plot_enr <- downloadHandler(
        filename = function() {
            paste0("EnrichR_plot",
                   Sys.Date(),
                   ".pdf")
        },
        content = function(file) {
            if (!input$advanced) {
                pdf(file)
            } else if (input$enr_method == "enrichR") {
                pdf(file,
                    width = plot_width(),
                    height = plot_height())
            } else {
                pdf(file,
                    width = fgsea_width(),
                    height = fgsea_height())
            }
            plot(my_enrich())
            dev.off()
        }
    )
    
    
    
    output$enrichres <- renderDataTable({
        if (any(input$enr_method == "enrichR",
                !input$advanced)) {
            res_enrich() %>% 
                datatable() %>% 
                formatSignif(
                    columns = c(
                        "Adjusted.P.value",
                        "P.value",
                        "Odds.Ratio",
                        "Combined.Score"
                    ),
                    digits = 3
                )
        } else {
            
            res_enrich()[, !c("log2err", "ES")] %>% 
                datatable() %>% formatSignif(columns = c("padj",
                                                        "pval",
                                                        "NES"),
                                            digits = 3) 
        }
        
    })
    
    output$dl_table_enr <- downloadHandler(
        filename = function() {
            fifelse(
                input$filename_enr == "",
                fifelse(
                    input$enr_method == "enrichR",
                    paste0("EnrichR_",
                           Sys.Date(),
                           ".csv"),
                    paste0("fgsea_",
                           Sys.Date(),
                           ".csv")
                ),
                paste0(
                    str_replace_all(
                        input$filename_enr,
                        " |\\.|\\\\|\\/|\\,|\\-",
                        "_"
                    ),
                    ".csv"
                )
            )
        },
        content = function(fname) {
            write.csv(res_enrich(),
                      fname,
                      row.names = F)
        }
    )
}


# Run the application 
shinyApp(ui = ui, server = server)
