library(data.table)
library(fst)

dataDir <- ('data')
resDir <- ("results")
shinyDir <- ('dge-app')


# Load Data ---------------------------------------------------------------

res <- read_rds(file.path(resDir, "res.rds"))
setDT(res)

res[, ':='(start = fcase(grepl('ctrl', treatment), 'NoT',
                             default = 'WT'),
               end = fcase(grepl('KO', treatment), 'cKO',
                           default = 'WT'))]
# Create Data -------------------------------------------------------------
design <- list(experiment = res[, unique(experiment)],
               celltype = res[, unique(celltype)])

write_rds(design, file.path(shinyDir, 'design.rds'))

cols <- c('coef', 'AveExpr', 't', 'B', 'direction')
res[, 
    keyby = .(celltype, organ, experiment), 
    write_fst(.SD, 
              path = file.path(shinyDir, 'data', paste0(
                  paste(.BY, collapse = '_'), 
                  '.fst')
                  ),
              compress = 0),
    .SDcols = !cols
    ]

list <- res[,.(rn, celltype, organ, experiment, 'direction2' = direction, treatment, t)]
write.fst(list, file.path(shinyDir, 'data', 'list_for_shiny.fst'))

