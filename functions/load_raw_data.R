#Load Raw Data function
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