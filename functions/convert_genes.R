#Convert Gene List Function
convert_gene_list <- function(gene_list, convert_to = c('mouse', 'human')) {
    
    convert_to <- match.arg(convert_to, c('mouse', 'human'))
    if(convert_to == 'mouse') convert_to <- 'mouse, laboratory'
    
    input <- fifelse(convert_to == 'human', 'mouse, laboratory', 'human')
    if(!exists('mhg_list')){
        mhg_list <<- fread("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt")
        setDT(mhg_list)
    }
    mouse_human_genes <- mhg_list[, .(`DB Class Key`,
                                      `Common Organism Name`,
                                      Symbol)]
    
    nomatch <- gene_list %notin% mouse_human_genes[`Common Organism Name` == input,
                                                   Symbol]
    
    keys <- mouse_human_genes[`Common Organism Name` == input
                              ][
                                  match(gene_list[!nomatch], Symbol),
                                  `DB Class Key`]
    
    mouse_human_genes <- mouse_human_genes[,
                                           dcast(
                                               .SD,
                                               `DB Class Key` ~ `Common Organism Name`,
                                               value.var = 'Symbol',
                                               fun.aggregate = paste,
                                               fill = NA,
                                               subset = .(`DB Class Key` %in% keys),
                                               drop = FALSE
                                           )]
    
    nomatch <- nomatch | gene_list %notin% mouse_human_genes[, 
                                                             input, 
                                                             env = list(input = input)]

    mouse_human_genes[is.na(human), 
                      human := str_to_upper(`mouse, laboratory`)]
    mouse_human_genes[is.na(`mouse, laboratory`), 
                      `mouse, laboratory` := str_to_title(human)]
    
    mouse_human_genes <- mouse_human_genes[match(gene_list, 
                                                 mouse_human_genes[,input,
                                                                   env = list(input = input)]
                                                 )]
    
    setnames(mouse_human_genes, convert_to, 'out')
    
    func <- fifelse(convert_to == 'human', 'str_to_upper', 'str_to_title') %>% 
        match.fun()

    mouse_human_genes[is.na(out), out := func(gene_list[nomatch])]

    mouse_human_genes[,out]
}
