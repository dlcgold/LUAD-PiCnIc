# GROUPS EXCLUSIVITY

## print groups from mutex
if(verbose){
  print(LUAD.mutex)
}

## oncoprint first two group
if(plot_verbose){
  grid.arrange(
    oncoprint(
      events.selection(LUAD,
                       filter.in.names = LUAD.mutex[[1]]), 
      title = paste("LUAD - Mutex group 1"),
      legend.cex = .3,
      font.row = 6,
      ann.hits = FALSE, 
      cellheight = 10,
      silent = TRUE,    
      gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      gtable = TRUE
    )$gtable,
    oncoprint(
      events.selection(LUAD,
                       filter.in.names = LUAD.mutex[[2]]), 
      title = paste("LUAD - Mutex group 2"),
      legend.cex = .3,
      silent = TRUE,
      font.row = 6,
      ann.hits = FALSE,
      cellheight = 10,
      gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      gtable = TRUE
    )$gtable,
    ncol = 1    
  )
}

## apriori knowledge
## TODO now quasi-random genes, RAF genes and MTOR genes
LUAD.raf <- c('KRAS', 'EGFR')
LUAD.mtor <-  c('PIK3CA', 'STK11')

## TODO plot not work
if(plot_verbose && FALSE){
  grid.arrange(
    oncoprint(
      events.selection(LUAD,
                       filter.in.names = LUAD.mtor),
      title = paste("LUAD - MTOR exclusivity (knowledge prior)"),
      legend.cex = .3,
      font.row = 6,
      ann.hits = FALSE,
      cellheight = 10,
      cellwidth = 3,
      silent = T,
      gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      gtable = TRUE
    )$gtable,
    oncoprint(
      events.selection(LUAD,
                       filter.in.names = LUAD.raf), 
      title = paste("LAUD - RAF KRAS/NRAS/BRAF exclusivity (knowledge prior)"),
      legend.cex = .3,
      font.row = 6,
      ann.hits = FALSE,
      cellheight = 10,
      cellwidth = 3,
      silent = T,
      gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      gtable = TRUE
    )$gtable,
    ncol = 1 
  )
}
