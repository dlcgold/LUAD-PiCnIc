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
## as in marker paper page 3
LUAD.raf <- c('KRAS', 'EGFR')
LUAD.megsa <- c('STK11', 'EGFR', 'ERBB2','U2AF1')
LUAD.megsa2 <- c('KRAS', 'EGFR', 'NF1','BRAF','MET')

## TODO does it make any sense?
if(plot_verbose){
  oncoprint(
    events.selection(LUAD,
                     filter.in.names = LUAD.raf), 
    title = paste("LAUD - RAF KRAS/EGFR exclusivity (knowledge prior)"),
    legend.cex = .3,
    font.row = 6,
    ann.hits = FALSE,
    cellheight = 10,
    cellwidth = 1,
    #silent = T,
    gene.annot = pathway.list,
    gene.annot.color = pathways.color,
    gtable = TRUE)
}


if(plot_verbose){
  oncoprint(
    events.selection(LUAD,
                     filter.in.names = LUAD.megsa), 
    title = paste("LAUD - MEGSA1 (knowledge prior)"),
    legend.cex = .3,
    font.row = 6,
    ann.hits = FALSE,
    cellheight = 10,
    cellwidth = 1,
    #silent = T,
    gene.annot = pathway.list,
    gene.annot.color = pathways.color,
    gtable = TRUE)
}

if(plot_verbose){
  oncoprint(
    events.selection(LUAD,
                     filter.in.names = LUAD.megsa2), 
    title = paste("LAUD - MEGSA2 (knowledge prior)"),
    legend.cex = .3,
    font.row = 6,
    ann.hits = FALSE,
    cellheight = 10,
    cellwidth = 1,
    #silent = T,
    gene.annot = pathway.list,
    gene.annot.color = pathways.color,
    gtable = TRUE)
}
