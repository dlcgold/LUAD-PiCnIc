# GROUPS EXCLUSIVITY

## print groups from mutex
if (verbose) {
  print(LUAD.mutex)
}

## oncoprint first two group of mutex file
if (plot_verbose) {
  grid.arrange(
    oncoprint(
      events.selection(LUAD,
                       filter.in.names = LUAD.mutex[[1]]),
      title = paste("LUAD - Mutex group 1"),
      legend.cex = .3,
      font.row = 6,
      cellheight = 10,
      ann.hits = FALSE,
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
## and as in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4478584/
LUAD.raf <- c('KRAS', 'EGFR')

## apriori mutual exlusions using MEGSA
## https://pubmed.ncbi.nlm.nih.gov/26899600/
## https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002162
LUAD.megsa <- c('STK11', 'EGFR', 'ERBB2', 'U2AF1')
LUAD.megsa2 <- c('KRAS', 'EGFR', 'NF1', 'BRAF', 'MET')

## oncoprint to see the mutual exclusions
if (plot_verbose) {
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
    gtable = TRUE
  )
}


if (plot_verbose) {
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
    gtable = TRUE
  )
}

if (plot_verbose) {
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
    gtable = TRUE
  )
}

par(.pardefault)