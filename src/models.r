#' Function to make model reconstruction for LUAD
#'
#' @param LUAD TRONCO object dataset
#' @param label label to identify subtype
#'
#' @return the model obtained with capri algorithm
model <- function(LUAD,
                  label) {
  ## model <- function(LUAD, gene.hypotheses, gene.sel, genes.compare, genes.to, label){
  
  ## select from LUAD with min freq, apriori knowledge and mutex genes
  LUAD.select <- select(LUAD,
                        min_freq,
                        unique(c(
                          LUAD.raf,
                          LUAD.megsa,
                          LUAD.megsa2,
                          unlist(LUAD.mutex)
                        )))
  LUAD.select <-
    annotate.description(LUAD.select, paste('LUAD', label, 'selection'))
  
  ## oncoprint of the selection
  if (plot_verbose) {
    oncoprint(
      LUAD.select,
      # legend.cex = .5,
      # cellwidth = 3,
      # cellheight = 10,
      gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      sample.id = TRUE
    )
  }
  
  ## consolidate dataset
  del <- consolidate.data(LUAD.select,
                          print = TRUE)
  ## remove ambiguous events
  if (length(del[["indistinguishable"]]) > 0) {
    for (i in seq_along(del[["indistinguishable"]])) {
      for (j in seq_len(nrow(del[["indistinguishable"]][[i]]))) {
        gene <- del[["indistinguishable"]][[i]][j,][2][[1]]
        type <- del[["indistinguishable"]][[i]][j,][1][[1]]
        LUAD.select <- delete.event(LUAD.select,
                                    gene = as.character(gene),
                                    type = as.character(type))
      }
    }
  }
  
  ## add hypotheses
  LUAD.hypo <- LUAD.select
  
  ## first hypotheses from mutex (using only available genes)
  if (!is.null(LUAD.mutex)) {
    for (group in LUAD.mutex) {
      group <- group[group %in% as.genes(LUAD.hypo)]
      if (length(group) >= 2) {
        print(group)
        LUAD.hypo <- hypothesis.add.group(
          LUAD.hypo,
          FUN = OR,
          group = group,
          dim.min = length(group)
        )
      }
    }
  }
  
  ## Add hypothes for RAF, checking if we have the genes in LUAD.select
  LUAD.raf.subtype <- LUAD.raf[LUAD.raf %in% as.genes(LUAD.hypo)]
  LUAD.hypo <- hypothesis.add.group(
    LUAD.hypo,
    FUN = XOR,
    group = LUAD.raf.subtype,
    dim.min = length(LUAD.raf.subtype)
  )
  
  ## Add hypothes for megsa, checking if we have the genes in LUAD.select
  LUAD.megsa.subtype <-
    LUAD.megsa[LUAD.megsa %in% as.genes(LUAD.hypo)]
  LUAD.hypo <- hypothesis.add.group(
    LUAD.hypo,
    FUN = OR,
    group = LUAD.megsa.subtype,
    dim.min = length(LUAD.megsa.subtype)
  )
  
  ## Add hypothes for megsa2, checking if we have the genes in LUAD.select
  LUAD.megsa2.subtype <-
    LUAD.megsa2[LUAD.megsa2 %in% as.genes(LUAD.hypo)]
  LUAD.hypo <- hypothesis.add.group(
    LUAD.hypo,
    FUN = OR,
    group = LUAD.megsa2.subtype,
    dim.min = length(LUAD.megsa2.subtype)
  )
  
  ## Added co-mutation(NF1, TP53) hypotesis as in marker paper
  if (label == 'proximal proliferative (PP, magnoid)') {
    LUAD.hypo <-
      hypothesis.add(LUAD.hypo, 'KRAS and STK11', AND('KRAS', 'STK11'))
  }
  
  ## Added co-mutation(NF1, TP53) hypotesis as in marker paper
  if (label == 'proximal inflammatory (PI, squamoid)') {
    LUAD.hypo <-
      hypothesis.add(LUAD.hypo, 'NF1 and TP53', AND('NF1', 'TP53'))
  }
  
  ## Add all the hypotheses related to homologous events
  LUAD.hypo <- hypothesis.add.homologous(LUAD.hypo)
  
  ## Edit annotation
  LUAD.hypo <- annotate.description(LUAD.hypo,
                                    as.description(LUAD.select))
  
  
  ## First use of CAPRI, with default parameters
  ## except for bootstrap iterations number
  LUAD.model <- tronco.capri(LUAD.hypo,
                             boot.seed = 42,
                             nboot = num_boot_iter)
  
  ## DAG of model with hypotheses
  if (plot_verbose) {
    tronco.plot(
      LUAD.model,
      pathways = pathway.list,
      edge.cex = 1.5,
      legend.cex = .35,
      scale.nodes = .3,
      confidence = c('tp', 'pr', 'hg'),
      pathways.color = pathways.color,
      disconnected = F,
      height.logic = .3,
    )
  }
  
  ## random test with a set of forced hypotheses
  ## gene.hypotheses <- c('KRAS', 'BRAF', 'ATM', 'STK11')
  
  ## TODO following lines are useless
  # alterations <- events.selection(as.alterations(LUAD.select),
  #                                 filter.freq = min_freq)
  # LUAD.hypo.clean <- events.selection(LUAD.select,
  #                                     filter.in.names = c(as.genes(alterations),
  #                                                         gene.hypotheses))
  # LUAD.hypo.clean <- annotate.description(LUAD.hypo.clean,
  #                                         paste(
  #                                           'LUAD forced hypos (selected events)',
  #                                           label))
  # if(plot_verbose){
  #   oncoprint(LUAD.hypo.clean,
  #             gene.annot = list(priors = gene.hypotheses),
  #             sample.id = TRUE)
  # }
  # if(plot_verbose){
  #   oncoprint(LUAD.hypo.clean,
  #             gene.annot = list(priors = gene.hypotheses),
  #             sample.id = TRUE,
  #             font.row=10,
  #             font.column=5,
  #             cellheight=5,
  #             cellwidth=1)
  # }
  
  ## save data
  ##save(LUAD.hypo,
  ##     file = "input/luadDefHypo.rda")
  ##save(LUAD.model,
  ##     file = "input/luadDefHypoModel.rda")
  
  
  ## dataframe with selective advantages, with fit probabilities, optimized
  LUAD.hypo.model.selfit <-
    as.selective.advantage.relations(LUAD.model)
  
  if (verbose) {
    print(paste("advatanges selection fit probabilities for", label))
    print(LUAD.hypo.model.selfit)
  }
  
  
  
  ## dataframe with selective advances, with prima facie, full set of edge
  LUAD.hypo.model.selpf <-
    as.selective.advantage.relations(LUAD.model,
                                     type = "pf")
  
  if (verbose) {
    print(paste("advatanges selection full set for", label))
    print(LUAD.hypo.model.selpf)
  }
  
  ## There are not hardexclusivity subgroups. Trust The Data, full stop!
  
  ## plot of final reconstruction model before statistical analysis
  if (plot_verbose) {
    tronco.plot(
      LUAD.model,
      pathways = pathway.list,
      edge.cex = 1.5,
      legend.cex = .35,
      scale.nodes = .6,
      confidence = c('tp', 'pr', 'hg'),
      pathways.color = pathways.color,
      disconnected = F,
      height.logic = .3,
    )
  }
  
  
  return(LUAD.model)
  
}
