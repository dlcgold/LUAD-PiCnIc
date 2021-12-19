# MODELS

## select from LUAD with min freq and apriori
LUAD.select <- select(LUAD, 
                      min_freq, 
                      unique(             
                        c(LUAD.raf,
                          LUAD.enrich,
                          unlist(LUAD.mutex))))
LUAD.select <- annotate.description(LUAD.select,
                                    'LUAD selection')

## oncoprint the selection
if(plot_verbose){
  oncoprint(LUAD.select, 
            legend.cex = .5,          
            cellwidth = 3,            
            cellheight = 10,
            gene.annot = pathway.list, 
            gene.annot.color = pathways.color, 
            sample.id = TRUE)
}

## consolidate dataset
del <- consolidate.data(LUAD.select, 
                        print = TRUE)
if(length(del[["indistinguishable"]]) > 0){
  for (i in 1:length(del[["indistinguishable"]])) {
    for (j in 1:nrow(del[["indistinguishable"]][[i]])){
      gene <- del[["indistinguishable"]][[i]][j,][2][[1]]
      type <- del[["indistinguishable"]][[i]][j,][1][[1]]
      LUAD.select <- delete.event(LUAD.select,
                                  gene = as.character(gene),
                                  type = as.character(type))
    }
  }
}


## add hypotheses

if(hypo_reload){
  LUAD.hypo <- LUAD.select
  
  ## first hypotheses from mutex (using only available genes)
  if (!is.null(LUAD.mutex)) {
    for (group in LUAD.mutex) {
      group <- group[group%in% as.genes(LUAD.hypo)]
      if(length(group) >= 2){
        print(group)
        LUAD.hypo <- hypothesis.add.group(LUAD.hypo,
                                          FUN = OR,
                                          group = group,
                                          dim.min = length(group))
      }
    }
  }
  
  
  ## Add hypothes for RAF, checking if we have the genes in LUAD.select
  LUAD.raf.subtype <- LUAD.raf[LUAD.raf%in% as.genes(LUAD.hypo)]
  LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
                                    FUN = XOR, 
                                    group = LUAD.raf.subtype, 
                                    dim.min = length(LUAD.raf.subtype)) 
  ## Add hypothes for entich, checking if we have the genes in LUAD.select
  LUAD.enrich.subtype <- LUAD.enrich[LUAD.enrich%in% as.genes(LUAD.hypo)]
  LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
                                    FUN = AND, 
                                    group = LUAD.enrich.subtype, 
                                    dim.min = length(LUAD.enrich.subtype)) 
  
  # ## then for MTOR group
  # LUAD.mtor.subtype <- LUAD.mtor[LUAD.mtor%in% as.genes(LUAD.hypo)]
  # LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
  #                                   FUN = OR, 
  #                                   group = LUAD.mtor.subtype, 
  #                                   dim.min = length(LUAD.mtor.subtype))
  
  ## add all the hypotheses related to homologou events
  LUAD.hypo <- hypothesis.add.homologous(LUAD.hypo)
  
  ## add annotation
  LUAD.hypo <- annotate.description(LUAD.hypo, 
                                    as.description(LUAD.select))
  
  
  ## first use of CAPRI
  LUAD.model <- tronco.capri(LUAD.hypo, 
                             boot.seed = 12345,
                             nboot = num_boot_iter)
  
  ## DAG of model with hypotheses
  if(plot_verbose){
    dev.off()
    tronco.plot(LUAD.model, 
                pathways = pathway.list,  
                edge.cex = 1.5,          
                legend.cex = .35, 
                scale.nodes = .5,        
                confidence = c('tp', 'pr', 'hg'), 
                pathways.color = pathways.color,  
                disconnected = F,        
                height.logic = .3)
  }
  
  ## random test with a set of forced hypotheses
  gene.hypotheses <- c('KRAS', 'BRAF', 'ATM', 'STK11')
  alterations <- events.selection(as.alterations(LUAD.select), 
                                  filter.freq = min_freq)
  LUAD.hypo.clean <- events.selection(LUAD.select,
                                      filter.in.names = c(as.genes(alterations), 
                                                          gene.hypotheses))
  LUAD.hypo.clean <- annotate.description(LUAD.hypo.clean,
                                          'LUAD forced hypos (selected events)')
  if(plot_verbose){
    oncoprint(LUAD.hypo.clean,
              gene.annot = list(priors = gene.hypotheses), 
              sample.id = TRUE)
    oncoprint(LUAD.hypo.clean, 
              gene.annot = list(priors = gene.hypotheses), 
              sample.id = TRUE,
              font.row=10,
              font.column=5,
              cellheight=15, 
              cellwidth=4)
  }
  
  ## save data
  save(LUAD.hypo, 
       file = "input/luadDefHypo.rda")
  save(LUAD.model, 
       file = "input/luadDefHypoModel.rda")
  
}else{
  LUAD.hypo <- loadRData("input/luadDefHypo.rda")
  LUAD.model <- loadRData("input/luadDefHypoModel.rda")
}

## dataframe with selective advanges, with fit probabilities, optimized
LUAD.hypo.model.selfit <- as.selective.advantage.relations(LUAD.model)

if(verbose){
  print("advatanges selection fit probabilities")
  print(LUAD.hypo.model.selfit)
}

## dataframe with selective advanges, with prima facie, full set of edge
LUAD.hypo.model.selpf <- as.selective.advantage.relations(LUAD.model,
                                                          type = "pf")

if(verbose){
  print("advatanges selection full set")
  print(LUAD.hypo.model.selpf)
}

## dataframe with selective advanges, with a subset of genes
## TODO make test with usefull subset of genes
LUAD.hypo.model.selsub <- as.selective.advantage.relations(LUAD.model,
                                                           events = as.events(LUAD.model, 
                                                                              genes = P53))
if(verbose){
  print("advatanges selection by pathway")
  print(LUAD.hypo.model.selsub)
}

## TODO add some graph regarding pattern
## such as these but working
## examples for hard exclusivity
if(plot_verbose){
  if(all_mut){
    tronco.pattern.plot(LUAD.model,
                        group = as.events(LUAD.model, genes=c('TP53', 
                                                              'ATM')),
                        to = c('KRAS', 
                               'Nonsense_Mutation'),
                        legend.cex=0.8,
                        label.cex=1.0)
  }else{
    tronco.pattern.plot(LUAD.model,
                        group = as.events(LUAD.model, genes=c('TP53', 
                                                              'ATM')),
                        to = c('KRAS', 
                               'Mutation'),
                        legend.cex=0.8,
                        label.cex=1.0)
  }
}

if(plot_verbose){
  if(all_mut){
    tronco.pattern.plot(LUAD.model,
                        group = as.events(LUAD.model, genes=c('TP53', 
                                                              'ATM')),
                        to = c('KRAS', 
                               'Nonsense_Mutation'),
                        legend.cex=0.8,
                        label.cex=1.0,
                        mode = "circos")
  }else{
    tronco.pattern.plot(LUAD.model,
                        group = as.events(LUAD.model, genes=c('TP53', 
                                                              'ATM')),
                        to = c('KRAS', 
                               'Mutation'),
                        legend.cex=0.8,
                        label.cex=1.0,
                        mode = "circos")
  }
}

## a first brutal plot after capri
if(plot_verbose){
  dev.off()
  tronco.plot(LUAD.model, 
              pathways = pathway.list,  
              edge.cex = 1.5,          
              legend.cex = .35,         
              scale.nodes = .6,        
              confidence = c('tp', 'pr', 'hg'), 
              pathways.color = pathways.color,  
              disconnected = F,        
              height.logic = .3)
}
