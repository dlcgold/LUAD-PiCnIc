library("xlsx")
library(TRONCO)

hypo_analysis <- function(LUAD){
    LUAD_less <- events.selection(LUAD, filter.freq = .05);
    oncoprint(LUAD_less);
    del <- consolidate.data(LUAD_less);
    if(length(del[["indistinguishable"]]) > 0){
      for (i in 1:length(del[["indistinguishable"]])) {
          type <- del[["indistinguishable"]][[i]][[1]];
          start <- as.integer(length(del[["indistinguishable"]][[i]]) / 2) + 1;
          for (j in start:length(del[["indistinguishable"]][[i]])){
              dataset <- delete.event(LUAD_less, 
                                      gene = del[["indistinguishable"]][[i]][[j]], 
                                      type = type);
              
          }
      }
    }
    oncoprint(LUAD_less);
    alterations = events.selection(as.alterations(LUAD_less),
                                   filter.freq = .05);
    gene.hypotheses = c('KRAS', 'TP53', 'TTN', 'PYR2', 'FLG', 'SPTA1', 'LRP1B');
    LUAD.clean = events.selection(LUAD,
                                  filter.in.names = c(as.genes(alterations), 
                                                      gene.hypotheses));
    LUAD.clean = annotate.description(LUAD.clean,
                                      'Lung cancer data from Cbio portal
                                    (selected events)');
    ## TODO understand why events with 0%
    oncoprint(LUAD.clean,
              gene.annot = list(priors = gene.hypotheses),
              sample.id = TRUE);
    #LUAD_hypo <- nhypotheses(LUAD);
    #print(LUAD_hypo);
    
}
