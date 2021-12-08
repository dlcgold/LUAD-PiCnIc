library("xlsx")
library(TRONCO)

hypo_analysis <- function(LUAD){
    LUAD_less <- events.selection(LUAD, 
                                  filter.freq = .05);
    oncoprint(LUAD_less);
    del <- consolidate.data(LUAD_less);
    # if(length(del[["indistinguishable"]]) > 0){
    #     for (i in 1:length(del[["indistinguishable"]])) {
    #         type <- del[["indistinguishable"]][[i]][[1]];
    #         beg <- as.integer(length(del[["indistinguishable"]][[i]]) / 2) + 1;
    #         for (j in beg:length(del[["indistinguishable"]][[i]])){
    #             LUAD_less <- delete.event(LUAD_less, 
    #                                       gene = del[["indistinguishable"]][[i]][[j]], 
    #                                       type = type);
    #             
    #         }
    #     }
    # }
    oncoprint(LUAD_less);
    alterations = events.selection(as.alterations(LUAD_less),
                                   filter.freq = .05);
    gene.hypotheses = c('KRAS', 'TP53', 'TTN', 'PYR2', 'FLG', 'SPTA1', 'LRP1B');
    LUAD.clean = events.selection(LUAD_less,
                                  filter.in.names = c(as.genes(alterations), 
                                                      gene.hypotheses),
                                  filter.freq = .05);
    LUAD.clean = annotate.description(LUAD.clean,
                                      'Lung cancer data from GDC portal
                                    (selected events)');
    ## TODO understand why events with 0%
    oncoprint(LUAD.clean,
              gene.annot = list(priors = gene.hypotheses),
              sample.id = TRUE);
                        #LUAD_hypo <- nhypotheses(LUAD);
                                        #print(LUAD_hypo);
                                        # RANDOM hypotesis
    LUAD.hypo = hypothesis.add(LUAD.clean,
                               'TP53 xor KRAS',
                               XOR('TP53', 'KRAS'));
    oncoprint(events.selection(LUAD.hypo,
                               filter.in.names = c('TP53', 'KRAS')),
              font.row = 8,
              ann.hits = FALSE);
    LUAD.hypo = hypothesis.add(LUAD.clean,
                               'BRAF xor KRAS',
                               XOR('BRAF', 'KRAS'));
    oncoprint(events.selection(LUAD.hypo,
                               filter.in.names = c('BRAF', 'KRAS')),
              font.row = 8,
              ann.hits = FALSE);
    # LUAD.hypo = hypothesis.add(LUAD.clean,
    #                            'ERBB2 xor KRAS',
    #                            XOR('ERBB2', 'KRAS'));
    # oncoprint(events.selection(LUAD.hypo,
    #                            filter.in.names = c('ERBB2', 'KRAS')),
    #           font.row = 8,
    #           ann.hits = FALSE);
    # LUAD.hypo = hypothesis.add(LUAD.clean,
    #                            'ERBB2 xor EGRF',
    #                            XOR('ERBB2', 'EGRF'));
    # oncoprint(events.selection(LUAD.hypo,
    #                            filter.in.names = c('ERBB2', 'EGRF')),
    #           font.row = 8,
    #           ann.hits = FALSE);
  
    
    LUAD.hypo <- hypothesis.add.homologous(LUAD.hypo);
    LUAD <- hypothesis.add.group(LUAD.clean, 
                                 XOR,
                                 group = c('KRAS', 
                                           'EGFR', 
                                           'ALK',
                                           'ERBB2',
                                           'BRAF'));
    oncoprint(LUAD.hypo, 
              gene.annot = list(priors = gene.hypotheses), 
              sample.id = TRUE,
              font.row = 10, 
              font.column = 5, 
              cellheight = 15, 
              cellwidth = 4);
    npatt <<- npatterns(LUAD);
    nhypo <<- nhypotheses(LUAD);
    
    tronco.pattern.plot(LUAD,
                        group = as.events(LUAD, 
                                          genes=c('TP53', 
                                                  'KRAS')),
                        to = c('LRP1B', 
                               'Missense_Mutation'),
                        legend.cex=0.8,
                        label.cex=1.0);
    
    tronco.pattern.plot(LUAD,
                        group = as.events(LUAD, 
                                          genes=c('TP53', 
                                                  'KRAS')),
                        to = c('LRP1B', 
                               'Missense_Mutation'),
                        legend.cex=0.8,
                        label.cex=1.0, 
                        mode = 'circos');
    
    model.capri <- tronco.capri(LUAD.hypo, 
                                boot.seed = 12345, 
                                nboot = 5);
    model.boot <- tronco.bootstrap(model.capri, 
                                   nboot = 3, 
                                   cores.ratio = 0);
    boot <<- as.bootstrap.scores(model.boot);
    info <<- as.data.frame(as.parameters(model.capri));
    hm <<- has.model(model.capri);
    marginal.prob <<- as.marginal.probs(model.capri);
    joint.prob <<- as.joint.probs(model.capri, 
                                  models='capri_bic');
    conditional.prob <<- as.conditional.probs(model.capri, 
                                              models='capri_bic');
    confidence <<- as.confidence(model.capri, 
                                 conf = c('tp', 'pr', 'hg'));
    advantages <<- as.selective.advantage.relations(model.capri);
    
    model.boot <- tronco.bootstrap(model.boot, 
                                   nboot = 3, 
                                   cores.ratio = 0, 
                                   type = 'statistical');
    scoreboot <<- as.bootstrap.scores(model.boot);
    tronco.plot(model.boot,
                fontsize = 12,
                scale.nodes = .6,
                confidence=c('sb', 'npb'),
                height.logic = 0.25,
                legend.cex = .35,
                pathways = list(priors= gene.hypotheses),
                label.edge.size=10);
    
    pheatmap(keysToNames(model.boot, 
                         as.confidence(model.boot, 
                                       conf = 'sb')$sb$capri_aic) * 100,
             main = 'Statistical bootstrap scores for AIC model',
             fontsize_row = 6,
             fontsize_col = 6,
             display_numbers = TRUE,
             number_format = "%d");
    model.boot <<- tronco.kfold.eloss(model.boot);
    model.boot <<- tronco.kfold.prederr(model.boot, 
                                       runs = 2,
                                       cores.ratio = 0);
    model.boot <- tronco.kfold.posterr(model.boot,
                                       runs = 2,
                                       cores.ratio = 0);
    tab <<- tabular(model.boot, 
                    'capri_bic');
    tronco.plot(model.boot,
                fontsize = 12,
                scale.nodes = .6,
                confidence = c('npb', 'eloss', 'prederr', 'posterr'),
                height.logic = 0.25,
                legend.cex = .35,
                pathways = list(priors= gene.hypotheses),
                label.edge.size=10);
    return(LUAD.hypo);
}

tabular <- function(obj, M){
    tab <- Reduce(
        function(...) merge(..., all = TRUE),
        list(as.selective.advantage.relations(obj, models = M),
             as.bootstrap.scores(obj, models = M),
             as.kfold.prederr(obj, models = M),
             as.kfold.posterr(obj,models = M)))
    ## merge reverses first with second column
    tab <- tab[, c(2,1,3:ncol(tab))]
    tab <- tab[order(tab[, paste(M,
                                 '.NONPAR.BOOT',
                                 sep='')],
                     na.last = TRUE,
                     decreasing = TRUE), ]
    return(tab)
}

