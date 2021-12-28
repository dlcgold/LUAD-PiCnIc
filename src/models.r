# model(LUAD, gene.hypotheses, gene.sel, genes.compare, genes.to, 'all')
model <- function(LUAD, gene.hypotheses, gene.sel, genes.compare, genes.to, label){

  ## select from LUAD with min freq and apriori
  LUAD.select <- select(LUAD, 
                        min_freq, 
                        unique(             
                          c(LUAD.raf,
                            LUAD.megsa,
                            LUAD.megsa2,
                            unlist(LUAD.mutex))))
  LUAD.select <- annotate.description(LUAD.select, paste('LUAD', label, 'selection'))
  
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
  
  ## Add hypothes for megsa, checking if we have the genes in LUAD.select
  LUAD.megsa.subtype <- LUAD.megsa[LUAD.megsa%in% as.genes(LUAD.hypo)]
  LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
                                    FUN = OR, 
                                    group = LUAD.megsa.subtype, 
                                    dim.min = length(LUAD.megsa.subtype))
  
  ## Add hypothes for megsa2, checking if we have the genes in LUAD.select
  LUAD.megsa2.subtype <- LUAD.megsa2[LUAD.megsa2%in% as.genes(LUAD.hypo)]
  LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
                                    FUN = OR, 
                                    group = LUAD.megsa2.subtype, 
                                    dim.min = length(LUAD.megsa2.subtype))
  
  ## Add all the hypotheses related to homologou events
  LUAD.hypo <- hypothesis.add.homologous(LUAD.hypo)
  
  ## Add annotation
  LUAD.hypo <- annotate.description(LUAD.hypo, 
                                    as.description(LUAD.select))
  
  
  ## First use of CAPRI
  LUAD.model <- tronco.capri(LUAD.hypo, 
                             boot.seed = 42,
                             nboot = num_boot_iter)
  
  ## DAG of model with hypotheses
  if(plot_verbose){
    tronco.plot(LUAD.model, 
                pathways = pathway.list,  
                edge.cex = 1.5,          
                legend.cex = .35, 
                scale.nodes = .3,        
                confidence = c('tp', 'pr', 'hg'), 
                pathways.color = pathways.color,  
                disconnected = F,        
                height.logic = .3,)
  }
  
  ## random test with a set of forced hypotheses
  ## gene.hypotheses <- c('KRAS', 'BRAF', 'ATM', 'STK11')
  
  alterations <- events.selection(as.alterations(LUAD.select), 
                                  filter.freq = min_freq)
  LUAD.hypo.clean <- events.selection(LUAD.select,
                                      filter.in.names = c(as.genes(alterations), 
                                                          gene.hypotheses))
  LUAD.hypo.clean <- annotate.description(LUAD.hypo.clean,
                                          paste(
                                            'LUAD forced hypos (selected events)',
                                            label))
  if(plot_verbose){
    oncoprint(LUAD.hypo.clean,
              gene.annot = list(priors = gene.hypotheses), 
              sample.id = TRUE)
  }
  if(plot_verbose){
    oncoprint(LUAD.hypo.clean, 
              gene.annot = list(priors = gene.hypotheses), 
              sample.id = TRUE,
              font.row=10,
              font.column=5,
              cellheight=5, 
              cellwidth=1)
  }
  
  ## save data
  ##save(LUAD.hypo, 
  ##     file = "input/luadDefHypo.rda")
  ##save(LUAD.model, 
  ##     file = "input/luadDefHypoModel.rda")
  
  
  ## dataframe with selective advanges, with fit probabilities, optimized
  LUAD.hypo.model.selfit <- as.selective.advantage.relations(LUAD.model)
  
  if(verbose){
    print(paste("advatanges selection fit probabilities for", label))
    print(LUAD.hypo.model.selfit)
  }
  
  ## dataframe with selective advanges, with prima facie, full set of edge
  LUAD.hypo.model.selpf <- as.selective.advantage.relations(LUAD.model,
                                                            type = "pf")
  
  if(verbose){
    print(paste("advatanges selection full set for", label))
    print(LUAD.hypo.model.selpf)
  }
  
  
  # interjection
  gene.sel <- gene.sel[gene.sel%in% as.genes(LUAD)]
  ## dataframe with selective advanges, with a subset of genes
  ## TODO make test with usefull subset of genes
  LUAD.hypo.model.selsub <- as.selective.advantage.relations(LUAD.model,
                                                             events = as.events(LUAD.model, 
                                                                                genes = gene.sel))
  if(verbose){
    print(paste("advatanges selection by for", label))
    print(LUAD.hypo.model.selsub)
  }
  
  # TODO add some graph regarding pattern
  # such as these but working
  # examples for hard exclusivity
  # plots for presentation

  if(plot_verbose){
    tronco.pattern.plot(LUAD.hypo.clean,
                        group = as.events(LUAD.hypo.clean, genes=genes.compare),
                        to = genes.to,
                        legend.cex=0.8,
                        label.cex=1.0,
                        mode="barplot")
  }
  par(.pardefault)

  if(plot_verbose){

    tronco.pattern.plot(LUAD.hypo.clean,
                        group = as.events(LUAD.hypo.clean, genes=genes.compare),
                        to = genes.to,
                        legend.cex=0.8,
                        label.cex=1.0,
                        mode = "circos")
  }
  ## a first brutal plot after capri
  if(plot_verbose){
    tronco.plot(LUAD.model, 
                pathways = pathway.list,  
                edge.cex = 1.5,          
                legend.cex = .35,         
                scale.nodes = .6,        
                confidence = c('tp', 'pr', 'hg'), 
                pathways.color = pathways.color,  
                disconnected = F,        
                height.logic = .3,)
  }
  
  ## STATISTICS
  
  # STATISTICS
  ## non-parametric bootstrap
  LUAD.model <- tronco.bootstrap(LUAD.model,
                                 nboot = num_boot_iter,
                                 cores.ratio = .5)
  
  ## statistical bootstrap
  LUAD.model <- tronco.bootstrap(LUAD.model,
                                 type = "statistical",
                                 nboot = num_boot_iter,
                                 cores.ratio = .5)
  
  ## DAG of the model above
  if(plot_verbose){
    tronco.plot(LUAD.model, 
                pathways = pathway.list,  
                edge.cex = 1.5,          
                legend.cex = .5,         
                scale.nodes = .6,        
                confidence = c('tp', 'pr', 'hg'), 
                pathways.color = pathways.color,  
                disconnected = F,        
                height.logic = .3,)
    
  }
  
  ## plot of bootstrap scores
  ## TODO sometimes not work
  if(plot_verbose){
    ## first non-parametric
    pheatmap(keysToNames(LUAD.model,
                         as.confidence(LUAD.model,
                                       conf = 'npb')$npb$capri_bic) * 100,
             main = paste("non-parametric bootstrap scores for BIC model for", label),
             fontsize_row = 6,
             fontsize_col = 6,
             display_numbers = T,
             number_format = "%f"
    )
  }
  if(plot_verbose){
    pheatmap(keysToNames(LUAD.model,
                         as.confidence(LUAD.model,
                                       conf = 'npb')$npb$capri_aic) * 100,
             main = paste("non-parametric bootstrap scores for AIC model for",label),
             fontsize_row = 6,
             fontsize_col = 6,
             display_numbers = T,
             number_format = "%f"
    )
  }
  if(plot_verbose){
  ## then parametric ones
    pheatmap(keysToNames(LUAD.model,
                         as.confidence(LUAD.model,
                                       conf = 'sb')$sb$capri_bic) * 100,
             main = paste("non-parametric bootstrap scores for BIC model for",label),
             fontsize_row = 6,
             fontsize_col = 6,
             display_numbers = T,
             number_format = "%f"
    )
  }
  if(plot_verbose){
    pheatmap(keysToNames(LUAD.model,
                         as.confidence(LUAD.model,
                                       conf = 'sb')$sb$capri_aic) * 100,
             main = paste("non-parametric bootstrap scores for AIC model for",label),
             fontsize_row = 6,
             fontsize_col = 6,
             display_numbers = T,
             number_format = "%f"
    )
  }
  
  ## table with bootstrap scores
  boot_tab <- as.bootstrap.scores(LUAD.model)
  
  if(verbose){
    print(boot_tab)
  }
  
  ## kfold
  ## k-fold cross validation, prediction error for each parent set X
  LUAD.model <- tronco.kfold.eloss(LUAD.model)
  kfold_eloss <- as.kfold.eloss(LUAD.model)
  
  ## plot for every fold
  ## TODO make it work
  if(plot_verbose){
    vioplot(LUAD.model$kfold$capri_bic$eloss,
            LUAD.model$kfold$capri_aic$eloss,
            col = 'red',
            lty = 1, rectCol="gray",
            colMed = 'black',
            names = c('BIC', 'AIC'), 
            pchMed = 15, 
            horizontal = T)
    title(main = paste('Entropy loss \n LUAD tumors -', label))
  }
  
  ## k-fold cross validation, prediction error for each parent set X
  LUAD.model <- tronco.kfold.prederr(LUAD.model)
  kfold_pred <- as.kfold.prederr(LUAD.model)
  
  
  ## k-fold cross validation, posterior classification error for each edge
  LUAD.model <- tronco.kfold.posterr(LUAD.model)
  kfold_post <- as.kfold.posterr(LUAD.model)
  
  ## visualize a table with all edge statistics
  tab_bic <- tabular(LUAD.model, 'capri_bic')
  tab_aic <- tabular(LUAD.model, 'capri_aic')
  
  if(verbose){
    print(label)
    print("table with all edge statistics using capri_bic")
    print(tab_bic)
    print("table with all edge statistics using capri_aic")
    print(tab_aic)
  }
  
  ## save model
  save(LUAD.model, 
       file = paste("output/luadDefModel_", label, ".rda", sep=''))
  
  
  ## last DAG
  if(plot_verbose){
    tronco.plot(LUAD.model, 
              pathways = pathway.list,  
              edge.cex = 1.5,          
              legend.cex = .35,         
              scale.nodes = .6,        
              confidence = c('tp', 'pr', 'hg'), 
              pathways.color = pathways.color,  
              disconnected = F,        
              height.logic = .3,)
  }
  # export.graphml(LUAD.model, 
  #                file = "output/LUADgraphml.xml",
  #                pathways = pathway.list,  
  #                edge.cex = 1.5,          
  #                legend.cex = .35,         
  #                scale.nodes = .6,        
  #                confidence = c('tp', 'pr', 'hg'), 
  #                pathways.color = pathways.color,  
  #                disconnected = F,        
  #                height.logic = .3)
  # 
  # igraph <- read.graph(file = "output/LUADgraphml.xml", format = "graphml")
  # lisg <- as_adj_list(igraph, mode = "out")
  # lisge <- as_adj_edge_list(igraph, mode = "out")
  # matrix <- as_adjacency_matrix(igraph, sparse = FALSE, attr = "weight")
  # #matrix <- as.data.frame(matrix)
  # matrix <- matrix[rownames(matrix) %in% pathway.genes,]
  # matrix <- matrix[, colnames(matrix) %in% pathway.genes]
  # matrix <- matrix[, colnames(matrix) %in% pathway.genes]
  
  ## TODO add fishplot
  ## edit data to obtain something like this:
  
  #provide a list of timepoints to plot
  #You may need to add interpolated points to end up with the desired
  #visualization. Example here was actually sampled at days 0 and 150
  # timepoints=c(0,30,75,150) 
  #provide a matrix with the fraction of each population
  #present at each timepoint
  # frac.table = matrix(
  #   c(100, 45, 00, 00,
  #     02, 00, 00, 00,
  #     02, 00, 02, 01,
  #     98, 00, 95, 40),
  #   ncol=length(timepoints))
  
  #provide a vector listing each clone's parent
  #(0 indicates no parent)
  # parents = c(0,1,1,3)
  
  #create a fish object
  # fish = createFishObject(frac.table,parents,timepoints=timepoints)
  
  #calculate the layout of the drawing
  # fish = layoutClones(fish)
  
  #draw the plot, using the splining method (recommended)
  #and providing both timepoints to label and a plot title
  # fishPlot(fish,shape="spline",title.btm="Sample1",
  #          cex.title=0.5, vlines=c(0,150), 
  #          vlab=c("day 0","day 150"))
  ## excel with all data
  excel.file = paste("output/LUAD_statistics_", label, ".xlsx", sep='')
  
  excel.wbook = createWorkbook()
  
  sheet.luad.bic <- createSheet(wb = excel.wbook, 
                                sheetName="LUAD-bic")
  sheet.luad.aic <- createSheet(wb = excel.wbook, 
                                sheetName="LUAD-aic")
  
  addDataFrame(x = tabular(LUAD.model, 
                           'capri_bic'),
               sheet = sheet.luad.bic,
               showNA = T,
               characterNA = 'NA')
  addDataFrame(x = tabular(LUAD.model, 
                           'capri_aic'),
               sheet = sheet.luad.aic,
               showNA = T,
               characterNA = 'NA')
  
  saveWorkbook(excel.wbook, 
               excel.file)
  
}
