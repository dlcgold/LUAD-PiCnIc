
# STATISTICS

## non-parametric bootstrap
if(boot_reload){
  LUAD.model <- tronco.bootstrap(LUAD.model,
                                 nboot = num_boot_iter,
                                 cores.ratio = .5)
  
  ## statistical bootstrap
  LUAD.model <- tronco.bootstrap(LUAD.model,
                                 type = "statistical",
                                 nboot = num_boot_iter,
                                 cores.ratio = .5)
  
  save(LUAD.model, file = "input/luadDefBoot.rda")
}else{
  LUAD.model <- loadRData("input/luadDefBoot.rda")
}

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
              height.logic = .3,       
              #file = "output/model_hypo_boot.pdf"
  )
}

## plot of bootstrap scores
## TODO sometimes not work
if(plot_verbose && FALSE){
  ## first non-parametric
  pheatmap(keysToNames(LUAD.model,
                       as.confidence(LUAD.model,
                                     conf = 'npb')$npb$capri_bic) * 100,
           main = "non-parametric bootstrap scores for BIC model",
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%f"
  )
  pheatmap(keysToNames(LUAD.model,
                       as.confidence(LUAD.model,
                                     conf = 'npb')$npb$capri_aic) * 100,
           main = "non-parametric bootstrap scores for AIC model",
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%f"
  )
  ## then parametric ones
  pheatmap(keysToNames(LUAD.model,
                       as.confidence(LUAD.model,
                                     conf = 'sb')$sb$capri_bic) * 100,
           main = "non-parametric bootstrap scores for BIC model",
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%f"
  )
  pheatmap(keysToNames(LUAD.model,
                       as.confidence(LUAD.model,
                                     conf = 'sb')$sb$capri_aic) * 100,
           main = "non-parametric bootstrap scores for AIC model",
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
if(kfold_reload){
  ## k-fold cross validation, prediction error for each parent set X
  LUAD.model <- tronco.kfold.eloss(LUAD.model)
  kfold_eloss <- as.kfold.eloss(LUAD.model)
  
  ## plot for every fold
  ## TODO make it work
  if(plot_verbose && FALSE){
    vioplot(LUAD.model$kfold$capri_bic$eloss,
            LUAD.model$kfold$capri_aic$eloss,
            col = 'red',
            lty = 1, rectCol="gray",
            colMed = 'black',
            names = c('BIC', 'AIC'), 
            pchMed = 15, 
            horizontal = T)
    title(main = 'Entropy loss \n LUAD tumors')
  }
  
  ## k-fold cross validation, prediction error for each parent set X
  LUAD.model <- tronco.kfold.prederr(LUAD.model)
  kfold_pred <- as.kfold.prederr(LUAD.model)
  
  
  ## k-fold cross validation, posterior classification error for each edge
  LUAD.model <- tronco.kfold.posterr(LUAD.model)
  kfold_post <- as.kfold.posterr(LUAD.model)
  save(LUAD.model, 
       file = "input/luadDefKfold.rda")
  
}else{
  LUAD.model <- loadRData("input/luadDefKfold.rda")
}

## visualize a table with all edge statistics
tab_bic <- tabular(LUAD.model, 'capri_bic')
tab_aic <- tabular(LUAD.model, 'capri_aic')

if(verbose){
  print("table with all edge statistics using capri_bic")
  print(tab_bic)
  print("table with all edge statistics using capri_aic")
  print(tab_aic)
}

## save model
save(LUAD.model, 
     file = "input/luadDefModel.rda")

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
              height.logic = .3,       
              #file = "output/model_hypo_boot.pdf"
  )
}

## excel with all data
excel.file = "output/LUAD_statistics.xlsx"

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
