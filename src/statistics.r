# STATISTICAL ANALYSIS
#' Function to make statistical analysis of model LUAD
#'
#' @param LUAD.model TRONCO object model
#' @param label label to identify subtype
#' 
#' @return LUAD.model TRONCO object model
statistics <- function(LUAD.model, label, label.short) {
  ## non-parametric bootstrap
  LUAD.model <- tronco.bootstrap(LUAD.model,
                                 nboot = num_boot_iter,
                                 cores.ratio = .5)

  ## statistical bootstrap
  LUAD.model <- tronco.bootstrap(LUAD.model,
                                 type = "statistical",
                                 nboot = num_boot_iter,
                                 cores.ratio = .5)

  save(LUAD.model, file = paste0("input/model_boostrap", label.short, ".rda"))

  ## DAG of the model above
  ## Equal as model tronco.plot
  if (plot_verbose) {
    par(.pardefault)
    tronco.plot(
      LUAD.model,
      pathways = pathway.list,
      edge.cex = 1.5,
      legend.cex = .35,
      scale.nodes = .6,
      confidence = c('npb', 'sb'),
      pathways.color = pathways.color,
      disconnected = F,
      height.logic = .3,
      title = paste(label, "- first bootstrap")
    )

  }

  ## plot of bootstrap scores
  if (plot_verbose) {
    par(.pardefault)
    ## first non-parametric
    pheatmap(
      keysToNames(
        LUAD.model,
        as.confidence(LUAD.model,
                      conf = 'npb')$npb$capri_bic
      ) * 100,
      main = paste("non-parametric bootstrap scores for BIC model for", label),
      fontsize_row = 6,
      fontsize_col = 6,
      display_numbers = T,
      number_format = "%f"
    )
  }
  if (plot_verbose) {
    par(.pardefault)
    pheatmap(
      keysToNames(
        LUAD.model,
        as.confidence(LUAD.model,
                      conf = 'npb')$npb$capri_aic
      ) * 100,
      main = paste("non-parametric bootstrap scores for AIC model for", label),
      fontsize_row = 6,
      fontsize_col = 6,
      display_numbers = T,
      number_format = "%f"
    )
  }
  if (plot_verbose) {
    ## then parametric ones
    par(.pardefault)
    pheatmap(
      keysToNames(
        LUAD.model,
        as.confidence(LUAD.model,
                      conf = 'sb')$sb$capri_bic
      ) * 100,
      main = paste("statistical bootstrap scores for BIC model for", label),
      fontsize_row = 6,
      fontsize_col = 6,
      display_numbers = T,
      number_format = "%f"
    )
  }
  if (plot_verbose) {
    par(.pardefault)
    pheatmap(
      keysToNames(
        LUAD.model,
        as.confidence(LUAD.model,
                      conf = 'sb')$sb$capri_aic
      ) * 100,
      main = paste("statistical bootstrap scores for AIC model for", label),
      fontsize_row = 6,
      fontsize_col = 6,
      display_numbers = T,
      number_format = "%f"
    )
  }

  ## table with bootstrap scores
  boot_tab <- as.bootstrap.scores(LUAD.model)

  if (verbose) {
    print(boot_tab)
  }

  ## kfold
  ## k-fold cross validation, prediction error for each parent set X
  LUAD.model <- tronco.kfold.eloss(LUAD.model)
  kfold_eloss <- as.kfold.eloss(LUAD.model)


  if (verbose) {
    print("kfold loss")
    print(kfold_eloss)
  }

  ## plot for every fold
  if (plot_verbose) {
    vioplot(
      LUAD.model$kfold$capri_bic$eloss,
      LUAD.model$kfold$capri_aic$eloss,
    
      col = 'red',
      lty = 1,
      rectCol = "gray",
      colMed = 'black',
      names = c('BIC', 'AIC'),
      pchMed = 15,
      horizontal = T
    )
    
    legend(legend =c(paste("mean: ", round(kfold_eloss$Mean[1],2)),
                     paste("%log: ", round(kfold_eloss$`%-of-logLik`[1],2)),
                     paste("std dev: ", round(kfold_eloss$Stdev[1],2)),
                     paste("logLik: ", round(LUAD.model$model$capri_aic$logLik, 2)),
                     paste("score:  ", round(LUAD.model$model$capri_aic$score, 2))),
          x='topright',
          title = "AIC")
    
    legend(legend =c(paste("mean: ", round(kfold_eloss$Mean[2],2)),
                     paste("%log: ", round(kfold_eloss$`%-of-logLik`[2],2)),
                     paste("std dev: ", round(kfold_eloss$Stdev[2],2)),
                     paste("logLik: ", round(LUAD.model$model$capri_bic$logLik, 2)),
                     paste("score:  ", round(LUAD.model$model$capri_bic$score, 2))),
           x='bottomleft',
           title = "BIC")    
    
    title(main = paste('Entropy loss \n LUAD tumors -', label))
  }

  ## k-fold cross validation, prediction error for each parent set X
  LUAD.model <- tronco.kfold.prederr(LUAD.model)
  kfold_pred <- as.kfold.prederr(LUAD.model)

  if (verbose) {
    print("kfold prediction")
    print(kfold_pred)
  }

  ## k-fold cross validation, posterior classification error for each edge
  LUAD.model <- tronco.kfold.posterr(LUAD.model)
  kfold_post <- as.kfold.posterr(LUAD.model)

  if (verbose) {
    print("kfold posterior")
    print(kfold_post)
  }
  ## visualize a table with all edge statistics
  tab_bic <- tabular(LUAD.model, 'capri_bic')
  tab_aic <- tabular(LUAD.model, 'capri_aic')

  if (verbose) {
    print(label)
    print("table with all edge statistics using capri_bic")
    print(tab_bic)
    print("table with all edge statistics using capri_aic")
    print(tab_aic)
  }

  ## save model
  save(LUAD.model,
       file = paste("output/luadDefModel_", label.short, ".rda", sep = ''))
  
  # barplots 
  pred_aic <- data.frame(kfold_pred$capri_aic)
  pred_bic <- data.frame(kfold_pred$capri_bic)
  
  pred_sel <- pred_aic
  pred_sel[c('TYPE', 'NAME')] <- str_split_fixed(pred_sel$SELECTED, ' ', 2)
  pred_sel <- pred_sel[order(pred_sel$MEAN.PREDERR, decreasing = FALSE),]
  
  pred_sel$COLOR <- with(pred_sel, ifelse(TYPE == 'Deletion', '#8FBCBB',
                                          ifelse(TYPE=='Amplification', '#81A1C1',
                                                 ifelse(TYPE == 'Mutation','#D3AECC', '#D8DEE9' ))))
  
  barplot(pred_sel$MEAN.PREDERR,
          horiz = TRUE,
          names.arg = pred_sel$NAME, 
          las=2,
          col = pred_sel$COLOR,
          cex.names = 0.3,
          main = paste("Prediction error of AIC scores for -",label))
  legend("bottomright",                                    # Add legend to barplot
         legend = c("Amplification", "Deletion", "Mutation", "Pattern"),
         fill = c("#81A1C1", "#8FBCBB", "#D3AECC", "#D8DEE9"))
  
  pred_sel <- pred_bic
  pred_sel[c('TYPE', 'NAME')] <- str_split_fixed(pred_sel$SELECTED, ' ', 2)
  pred_sel <- pred_sel[order(pred_sel$MEAN.PREDERR, decreasing = FALSE),]
  
  pred_sel$COLOR <- with(pred_sel, ifelse(TYPE == 'Deletion', '#8FBCBB',
                                          ifelse(TYPE=='Amplification', '#81A1C1',
                                                 ifelse(TYPE == 'Mutation','#D3AECC', '#D8DEE9' ))))
  
  barplot(pred_sel$MEAN.PREDERR,
          horiz = TRUE,
          names.arg = pred_sel$NAME, 
          las=2,
          col = pred_sel$COLOR,
          cex.names = 0.3,
          main = paste("Prediction error of BIC scores for -",label))
  legend("bottomright",                                    # Add legend to barplot
         legend = c("Amplification", "Deletion", "Mutation", "Pattern"),
         fill = c("#81A1C1", "#8FBCBB", "#D3AECC", "#D8DEE9"))

  ## excel with all data
  excel.file <- paste0("output/LUAD_statistics_", label, ".xlsx")

  excel.wbook <- createWorkbook()

  sheet.luad.bic <- createSheet(wb = excel.wbook,
                                sheetName = "LUAD-bic")
  sheet.luad.aic <- createSheet(wb = excel.wbook,
                                sheetName = "LUAD-aic")

  addDataFrame(
    x = tabular(LUAD.model,
                'capri_bic'),
    sheet = sheet.luad.bic,
    showNA = T,
    characterNA = 'NA'
  )
  addDataFrame(
    x = tabular(LUAD.model,
                'capri_aic'),
    sheet = sheet.luad.aic,
    showNA = T,
    characterNA = 'NA'
  )

  saveWorkbook(excel.wbook,
               excel.file)

  return(LUAD.model)

}
