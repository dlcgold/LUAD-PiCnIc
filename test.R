
par(.pardefault)
troncomodel <- loadRData(file='input/saved_modelPI.rda')
tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.2,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD PI subtype - capri"
)


troncomodel <- loadRData(file='input/model_boostrapTRU.rda')
tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg'),
  pathways.color = pathways.color,
  disconnected = T,
  height.logic = .3,
  title = "LUAD TRU subtype - capri"
)



troncomodel <- loadRData(file='input/model_boostrapPI.rda')
tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD PI subtype - capri"
)

troncomodel <- loadRData(file='input/model_boostrapPP.rda')
tronco.plot(
  troncomodel,
  pathways = pathway.list,
  edge.cex = 1.5,
  legend.cex = .4,
  scale.nodes = .6,
  confidence = c('tp', 'pr', 'hg'),
  pathways.color = pathways.color,
  disconnected = F,
  height.logic = .3,
  title = "LUAD PP subtype - capri"
)



install_github('hdng/clonevol')
install.packages('gridBase')
install.packages('gridExtra')
install.packages('ggplot2')
install.packages('igraph')
install.packages('packcircles')
install_github('hdng/trees')

troncomodel <- loadRData(file='input/model_boostrapTRU.rda')
statistics(troncomodel, "TRU subtype", "TRU")



install_github("phillipnicol/OncoBN")
pkgbuild::check_build_tools(debug = TRUE)

devtools::load_all()

options(buildtools.check = function(action) TRUE )








### -------------- VIOPLOTS




troncomodel <- loadRData("input/model_boostrapPP.rda")
label <- label.pp
label.short <- label.pp
LUAD.model <- troncomodel

# LUAD.model <- loadRData(paste0("input/model_boostrap", label.short, ".rda"))
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
## TODO sometimes not work
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
## TODO make it work
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

pred_aic <- data.frame(kfold_pred$capri_aic)
pred_bic <- data.frame(kfold_pred$capri_bic)

pred_sel <- pred_aic
library(stringr)

#split 'player' column using '_' as the separator
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
        main = "PREDICTION ERROR FOR AIC SCORES",)
legend("bottomright",                                    # Add legend to barplot
       legend = c("Amplification", "Deletion", "Mutation", "Pattern"),
       fill = c("#81A1C1", "#8FBCBB", "#D3AECC", "#D8DEE9"))





