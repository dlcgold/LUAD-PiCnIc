
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
              height.logic = .3, 
              create.new.dev = TRUE)
              
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
  pheatmap(keysToNames(LUAD.model,
                       as.confidence(LUAD.model,
                                     conf = 'npb')$npb$capri_aic) * 100,
           main = paste("non-parametric bootstrap scores for AIC model for",label),
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%f"
  )
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
  title(main = 'Entropy loss \n LUAD tumors')
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
tronco.plot(LUAD.model, 
            pathways = pathway.list,  
            edge.cex = 1.5,          
            legend.cex = .35,         
            scale.nodes = .6,        
            confidence = c('tp', 'pr', 'hg'), 
            pathways.color = pathways.color,  
            disconnected = F,        
            height.logic = .3,
            create.new.dev = TRUE)

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
