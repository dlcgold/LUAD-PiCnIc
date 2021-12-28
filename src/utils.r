#' Function to load R object, in .rda format, into current workspace 
#'
#' @param fileName filename to load
#'
#' @return an R object
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#' Function to merge together all tables regarding edge statistics
#'
#' @param obj the TRONCO model object
#' @param M type of model
#'
#' @return the table with all the statistics
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


#' Function to reduce dataset to a selected number of samples
#'
#' @param dataset a TRONCO dataset
#' @param n_samples the number of samples desired
#'
#' @return the dataset reduced, as a TRONCO dataset
reduce_samples <- function(dataset, n_samples){
  uni_samples <- unique(dataset$Tumor_Sample_Barcode)
  if(length(uni_samples) < n_samples){
    n_samples <- length(uni_samples)
  }
  samples <- uni_samples[1:n_samples]
  datasetreduced <- subset(dataset, Tumor_Sample_Barcode %in% samples)
  return(datasetreduced)
}

#' Function to select a set of events of a TRONCO object based on conditions
#'
#' @param x a TRONCO object
#' @param min.freq minumun frequency for events
#' @param forced.genes a genes selection
#'
#' @return TRONCO onject with selected events
select <- function(x, min.freq, forced.genes) {
  ## Collapse multiple events per gene in one unique event 
  x.sel <- as.alterations(x)
  ## Events selection based on frequency and genes
  x.sel <- events.selection(x.sel, 
                            filter.freq = min.freq, 
                            filter.in.names = forced.genes)
  
  ## Subset input
  x <- events.selection(x, 
                        filter.in.names = as.genes(x.sel))
  return(x)
}