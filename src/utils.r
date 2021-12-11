loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
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

reduce_samples <- function(LUAD, n_samples){
  uni_samples <- unique(LUAD$Tumor_Sample_Barcode)
  if(length(uni_samples) < n_samples){
    n_samples <- length(uni_samples)
  }
  
  samples <- uni_samples[1:n_samples]
  LUADreduced <- subset(LUAD, Tumor_Sample_Barcode %in% samples)
  return(LUADreduced)
}

select = function(x, min.freq, forced.genes) {
  x.sel = as.alterations(x)
  x.sel = events.selection(x.sel, 
                           filter.freq = min.freq, 
                           filter.in.names = forced.genes)
  

  x = events.selection(x, 
                       filter.in.names = as.genes(x.sel))
  return(x)
}

