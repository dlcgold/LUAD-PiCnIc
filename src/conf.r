## files (mutex, clinical, genes drivers from IMCDriver)
file.mutex <- "input/LUAD_mutex.txt"
file.clinical <- "input/LUAD_clinical.txt"
## not used in final version
file_drivers <- "input/gene_drivers.xlsx"

## bool for print some textual results and plot
verbose <- TRUE
plot_verbose <- TRUE

## bool for not distinguish mutation 
## and eventually know how to distinguish mutations
all_mut <- FALSE
if(all_mut){
  mut <- 'Nonsense_Mutation'
}else{
  mut <- 'Mutation'
}

## bool for data reloading (maf, clinical and gistic)
maf_reload <- TRUE
clinic_reload <- maf_reload
gistic_reload <- maf_reload


## min frequency for events
min_freq <- 0.03

## bootstrap iteration, should be around 100
num_boot_iter <- 5