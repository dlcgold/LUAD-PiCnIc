## files
file.mutex <- "input/LUAD_mutex.txt"
file.clinical <- "input/LUAD_clinical.txt"
file_drivers <- "input/gene_drivers.xlsx"

## bool for print and plot
verbose <- TRUE
plot_verbose <- TRUE

## bool for not distinguish mutation
all_mut <- FALSE
if(all_mut){
  mut <- 'Nonsense_Mutation'
}else{
  mut <- 'Mutation'
}

## bool for data reloading ()
maf_reload <- TRUE
clinic_reload <- maf_reload
gistic_reload <- maf_reload


## min frequency
min_freq <- 0.03

## boot iteration, should be around 100
num_boot_iter <- 10
