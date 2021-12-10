## files
file.mutex <- "input/LUAD_mutex.txt";
file.clinical <- "input/LUAD_clinical.txt";
file_drivers <- "input/gene_drivers.xlsx";

## bool for print and plot
verbose <- TRUE;
plot_verbose <- TRUE;

## bool for data reloading ()
maf_reload <- TRUE;
clinic_reload <- maf_reload;
hypo_reload <- maf_reload;
gistic_reload <- maf_reload;
intersect_reload <- maf_reload;
boot_reload <- maf_reload;
## for now FALSE do to error
# gistic_reload <- FALSE;
# intersect_reload <- FALSE;

## min frequency
min_freq <- 0.02;

## boot iteration, should be around 100
num_boot_iter <- 2
