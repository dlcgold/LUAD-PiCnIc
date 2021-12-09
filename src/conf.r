## files
file.mutex <- "input/LUAD_mutex.txt";
file.clinical <- "input/LUAD_clinical.txt";
file_drivers <- "input/gene_drivers.xlsx";

## bool for print and plot
verbose = TRUE;
plot_verbose = TRUE;
gistic_check = FALSE;
hypo_reload = TRUE;

## color
pathways.color = c('firebrick1', 
                   'darkblue', 
                   'darkgreen',
                   'darkmagenta', 
                   'darkorange');
## min frequency
min_freq <- 0.05;

## boot iteration, should be around 100
num_boot_iter <- 2;
