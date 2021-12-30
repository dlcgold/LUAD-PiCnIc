## sourcing every Library used
## TODO remove useless ones
library(TRONCO)
library(RTCGAToolbox)
library(TCGAbiolinks)
library(readr)
library(dtutils)
library(devtools)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(vioplot)
library(openxlsx)
library(dplyr)
library(plyr)
library(maftools)
library(RColorBrewer)
library(GenVisR)
library(fishplot)
library(igraph)
library("xlsx")
library(rWikiPathways)

### PIPELINE CONFIGURATION
## files (mutex, clinical, genes drivers from IMCDriver)
file.mutex <- "input/LUAD_mutex.txt"
file.clinical <- "input/LUAD_clinical.txt"
## not used in final version
file_drivers <- "input/gene_drivers.xlsx"

## flags for textual and visual output
verbose <- TRUE
plot_verbose <- TRUE
histological_verbose <- FALSE

## flag for not distinguish mutation
## and eventually know how to distinguish mutations
all_mut <- FALSE
if (all_mut) {
  mut <- 'Nonsense_Mutation'
} else{
  mut <- 'Mutation'
}

## flag for data reload (maf, clinical and gistic)
maf_reload <- TRUE
clinic_reload <- maf_reload
gistic_reload <- maf_reload

## min frequency for events
min_freq <- 0.03

## bootstrap iteration, should be around 100
num_boot_iter <- 5

## workaround for plots
.pardefault <- par()
par(.pardefault)
## setwd('~/DCB-project/')
### END CONF

#source file with some useful functions
print("loading useful functions")
## various functions
source("src/utils.r")

## function for model reconstruction
source("src/models.r")

## function for statistical analysis
source("src/statistics.r")

## source file with the pipeline config
print("loading selected configurations")
source("src/conf.r")

print("loading selected genes")
## source file with the input genes and the pathways config
source("src/genes.r")

## source file with the data loading
## (MUTEX, MAF, CLINICAL, GISTIC and MAF-GISTIC intersection)
print("loadind data")
source("src/load_data.r")

## source file with the subtyping analysis
print("subtyping analysis")
source("src/subtyping.r")

## source file with the first group exclusivity analysis
print("group exclusivity analysis")
source("src/group_exclusivity.r")


# models list for analysis (the first is dataset without subtype selection)
models <- list(LUAD,
               LUAD.TRU,
               LUAD.PI,
               LUAD.PP)
# LUAD.acinar,
# LUAD.nonmucinous,
# LUAD.papillary,
# LUAD.mucinous,

## labels for every subtype (the first is dataset without subtype selection)
labels <- c(
  'all subtypes',
  'terminal respiratory unit (TRU, branchoid)',
  'proximal inflammatory (PI, squamoid)',
  'proximal proliferative (PP, magnoid)')
  # 'acinar',
  # 'nonmucinous',
  # 'papillary',
  # 'mucinous',
## mucinous subgroup is too small!

## model reconstruction parametes
## TODO they are random at the moment
## gene.hypotheses <- c('KRAS', 'BRAF', 'ATM', 'STK11')



## make analysis for every subtype
i <- 1
for (m in models) {
  ## MODELS CONF

  ## model reconstruction
  troncomodel <- model(m, labels[i])
  ## statistical analysis
  statistics(troncomodel, labels[i])
  i <- i + 1
}

print("END OF PiCnIc")


for (pw in pathway.list) {
  print("pathway")
  print(pw)
  my.pathways <-
    findPathwaysByText(paste(pw, collapse = ' '))
  my.hs.pathways <-
    my.pathways[my.pathways$species == 'Homo sapiens',]
  my.wpids <- my.hs.pathways$id
  # pw.title <- my.hs.pathways[1]$name
  print(my.hs.pathways[1]$name)
  print(my.hs.pathways[2]$name)
  print(my.hs.pathways[3]$name)
  pw.genes <- getXrefList(my.wpids[1], "H")
  # browseURL(as.character(getPathwayInfo(my.wpids[1])[2]))
  print("-------------------------------------------------------------")
}
