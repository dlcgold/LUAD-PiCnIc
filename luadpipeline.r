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

## workaround for plots
.pardefault <- par()

# setwd('~/DCB-project/')

## source file with some useful functions
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
               LUAD.acinar, 
               LUAD.nonmucinous, 
               LUAD.papillary, 
               LUAD.mucinous)

## labels for every subtype (the first is dataset without subtype selection)
labels <- c('all',
            'acinar',
            'nonmucinous',
            'papillary',
            'mucinous')
## mucinous subgroup is too small!
## excluded <- c('mucinous')
excluded <- 'mucinous'

## model reconstruction parametes
## TODO they are random at the moment
## gene.hypotheses <- c('KRAS', 'BRAF', 'ATM', 'STK11')
gene.sel <- P53
genes.compare <- c('TP53', 'ATM')
genes.to <- c('KRAS', mut)

## make analysis for every subtype
i <- 1
for(m in models){
  if (labels[i]!=excluded){
    ## model reconstruction
    troncomodel <- model(m, gene.sel, genes.compare,
                         genes.to, labels[i])
    ## statistical analysis
    statistics(troncomodel, labels[i])
  }
  i <- i + 1
}

print("END")
