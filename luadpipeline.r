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

print("subtyping analysis")
## source file with the subtyping analysis
source("src/subtyping.r")

print("group exclusivity analysis")
## source file with the first group exclusivity analysis
source("src/group_exclusivity.r")


# CONF For Models
models <- list(LUAD,
               LUAD.acinar, 
               LUAD.nonmucinous, 
               LUAD.papillary, 
               LUAD.mucinous)
labels <- c('',
            'acinar',
            'nonmucinous',
            'papillary',
            'mucinous')
excluded <- c('mucinous')


# model selection
gene.hypotheses <- c('KRAS', 'BRAF', 'ATM', 'STK11')
gene.sel <- P53
genes.compare <- c('TP53', 'ATM')
genes.to <- c('KRAS', mut)

## mucinous subgroup is too small!
i <- 1
for(m in models){
  if (labels[i]!=excluded){
    ##print(nsamples(m))
    troncomodel <- model(m, gene.hypotheses, gene.sel, genes.compare, 
                         genes.to, labels[i])
    statistics(troncomodel, labels[i])
  }
  i <- i + 1
}
