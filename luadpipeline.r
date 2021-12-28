library(TRONCO)
library(RTCGAToolbox)
library(TCGAbiolinks)
library(readr)
library(dtutils)
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

# setwd('~/DCB-project/')

print("loading useful functions")
## source file with some useful functions
source("src/utils.r")

print("loading model function for hypotheses analysis")
# source file with the model reconstruction based on hypotheses
source("src/models.r")

print("loading selected configurations")
## source file with the pipeline config
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
    print(nsamples(m))
    model(m, gene.hypotheses, gene.sel, genes.compare, genes.to, labels[i])
    }
  i <- i + 1
}
