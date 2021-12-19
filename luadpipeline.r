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
library("xlsx")

print("loading useful functions")
## source file with some useful functions
source("src/utils.r")

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

models <- list(LUAD,
               LUAD.acinar, 
               LUAD.nonmucinous, 
               LUAD.papillary, 
               LUAD.mucinous)

## TODO for every model make models.r and statistics.r
for(model in models){
  LUADtest <- model
  print(nsamples(LUADtest))
}

print("model with hypotheses analysis")
## source file with the model reconstruction based on hypotheses
source("src/models.r")

## source file with the statistical analysis
source("src/statistics.r")

