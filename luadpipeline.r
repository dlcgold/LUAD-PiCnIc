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
library(maftools)
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

print("model with hypotheses analysis")
## source file with the model reconstruction based on hypotheses
source("src/model.r")

## source file with the statistical analysis
source("src/statistics.r")

