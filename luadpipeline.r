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
library(fishplot)
library(stringr)

## source file with the pipeline config
print("loading selected configurations")
source("src/conf.r")

## source file with some useful functions
print("loading useful functions")
## various functions
source("src/utils.r")
## function for model reconstruction
source("src/models.r")
## function for statistical analysis
source("src/statistics.r")

## source file with the input genes and the pathways config
print("loading selected genes")
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

## make analysis for every subtype
i <- 1
for (m in models) {

  ## model reconstruction
  troncomodel <- model(m, labels[i], labels.short[i])
  ## statistical analysis
  finalmodel <- statistics(troncomodel, labels[i], labels.short[i])

  if(plot_verbose)
    tronco.plot(finalmodel, 
                pathways = pathway.list, 
                fontsize = 15,
                edge.cex = 1.5,
                legend.cex = .7,
                scale.nodes = .6,
                confidence = c('tp', 'pr', 'hg', 'sb', 'npb'), 
                pathways.color = pathways.color,
                label.edge.size = 9,
                disconnected = F, 
                height.logic = .3,
                title = paste("Final Model -", labels[i]))

  i <- i + 1
}

print("END OF PiCnIc")

## pathway infomations
for (pw in pathway.list) {
  print("pathway")
  print(pw)
  my.pathways <- findPathwaysByText(paste(pw, collapse = ' '))
  my.hs.pathways <- my.pathways[my.pathways$species == 'Homo sapiens',]
  my.wpids <- my.hs.pathways$id
  print(my.hs.pathways[1]$name)
  print(my.hs.pathways[2]$name)
  print(my.hs.pathways[3]$name)
  pw.genes <- getXrefList(my.wpids[1], "H")
  ## open tabs in your favourite browser
  # browseURL(as.character(getPathwayInfo(my.wpids[1])[2]))
  print("-------------------------------------------------------------")
}

## manually cutared fishplots

## branching
par(.pardefault)
timepoints <- c(0, 30, 75, 150)
frac.table <- matrix(
  c(40, 10, 0,
    60, 30, 10,
    80, 70, 30,
    100, 90, 70),
  ncol = length(timepoints))
parents <- c(0, 1, 2)
fish <- createFishObject(frac.table,
                         parents,
                         timepoints = timepoints,
                         clone.annots = c("KEAP1", "RIT1", "ATM"),
                         clone.annots.angle = 30)
fish <- layoutClones(fish)
fish <- setCol(fish,
               col = c("#b48ead", "#a3be8c", "#8fbcbb"))
sample.times <- c(0, 150)
fishPlot(fish,
         shape = "spline",
         title.btm = "Sample1",
         cex.title = 1,
         vlab = c("day 0", "day 150"),
         bg.col = c("#ebcb8b", "#d08770", "#bf616a"))
par(.pardefault)
timepoints <- c(0, 30, 75, 150)
frac.table <- matrix(
  c(40, 25,
    60, 50,
    80, 75,
    100, 90),
  ncol = length(timepoints))
parents <- c(0, 1)
fish <- createFishObject(frac.table,
                         parents,
                         timepoints = timepoints,
                         clone.annots = c("KEAP1", "RIT1"),
                         clone.annots.angle = 30)
fish <- layoutClones(fish)
fish <- setCol(fish,
               col = c("#b48ead", "#5e81ac"))

sample.times <- c(0, 150)
fishPlot(fish,
         shape = "spline",
         title.btm = "Sample2",
         cex.title = 1,
         vlab = c("day 0", "day 150"),
         bg.col = c("#ebcb8b", "#d08770", "#bf616a"))

# confluence
par(.pardefault)
timepoints <- c(0, 30, 75, 150)
frac.table <- matrix(
  c(30, 20, 0,
    60, 40, 15,
    80, 80, 50,
    100, 95, 90),
  ncol = length(timepoints))
parents <- c(0, 1, 2)
fish <- createFishObject(frac.table,
                         parents,
                         timepoints = timepoints,
                         clone.annots = c("TP53", "KEAP1", "ARID2"),
                         clone.annots.angle = 30)
fish <- layoutClones(fish)
fish <- setCol(fish,
               col = c("#8fbcbb", "#a3be8c", "#5e81ac"))
sample.times <- c(0, 150)
fishPlot(fish,
         shape = "spline",
         title.btm = "Sample3",
         cex.title = 1,
         vlab = c("day 0", "day 150"),
         bg.col = c("#ebcb8b", "#d08770", "#bf616a"))
par(.pardefault)
timepoints <- c(0, 30, 75, 150)
frac.table <- matrix(
  c(40, 25, 25,
    60, 50, 50,
    80, 75, 75,
    100, 90, 90),
  ncol = length(timepoints))
parents <- c(0, 1)
fish <- createFishObject(frac.table,
                         parents,
                         timepoints = timepoints,
                         clone.annots = c("KEAP1", "ARID2"),
                         clone.annots.angle = 30)
fish <- layoutClones(fish)
fish <- setCol(fish,
               col = c("#b48ead", "#5e81ac"))

sample.times <- c(0, 150)
fishPlot(fish,
         shape = "spline",
         title.btm = "Sample4",
         cex.title = 1,
         vlab = c("day 0", "day 150"),
         bg.col = c("#ebcb8b", "#d08770", "#bf616a"))
timepoints <- c(0, 30, 75, 150)
frac.table <- matrix(
  c(40, 25,
    60, 50,
    80, 75,
    100, 90),
  ncol = length(timepoints))
parents <- c(0, 1, 2)
fish <- createFishObject(frac.table,
                         parents,
                         timepoints = timepoints,
                         clone.annots = c("TP53", "KEAP1", "ARID2"),
                         clone.annots.angle = 30)
fish <- layoutClones(fish)
fish <- setCol(fish,
               col = c("#8fbcbb", "#a3be8c", "#5e81ac"))
sample.times <- c(0, 150)
fishPlot(fish,
         shape = "spline",
         title.btm = "Sample3",
         cex.title = 1,
         vlab = c("day 0", "day 150"),
         bg.col = c("#ebcb8b", "#d08770", "#bf616a"))
par(.pardefault)
