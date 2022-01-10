# LUAD-PiCnIc
## Lung Adenocarcinoma (_LUAD_) model progression

A cancer progression model from cross-sectional lung adenocarcinoma data, using
the PiCnIc pipeline.

## Prerequisites
All the prerequisites could be satisfied running `import.r` script

### CRAN
- readr
- pheatmap
- ggplot2
- gridExtra
- vioplot
- openxlsx
- dplyr
- plyr
- xlsx

### Bioconductor
- TRONCO
- RTCGAToolbox
- TCGAbiolinks
- maftools
- GenVisR

### GitHub

- chrisamiller/fishplot
- and3k/dtutils

## Execution
In order to execute this script internet connection is required.
First of all, in _R_ console, set the working directory to this one:
```R 
 setwd("~/DCB-project")
 ```
then simply source `luadpipe.r`, making sure to adjust your configuration in
`src/conf.R` before.

With default bootstrap iterations (100) the script will take several tens of 
minutes to finish the computation. 

In `/output` directory will be saved the final statistical results and models,
while in `/input` directory, in addition to some input files (such as those
for mutex), partial models obtained during computation will also be saved.