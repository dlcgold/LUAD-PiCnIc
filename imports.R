## conda create -n luad-proj -y python=3.9 r-essentials
## conda activate luad-proj
## conda install -y -c conda-forge r-gert r-rjava

## run inside R
install.packages("devtools")
install.packages("BiocManager")
library("devtools")

BiocManager::install(c("TRONCO",
                       "RTCGAToolbox",
                       "TCGAbiolinks",
                       "maftools",
                       "GenVisR"))
install.packages(c("readr",
                   "pheatmap",
                   "ggplot2",
                   "gridExtra",
                   "vioplot",
                   "openxlsx",
                   "dplyr",
                   "plyr",
                   "xlsx"))
install_github("chrisamiller/fishplot")
install_github("and3k/dtutils")

