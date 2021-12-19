# gene selection

## use pathway as in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4231481/
## complete genes
P53 <- c("TP53", "ATM", "MDM2")
MAPK <- c("KRAS", "NRAS", "HRAS", "RIT1", "NF1", "BRAF", "MAP2K1",
         "EGFR", "ERBB2", "MET", "ALK", "RET", "ROS1")
MTOR <- c("PTEN", "PIK3CA", "PIK3R1", "STK11", "AKT1", "AMPK",
          "TSC1", "TSC2", "MTOR")
OXI <- c("KEAP1", "CUL3", "NFE2L2")
PROG <- c("CDKN2A", "CCND1", "CDK4", "CCNE1", "RB1")
REMO <- c("ARID1A", "ARID1B", "ARID2", "SMARCA4")
HIME <- c("SETD2")
RNASPL <- c("RBM10", "U2AF1")


## genes > 1%
## TODO use correct pathway names (search in papers)
# - P53: proliferation and cell cycle progression
# - RTK: RTK signalling, proliferation, cell survival, translation
# - MTOR: mTOR signalling, proliferation, cell survival, translation
# - OXI: oxidative stress response
# - PROG: cell cycle progression
# - REMO: nucleosome remodelling
# - HIME: histone methylation
# - RNASPL: RNA splicing/processing
# P53 <- c("TP53", "ATM", "MDM2")
# MAPK <- c("KRAS", "RIT1", "NF1", "BRAF", "EGFR", "ERBB2", "MET", "MAP2K1")
# MTOR <- c("PTEN", "PIK3CA", "STK11", "TSC1", "TSC2", "AKT1", "AMPK", "MTOR")
# OXI <- c("KEAP1", "NFE2L2")
# PROG <- c("CDKN2A", "CCND1", "CDK4", "CCNE1", "RB1")
# REMO <- c("ARID1A", "ARID1B", "ARID2", "SMARCA4")
# HIME <- c("SETD2")
# RNASPL <- c("RBM10", "U2AF1")

## extracted from maftools
## understand if they are usefull
## HFREQ <- c("TTN", "MUC16", "RYR2", "CSMD3", "LRP1B")

## create pathways
pathway.genes <- c(P53, MAPK, MTOR, OXI, PROG, REMO, HIME, RNASPL)
pathway.genes <- unique(pathway.genes)
pathway.names <- c("P53", "MAPK", "MTOR", "OXI", "PROG", "REMO", 
                   "HIME", "RNASPL")
pathway.list <- list(P53 = P53,
                     MAPK = MAPK,
                     MTOR = MTOR, 
                     OXI = OXI, 
                     PROG = PROG, 
                     REMO = REMO, 
                     HIME = HIME, 
                     RNASPL = RNASPL)

## colors for pathways
alteration.color = 'dimgray'
pathways.color = c('firebrick1', 
                   'darkblue', 
                   'darkgreen',
                   'darkmagenta', 
                   'darkorange',
                   'dodgerblue4',
                   'darkorchid1',
                   'darksalmon')