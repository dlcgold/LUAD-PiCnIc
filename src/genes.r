# gene selection

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4231481/
## In TCGA marker paper we have 40 gens drivers divided in 8 pathways,
## detected using MutSiGCV tool
## We ignore every other genes even if very frequent
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

## extracted from maftools
## understand if they are usefull
## HFREQ <- c("TTN", "MUC16", "RYR2", "CSMD3", "LRP1B")

## create pathways for various plots using TRONCO library
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

## colors for pathways in the plots 
## (dark colors because they will be used for writing names)
alteration.color <- 'dimgray'
pathways.color <- c('darkslategray',
                   'darkblue', 
                   'darkgreen',
                   'darkmagenta', 
                   'firebrick4',
                   'dodgerblue4',
                   'darkorchid1',
                   'darkorange4')