library(TRONCO)
library(xlsx)
source("src/input.r")
source("src/stages.r")
source("src/hypo_analysis.r")

LUAD <- input_genes_gdc(interactive = FALSE);
# LUADGistic <- input_gistic(interactive = TRUE);
oncoprint(LUAD);
# tronco.plot(
#   tronco.capri(LUAD),
#   scale.nodes = .9,
#   legend.cex = .9,
#   legend.pos = 'top',
#   confidence = c('tp', 'pr', 'hg'))
LUAD_alt <- as.alterations(LUAD);
LUAD_genes <- ngenes(LUAD);
LUAD_events <- nevents(LUAD);
LUAD_samples <- nsamples(LUAD);
LUAD_types <- ntypes(LUAD);
LUAD_patterns <- npatterns(LUAD);
# TODO vedere se si possono clusterizzare alcune mutazioni
# TODO mettere vere ipotesi
LUAD.hypo <- hypo_analysis(LUAD);
# TODO stages 
stagesStudy(LUAD, LUAD_samples, 30);

