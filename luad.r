library(TRONCO)

source("src/input.r")
source("src/stages.r")
source("src/hypo_analysis.r")

LUAD <- input_genes(interactive = TRUE);
oncoprint(LUAD);
# tronco.plot(
#   tronco.capri(LUAD), 
#   scale.nodes = .6, 
#   legend.cex = .6, 
#   legend.pos = 'top', 
#   confidence = c('tp', 'pr', 'hg'))
#LUAD_alt <- as.alterations(LUAD);
LUAD_genes <- ngenes(LUAD);
LUAD_events <- nevents(LUAD);
LUAD_samples <- nsamples(LUAD);
LUAD_types <- ntypes(LUAD);
LUAD_patterns <- npatterns(LUAD);
# TODO vedere se si possono clusterizzare alcune mutazioni
# TODO stages stagesStudy(LUAD, LUAD_samples, 3);
hypo_analysis(LUAD);