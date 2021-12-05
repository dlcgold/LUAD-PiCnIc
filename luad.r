library(TRONCO)

source("src/input.r")
source("src/stages.r")

LUAD <- input_genes(interactive = FALSE);
oncoprint(LUAD);
LUAD_alt <- as.alterations(LUAD);
LUAD_genes <- ngenes(LUAD);
LUAD_events <- nevents(LUAD);
LUAD_samples <- nsamples(LUAD);
LUAD_types <- ntypes(aCML);
LUAD_patterns <- npatterns(aCML);
# TODO vedere se si possono clusterizzare alcune mutazioni
stagesStudy(LUAD, LUAD_samples, 3);