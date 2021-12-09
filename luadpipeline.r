library(TRONCO)
library(RTCGAToolbox)
library(TCGAbiolinks)
library(readr)
library(dtutils)
library(ggplot2)
library(gridExtra)
library("xlsx")
source("src/utils.r")

# LOAD DATA

## files
file.mutex <- "input/LUAD_mutex.txt";
file.clinical <- "input/LUAD_clinical.txt";
file_drivers <- "input/gene_drivers.xlsx";

## bool for print and plot
verbose = TRUE;
plot_verbose = TRUE;

## color
pathways.color = c('firebrick1', 
                   'darkblue', 
                   'darkgreen',
                   'darkmagenta', 
                   'darkorange');

## gene selection
## TODO hardcode this part with correct pathway

LUAD.mutex <- import.mutex.groups(file.mutex);
genes <- read.xlsx(file_drivers, 
                   sheetIndex = 1, 
                   header = TRUE);
genes <- genes$LUAD;
for (group in LUAD.mutex){
  genes <- c(genes, group[[1]]);
}
genes <- append(genes, "ERBB2");
genes <- append(genes, "ALK");
genes <- append(genes, "NRAS");
genes <- append(genes, "PIK3CA");
genes <- append(genes, "CTNNB1");
genes <- append(genes, "STK11");
genes <- append(genes, "CDKN2A");
genes <- append(genes, "APC");
genes <- append(genes, "RB1");
genes <- append(genes, "NRAS");
genes <- append(genes, "WRN");
genes <- append(genes, "RHOC");
genes <- unique(genes);

## load data if not exist a rda
if(!file.exists("input/luadDef.rda")){
  
  ## download MAF file from TCGA
  LUAD.maf <- GDCquery_Maf(tumor = "LUAD", 
                           pipelines = "mutect2");
  LUAD.maf <- as.data.frame(LUADmut);
  
  ## select only desired genes
  ## LUADclear <- subset(LUADmut, 
  ##                     Hugo_Symbol %in% genes);
  
  ## create TRONCO object
  LUAD = import.MAF(LUAD.maf,
                    is.TCGA = TRUE,
                    merge.mutation.types = FALSE,
                    filter.fun = function(x) {
                      return(x['Hugo_Symbol'] %in% genes)
                    });
  LUAD <- annotate.description(LUAD, 
                               "Lung cancer data from GDC portal");
  save(LUAD, file = "input/luadDef.rda");
}else{
  LUAD <- loadRData("input/luadDef.rda");
}

## a first and brutal oncoprint
if(plot_verbose){
  oncoprint(LUAD);
}

## print some informations about genes
if(verbose){
  print(paste("number of LUAD genes:",
              ngenes(LUAD)));    
  print("LUAD genes:");
  print(as.genes(LUAD))                          
  print("LUAD genes for type");
  for (type in as.types(LUAD)) {
    print(paste("gene for type:", 
                type));
    print(as.genes(LUAD,
                   types = type)); 
  }
}

## print some informations about events
if(verbose){
  print(paste("number of LUAD events:",
              nevents(LUAD)));    
  print("LUAD events:");
  print(as.events(LUAD));                       
  ## TODO add other similar events that can be queried together
  print("similar events:");
  print(as.events(LUAD, genes = c('KRAS', 
                                  'TP53')));
  print("LUAD alteration types compacted:");
  print(as.events(LUAD),
        keysToNames = TRUE);
}

## print some informations about types of alteration
if(verbose){
  print(paste("number of LUAD alteration types:", 
              ntypes(LUAD)));    
  print("LUAD alteration types:");
  print(as.types(LUAD)); 
  ## TODO add other similar events that can be queried together
  print("alteration types for similar events:");
  print(as.types(LUAD, genes = c('KRAS', '
                                 TP53')));
}


## print samples information
if(verbose){
  print(paste("number of LUAD samples:", 
              nsamples(LUAD)));    
  print("LUAD samples:");
  print(as.samples(LUAD)); 
  ## TODO add more of these
  print(which.samples(LUAD, 
                      gene = 'KRAS', 
                      type = 'Missense_Mutation'));
}

## clinical data

if(!file.exists("input/luadClinical.rda")){
  data <- getFirehoseData("LUAD");
  clinical <- getData(data, "clinical");
  
  df <- as.data.frame(clinical);
  for (i in 1:length(rownames(df))) {
    row.names(df)[i] <- toupper(gsub("\\.", 
                                     "-", 
                                     row.names(df)[i]));
    df$pathologic_stage[i] <- gsub("STAGE", 
                                   "Stage", 
                                   toupper(df$pathologic_stage[i]));
  }
  
  dtutils::write_tsv(df, 
                     file = file.clinical,
                     row_names = TRUE);
  
  clinical.data <- TCGA.map.clinical.data(file = file.clinical,
                                          column.samples = 'rn',
                                          column.map = 'pathologic_stage');
  save(clinical.data, file = "input/luadClinical.rda");
}else{
  clinical.data <- loadRData("input/luadClinical.rda");
}

if(verbose){
  print(head(clinical.data));
}

## match samples and stages
LUAD <- annotate.stages(LUAD, 
                        clinical.data, 
                        match.TCGA.patients = TRUE);

## clear multiple stages
LUAD <- TCGA.remove.multiple.samples(LUAD);
LUAD <- TCGA.shorten.barcodes(LUAD);
LUAD <- annotate.stages(LUAD, 
                        clinical.data);

## another brutal oncoprint
if(plot_verbose){
  oncoprint(LUAD);
}


## load gistic 
## TODO it doesn't work
if(!file.exists("input/luadDefGistic.rda")){
  AnalyseDate <- getFirehoseAnalyzeDates(1);
  data <- getFirehoseData("LUAD",
                          gistic2_Date = lastAnalyseDate, 
                          GISTIC = TRUE);
  LUADGistic <- getData(data,
                           type = "GISTIC",
                           platform = "ThresholdedByGene");
  LUADGistic <- as.data.frame(LUADGistic);
  
  LUADGistic <- import.GISTIC(LUADGistic,
                              #filter.genes = genes
                              );
  LUADGistic <- annotate.description(LUADGistic,
                                     label = "LUAD data for driver genes");
  save(LUADGistic, 
       file = "input/luadDefGistic.rda");
  
}else{
  LUADGistic <- loadRData("input/luadDefGistic.rda");
}

## print types for GISTIC
if(verbose){
  print(as.types(LUADGistic));
}

## oncoprint for GISTIC
## TODO wtf
## if(plot_verbose){
##   oncoprint(LUADGistic);
## }

## some changes in LUADGistic
LUADGistic <- delete.type(LUADGistic, 
                          'Heterozygous Loss');
LUADGistic <- delete.type(LUADGistic, 
                          'Low-level Gain');
LUADGistic <- rename.type(LUADGistic, 
                          'Homozygous Loss', 
                          'Deletion');
LUADGistic <- rename.type(LUADGistic, 
                          'High-level Gain', 
                          'Amplification');

## stage annnotation in GISTIC
LUADGistic <- annotate.stages(LUADGistic, 
                              clinical.data)

## print informations for cleared GISTIC
if(verbose){
  print(paste("number of LUADGistic samples:", 
              nsamples(LUADGistic)));    
  print("LUADGistic samples:");
  print(as.samples(LUADGistic)); 
}

## Difference from MAF and GISTIC
if(verbose){
  print("LUAD MAF and GISTIC check:");
  print(as.samples(LUAD)); 
}

## intersect MAF and GISTIC
## TODO if GISTIC not work intersect not work
if(!file.exists("input/luadDefInt.rda")){
  LUADInt <- intersect.datasets(LUADGistic,
                                LUAD, 
                                intersect.genomes = FALSE);
  
  save(LUADInt, 
       file = "input/luadDefInt.rda");
  
}else{
  LUADInt <- loadRData("input/luadDefInt.rda");
}

## oncoprint of intersect
if(plot_verbose){
  oncoprint(LUADInt);
}

# SUBTYPING

## TODO understand if we have cluster (type of tumor like MSI/MSS for CRC)
## and if we can extract them from clinical sheet
## maybe from cBIO
## luad_cbio <- cbio.query(
##   genes = genes[!is.na(genes)],
##   cbio.study = 'luad_tcga',
##   cbio.dataset = 'luad_tcga_3way_complete',
##   cbio.profile = 'luad_tcga_mutations');
## clinical_bio <- luad_cbio$clinical;

# GROUPS EXCLUSIVITY

## print groups from mutex
if(verbose){
  print(LUAD.mutex);
}

## oncoprint first two group
## TODO add pathway feature
if(plot_verbose){
  grid.arrange(
    oncoprint(
      events.selection(LUAD,
                       filter.in.names = LUAD.mutex[[1]]), 
      title = paste("LUAD - Mutex group 1"),
      legend.cex = .3,
      font.row = 6,
      ann.hits = FALSE, 
      cellheight = 10,
      silent = TRUE,    
      # gene.annot = pathway.list,
       gene.annot.color = pathways.color,
      gtable = TRUE
    )$gtable,
    oncoprint(
      events.selection(LUAD,
                       filter.in.names = LUAD.mutex[[2]]), 
      title = paste("LUAD - Mutex group 2"),
      legend.cex = .3,
      silent = TRUE,
      font.row = 6,
      ann.hits = FALSE,
      cellheight = 10,
      #gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      gtable = TRUE
    )$gtable,
    ncol = 1    
  )
}

## apriori knowledge
## TODO now quasi-random genes, RAF genes
LUAD.raf <- c('KRAS', 'NRAS', 'BRAF');

## TODO add MEMo genes (repo offline), now random genes
LUAD.memo <-  c('ERBB2', 'PIK3CA');

## TODO plot not work
if(plot_verbose){
  grid.arrange(
    oncoprint(
      events.selection(LUAD,
                       filter.in.names = LUAD.memo),
      title = paste("LUAD - MEMO exclusivity (knowledge prior)"),
      legend.cex = .3,
      font.row = 6,
      ann.hits = FALSE,
      cellheight = 10,
      cellwidth = 3,
      silent = T,
      #gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      gtable = TRUE
    )$gtable,
    oncoprint(
      events.selection(LUAD,
                       filter.in.names = LUAD.raf), 
      title = paste("LAUD - RAF KRAS/NRAS/BRAF exclusivity (knowledge prior)"),
      legend.cex = .3,
      font.row = 6,
      ann.hits = FALSE,
      cellheight = 10,
      cellwidth = 3,
      silent = T,
      #gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      gtable = TRUE
    )$gtable,
    ncol = 1 
  )
}

# MODELS

## select from LUAD with min freq and apriori
freq <- 0.05;
LUAD.select = select(LUAD, 
                     freq, 
                     unique(             
                       c(LUAD.memo,
                         LUAD.raf,
                         unlist(LUAD.mutex))));
LUAD.select <- annotate.description(LUAD.select,
                                    'LUAD selection');

## oncoprint the selection
## TODO fix pathway
if(plot_verbose){
  oncoprint(LUAD.select, 
            legend.cex = .5,          
            cellwidth = 3,            
            cellheight = 10,
            # gene.annot = pathway.list, 
            gene.annot.color = pathways.color, 
            sample.id = TRUE);
}