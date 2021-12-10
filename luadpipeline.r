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
source("src/utils.r")
source("src/conf.r")

# LOAD DATA

## gene selection

## use pathway as in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4231481/
## complete genes
# P53 <- c("TP53", "ATM", "MDM2");
# RAS <- c("KRAS", "NRAS", "HRAS", "RIT1", "NF1", "BRAF", "MAP2K1",
#          "EGFR", "ERBB2", "MET", "ALK", "RET", "ROS1");
# MTOR <- c("PTEN", "PIK3CA", "PIK3R1", "STK11", "AKT1", "AMPK",
#           "TSC1", "TSC2", "MTOR");
# KEAP1 <- c("KEAP1", "CUL3", "NFE2L2");
# CDKN2A <- c("CDKN2A", "CCND1", "CDK4", "CCNE1", "RB1");
# ARID <- c("ARID1A", "ARID1B", "ARID2", "SMARCA4");
# SETD2 <- c("SETD2");
# RBM10 <- c("RBM10", "U2AF1");
## genes > 1%
P53 <- c("TP53", "ATM", "MDM2");
RAS <- c("KRAS", "RIT1", "NF1", "BRAF", "EGFR", "ERBB2", "MET");
MTOR <- c("PTEN", "PIK3CA", "STK11", "TSC1", "TSC2", "MTOR");
KEAP1 <- c("KEAP1", "NFE2L2");
CDKN2A <- c("CDKN2A", "CCND1", "CDK4", "CCNE1", "RB1");
ARID <- c("ARID1A", "ARID1B", "ARID2", "SMARCA4");
SETD2 <- c("SETD2");
RBM10 <- c("RBM10", "U2AF1");


## create pathways
pathway.genes <- c(P53, RAS, MTOR, KEAP1, CDKN2A, ARID, SETD2, RBM10);
pathway.genes <- unique(pathway.genes);
pathway.names <- c("P53", "RAS", "MTOR", "KEAP1", "CDKN2A", "ARID", 
                   "SETD2", "RBM10")
pathway.list <- list(P53 = P53,
                     RAS = RAS,
                     MTOR = MTOR, 
                     KEAP1 = KEAP1, 
                     CDKN2A = CDKN2A, 
                     ARID = ARID, 
                     SETD2 = SETD2, 
                     RBM10 = RBM10);

## colors for pathways
alteration.color = 'dimgray';
pathways.color = c('firebrick1', 
                   'darkblue', 
                   'darkgreen',
                   'darkmagenta', 
                   'darkorange',
                   'dodgerblue4',
                   'darkorchid1',
                   'darksalmon');

## import mutex data
LUAD.mutex <- import.mutex.groups(file.mutex);

## load MAF, reload if required
if(maf_reload){
  
  ## download MAF file from TCGA
  LUAD.maf <- GDCquery_Maf(tumor = "LUAD", 
                           pipelines = "mutect2");
  ## TODO fix errors
  if(plot_verbose && FALSE){
    maftools.input <- LUAD.maf %>% read.maf;
    plotmafSummary(maf = maftools.input, 
                   rmOutlier = TRUE, 
                   addStat = 'median', 
                   dashboard = TRUE);
  }
  LUAD.mafdf <- as.data.frame(LUAD.maf);
  
  ## create TRONCO object
  LUAD <- import.MAF(LUAD.mafdf,
                    is.TCGA = TRUE,
                    merge.mutation.types = FALSE,
                    filter.fun = function(x) {
                      return(x['Hugo_Symbol'] %in% pathway.genes)
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

##  oncoprint after deletions
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
if(clinic_reload){
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
if(gistic_reload){
  gistic.query <- GDCquery(project = "TCGA-LUAD",
                           data.category = "Copy Number Variation",
                           data.type = "Gene Level Copy Number Scores",
                           access = "open");
  GDCdownload(gistic.query);
  gistic <- GDCprepare(gistic.query);
  gist <- getGistic("LUAD", type = "thresholded");
  LUAD.gistic <- as.data.frame(gist);
  LUAD.gistic <- LUAD.gistic[LUAD.gistic$`Gene Symbol` %in% pathway.genes,];
  LUAD.gistic <- LUAD.gistic[ , ! names(LUAD.gistic) %in% c("Locus ID", 
                                                            "Cytoband")];
  LUAD.gistic <- t(LUAD.gistic);
  colnames(LUAD.gistic) <- lapply(LUAD.gistic[1, ], as.character);
  LUAD.gistic <- LUAD.gistic[-1,];
  LUADGistic <- import.GISTIC(LUAD.gistic,
                              trim = FALSE);
  LUADGistic <- annotate.description(LUADGistic,
                                     label = "LUAD GISTIC for driver genes");
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
if(plot_verbose){
  oncoprint(LUADGistic);
}

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
                              clinical.data,
                              match.TCGA.patients = TRUE);

## clear multiple stages
LUADGistic <- TCGA.remove.multiple.samples(LUADGistic);
LUADGistic <- TCGA.shorten.barcodes(LUADGistic);
LUADGistic <- annotate.stages(LUADGistic, 
                        clinical.data);

## oncoprint for GISTIC with stages
if(plot_verbose){
  oncoprint(LUADGistic);
}

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
  print(as.samples(LUAD) %in% as.samples(LUADGistic)); 
}

## intersect MAF and GISTIC and save in LUAD object
if(intersect_reload){
  LUAD <- intersect.datasets(LUADGistic,
                             LUAD, 
                             intersect.genomes = FALSE);
  LUAD <- trim(LUAD);
  LUAD <- annotate.stages(LUAD, 
                          clinical.data);
  LUAD <- annotate.description(x = LUAD,
                               label = "LUAD MAF/CNA data for driver genes");
  save(LUAD, 
       file = "input/luadDefInt.rda");
  ## oncoprint of intersect
  if(plot_verbose){
    oncoprint(LUAD);
  }
  
}else{
  LUAD <- loadRData("input/luadDefInt.rda");
}

## join mutations
LUAD <- join.types(LUAD,
                   "Frame_Shift_Del",
                   "In_Frame_Del",
                   "Frame_Shift_Ins",
                   "In_Frame_Ins",
                   new.type = "Del/Ins_Mutations");
## remove quasi-useless mutation type mutations
# LUAD <- delete.type(LUAD,
#                     type = "In_Frame_Ins");
# LUAD <- delete.type(LUAD,
#                     type = "In_Frame_Del");
# LUAD <- delete.type(LUAD,
#                     type = "Frame_Shift_Ins");
# LUAD <- delete.type(LUAD,
#                     type = "Frame_Shift_Del");
LUAD <- delete.type(LUAD,
                    type = "3'UTR");
LUAD <- delete.type(LUAD,
                    type = "5'UTR");
LUAD <- delete.type(LUAD,
                    type = "Splice_Region");
LUAD <- delete.type(LUAD,
                    type = "Splice_Site");
LUAD <- delete.type(LUAD,
                    type = "Intron");
LUAD <- delete.type(LUAD,
                    type = "Silent");

## select only event with a minfreq
LUAD <- events.selection(LUAD, 
                         filter.freq = min_freq);

## oncoprint of intersect with selection
if(plot_verbose){
  oncoprint(LUAD);
}

# SUBTYPING

## TODO understand if we have cluster (type of tumor like MSI/MSS for CRC)
## and if we can extract them from clinical sheet
## maybe from cBIO
## luad_cbio <- cbio.query(
##   pathway.genes = pathway.genes[!is.na(pathway.genes)],
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
      gene.annot = pathway.list,
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
      gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      gtable = TRUE
    )$gtable,
    ncol = 1    
  )
}

## apriori knowledge
## TODO now quasi-random genes, RAF genes and MTOR genes
LUAD.raf <- c('KRAS', 'EGFR');
LUAD.mtor <-  c('PIK3CA', 'STK11');

## TODO plot not work
if(plot_verbose && FALSE){
  grid.arrange(
    oncoprint(
      events.selection(LUAD,
                       filter.in.names = LUAD.mtor),
      title = paste("LUAD - MEMO exclusivity (knowledge prior)"),
      legend.cex = .3,
      font.row = 6,
      ann.hits = FALSE,
      cellheight = 10,
      cellwidth = 3,
      silent = T,
      gene.annot = pathway.list,
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
      gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      gtable = TRUE
    )$gtable,
    ncol = 1 
  )
}

# MODELS

## select from LUAD with min freq and apriori
LUAD.select = select(LUAD, 
                     min_freq, 
                     unique(             
                       c(LUAD.mtor,
                         LUAD.raf,
                         unlist(LUAD.mutex))));
LUAD.select <- annotate.description(LUAD.select,
                                    'LUAD selection');

## oncoprint the selection
if(plot_verbose){
  oncoprint(LUAD.select, 
            legend.cex = .5,          
            cellwidth = 3,            
            cellheight = 10,
            gene.annot = pathway.list, 
            gene.annot.color = pathways.color, 
            sample.id = TRUE);
}

## consolidate dataset
del <- consolidate.data(LUAD.select, 
                        print = TRUE);
if(length(del[["indistinguishable"]]) > 0){
  for (i in 1:length(del[["indistinguishable"]])) {
    for (j in 1:nrow(del[["indistinguishable"]][[i]])){
      gene <- del[["indistinguishable"]][[i]][j,][2][[1]];
      type <- del[["indistinguishable"]][[i]][j,][1][[1]];
      LUAD.select <- delete.event(LUAD.select,
                                  gene = as.character(gene),
                                  type = as.character(type));
    }
  }
}


## add hypotheses

if(hypo_reload){
  LUAD.hypo <- LUAD.select;
  
  ## first hypotheses from mutex (using only available genes)
  if (!is.null(LUAD.mutex)) {
    for (group in LUAD.mutex) {
      group <- group[group%in% as.genes(LUAD.hypo)];
      if(length(group) >= 2){
        LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
                                          FUN = OR,  
                                          group = group,
                                          dim.min = length(group));
      }
      
    }
  }
  
  
  ## ADD hypothes for RAF, checking if we have the genes in LUAD.select
  LUAD.raf.subtype <- LUAD.raf[LUAD.raf%in% as.genes(LUAD.hypo)];
  LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
                                    FUN = OR, 
                                    group = LUAD.raf.subtype, 
                                    dim.min = length(LUAD.raf.subtype)); 
  
  ## then for MTOR group
  LUAD.mtor.subtype <- LUAD.mtor[LUAD.mtor%in% as.genes(LUAD.hypo)];
  LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
                                    FUN = OR, 
                                    group = LUAD.mtor.subtype, 
                                    dim.min = length(LUAD.mtor.subtype));
  
  ## add all the hypotheses related to homologou events
  LUAD.hypo <- hypothesis.add.homologous(LUAD.hypo);
  
  ## add annotation
  LUAD.hypo <- annotate.description(LUAD.hypo, 
                                    as.description(LUAD.select));
  
  
  ## first use of CAPRI
  LUAD.model <- tronco.capri(LUAD.hypo, 
                             boot.seed = 12345,
                             nboot = num_boot_iter);
  
  ## DAG of model with hypotheses
  if(plot_verbose){
    tronco.plot(LUAD.model, 
                pathways = pathway.list,  
                edge.cex = 1.5,          
                legend.cex = .35, 
                scale.nodes = .5,        
                confidence = c('tp', 'pr', 'hg'), 
                pathways.color = pathways.color,  
                disconnected = F,        
                height.logic = .3);
  }
  
  ## save data
  save(LUAD.hypo, 
       file = "input/luadDefHypo.rda");
  save(LUAD.model, 
       file = "input/luadDefHypoModel.rda");
  
}else{
  LUAD.hypo <- loadRData("input/luadDefHypo.rda");
  LUAD.model <- loadRData("input/luadDefHypoModel.rda");
}

## dataframe with selective advanges, with fit probabilities, optimized
LUAD.hypo.model.selfit <- as.selective.advantage.relations(LUAD.model);

if(verbose){
  print("advatanges selection fit probabilities");
  print(LUAD.hypo.model.selfit);
}

## dataframe with selective advanges, with prima facie, full set of edge
LUAD.hypo.model.selpf <- as.selective.advantage.relations(LUAD.model,
                                                          type = "pf");

if(verbose){
  print("advatanges selection full set");
  print(LUAD.hypo.model.selpf);
}

## dataframe with selective advanges, with a subset of genes
## TODO make test with usefull subset of genes
LUAD.hypo.model.selsub <- as.selective.advantage.relations(LUAD.model,
                                                           events = as.events(LUAD.model, 
                                                                              genes = P53));
if(verbose){
  print("advatanges selection by pathway");
  print(LUAD.hypo.model.selsub);
}

# STATISTICS

## non-parametric bootstrap
if(boot_reload){
  LUAD.model <- tronco.bootstrap(LUAD.model,
                                 nboot = num_boot_iter,
                                 cores.ratio = .5);
  
  ## statistical bootstrap
  LUAD.model <- tronco.bootstrap(LUAD.model,
                                 type = "statistical",
                                 nboot = num_boot_iter,
                                 cores.ratio = .5);
  
  save(LUAD.model, file = "input/luadDefBoot.rda");
}else{
  LUAD.model <- loadRData("input/luadDefBoot.rda");
}

## DAG of the model above
if(plot_verbose){
  tronco.plot(LUAD.model, 
              pathways = pathway.list,  
              edge.cex = 1.5,          
              legend.cex = .5,         
              scale.nodes = .6,        
              confidence = c('tp', 'pr', 'hg'), 
              pathways.color = pathways.color,  
              disconnected = F,        
              height.logic = .3,       
              #file = "output/model_hypo_boot.pdf"
              );
}

## plot of bootstrap scores
if(plot_verbose){
  ## first non-parametric
  pheatmap(keysToNames(LUAD.model,
                       as.confidence(LUAD.model,
                                     conf = 'npb')$npb$capri_bic) * 100,
           main = "non-parametric bootstrap scores for BIC model",
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%f"
  );
  pheatmap(keysToNames(LUAD.model,
                       as.confidence(LUAD.model,
                                     conf = 'npb')$npb$capri_aic) * 100,
           main = "non-parametric bootstrap scores for AIC model",
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%f"
  );
  ## then parametric ones
  pheatmap(keysToNames(LUAD.model,
                       as.confidence(LUAD.model,
                                     conf = 'sb')$sb$capri_bic) * 100,
           main = "non-parametric bootstrap scores for BIC model",
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%f"
  );
  pheatmap(keysToNames(LUAD.model,
                       as.confidence(LUAD.model,
                                     conf = 'sb')$sb$capri_aic) * 100,
           main = "non-parametric bootstrap scores for AIC model",
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%f"
  );
}

## table with bootstrap scores
boot_tab <- as.bootstrap.scores(LUAD.model);

if(verbose){
  print(boot_tab);
}

## kfold
if(kfold_reload){
  ## k-fold cross validation, prediction error for each parent set X
  LUAD.model <- tronco.kfold.eloss(LUAD.model);
  kfold_eloss <- as.kfold.eloss(LUAD.model);
  
  ## plot for every fold
  ## TODO make it work
  if(plot_verbose && FALSE){
    vioplot(LUAD.model$kfold$capri_bic$eloss,
            LUAD.model$kfold$capri_aic$eloss,
            col = 'red',
            lty = 1, rectCol="gray",
            colMed = 'black',
            names = c('BIC', 'AIC'), 
            pchMed = 15, 
            horizontal = T)
    title(main = 'Entropy loss \n LUAD tumors');
  }
  
  ## k-fold cross validation, prediction error for each parent set X
  LUAD.model <- tronco.kfold.prederr(LUAD.model);
  kfold_pred <- as.kfold.prederr(LUAD.model);
  
  
  ## k-fold cross validation, posterior classification error for each edge
  LUAD.model <- tronco.kfold.posterr(LUAD.model);
  kfold_post <- as.kfold.posterr(LUAD.model);
  save(LUAD.model, 
       file = "input/luadDefKfold.rda");
  
}else{
  LUAD.model <- loadRData("input/luadDefKfold.rda");
}

## visualize a table with all edge statistics
tab_bic <- tabular(LUAD.model, 'capri_bic');
tab_aic <- tabular(LUAD.model, 'capri_aic');

if(verbose){
  print("table with all edge statistics using capri_bic");
  print(tab_bic);
  print("table with all edge statistics using capri_aic");
  print(tab_aic)
}

## save model
save(LUAD.model, 
     file = "input/luadDefModel.rda");

## TODO add some graph regarding pattern

## last DAG
if(plot_verbose){
  tronco.plot(LUAD.model, 
              pathways = pathway.list,  
              edge.cex = 1.5,          
              legend.cex = .35,         
              scale.nodes = .6,        
              confidence = c('tp', 'pr', 'hg'), 
              pathways.color = pathways.color,  
              disconnected = F,        
              height.logic = .3,       
              #file = "output/model_hypo_boot.pdf"
  );
}

## excel with all data
excel.file = "output/LUAD_statistics.xlsx";

excel.wbook = createWorkbook();

sheet.luad.bic <- createSheet(wb = excel.wbook, 
                              sheetName="LUAD-bic");
sheet.luad.aic <- createSheet(wb = excel.wbook, 
                              sheetName="LUAD-aic");

addDataFrame(x = tabular(LUAD.model, 
                         'capri_bic'),
             sheet = sheet.luad.bic,
             showNA = T,
             characterNA = 'NA');
addDataFrame(x = tabular(LUAD.model, 
                         'capri_aic'),
             sheet = sheet.luad.aic,
             showNA = T,
             characterNA = 'NA');

saveWorkbook(excel.wbook, 
             excel.file)
