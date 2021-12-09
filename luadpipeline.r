library(TRONCO)
library(RTCGAToolbox)
library(TCGAbiolinks)
library(readr)
library(dtutils)
library(ggplot2)
library(gridExtra)
library("xlsx")
source("src/utils.r")
source("src/conf.r")

# LOAD DATA

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

## TODO use pathway such as
P53 = c("TP53", "ATM")

## load data if not exist a rda
if(!file.exists("input/luadDef.rda")){
  
  ## download MAF file from TCGA
  LUAD.maf <- GDCquery_Maf(tumor = "LUAD", 
                           pipelines = "mutect2");
  LUAD.maf <- as.data.frame(LUAD.maf);
  
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

if(gistic_check){
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
LUAD.select = select(LUAD, 
                     min_freq, 
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
  
  ## first hypotheses from mutex
  if (!is.null(LUAD.mutex)) {
    for (group in LUAD.mutex) {
      LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
                                        FUN = OR,  
                                        group = group,
                                        dim.min = length(group));
      
    }
  }
  
  
  ## the raf,, checking if we have the genes in LUAD.select
  LUAD.raf.subtype <- LUAD.raf[LUAD.raf%in% as.genes(LUAD.hypo)];
  LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
                                    FUN = OR, 
                                    group = LUAD.raf.subtype, 
                                    dim.min = length(LUAD.raf.subtype)); 
  
  ## then for MEMO group
  LUAD.memo.subtype <- LUAD.memo[LUAD.memo%in% as.genes(LUAD.hypo)];
  LUAD.hypo <- hypothesis.add.group(LUAD.hypo, 
                                    FUN = OR, 
                                    group = LUAD.memo.subtype, 
                                    dim.min = length(LUAD.memo.subtype));
  
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
  ## TODO fix pathways
  if(plot_verbose){
    tronco.plot(LUAD.model, 
                #pathways = pathway.list,  
                edge.cex = 1.5,          
                legend.cex = .5,         
                scale.nodes = .6,        
                confidence = c('tp', 'pr', 'hg'), 
                pathways.color = pathways.color,  
                disconnected = F,        
                height.logic = .3,       
                file = "output/model_hypo.pdf");
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
LUAD.hypo.selfit <- as.selective.advantage.relations(LUAD.model);

## dataframe with selective advanges, with prima facie, full set of edge
LUAD.hypo.selpf <- as.selective.advantage.relations(LUAD.model,
                                                    type = "pf");

## dataframe with selective advanges, with a subset of genes
## TODO make test with usefull subset of genes
LUAD.hypo.selsub <- as.selective.advantage.relations(LUAD.model,
                                                     events = as.events(LUAD.model, 
                                                                        genes = P53));
## DAG of the model above
## TODO fix pathways
if(plot_verbose){
  tronco.plot(LUAD.model.hypo.selfit, 
              #pathways = pathway.list,  
              edge.cex = 1.5,          
              legend.cex = .5,         
              scale.nodes = .6,        
              confidence = c('tp', 'pr', 'hg'), 
              pathways.color = pathways.color,  
              disconnected = F,        
              height.logic = .3,       
              file = "output/model_hypo_selfit.pdf");
  tronco.plot(LUAD.model.hypo.selpf, 
              #pathways = pathway.list,  
              edge.cex = 1.5,          
              legend.cex = .5,         
              scale.nodes = .6,        
              confidence = c('tp', 'pr', 'hg'), 
              pathways.color = pathways.color,  
              disconnected = F,        
              height.logic = .3,       
              file = "output/model_hypo_selpf.pdf");
  tronco.plot(LUAD.model.hypo.selsub, 
              #pathways = pathway.list,  
              edge.cex = 1.5,          
              legend.cex = .5,         
              scale.nodes = .6,        
              confidence = c('tp', 'pr', 'hg'), 
              pathways.color = pathways.color,  
              disconnected = F,        
              height.logic = .3,       
              file = "output/model_hypo_selsub.pdf");
}

# STATISTICS

## non-parametric bootstrap
LUAD.model <- tronco.bootstrap(LUAD.model,
                                    nboot = num_boot_iter,
                                    cores.ratio = .5);

## statistical bootstrap
LUAD.model <- tronco.bootstrap(LUAD.model,
                                    type = "statistical",
                                    nboot = num_boot_iter,
                                    cores.ratio = .5);

## DAG of the model above
## TODO fix pathways
if(plot_verbose){
  tronco.plot(LUAD.model.hypo, 
              #pathways = pathway.list,  
              edge.cex = 1.5,          
              legend.cex = .5,         
              scale.nodes = .6,        
              confidence = c('tp', 'pr', 'hg'), 
              pathways.color = pathways.color,  
              disconnected = F,        
              height.logic = .3,       
              file = "output/model_hypo_boot.pdf");
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
                                     conf = 'sb')$npb$capri_bic) * 100,
           main = "non-parametric bootstrap scores for BIC model",
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = T,
           number_format = "%f"
  );
  pheatmap(keysToNames(LUAD.model,
                       as.confidence(LUAD.model,
                                     conf = 'sb')$npb$capri_aic) * 100,
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

## k-fold cross validation, prediction error for each parent set X
LUAD.model <- tronco.kfold.eloss(LUAD.model);
kfold_eloss<- as.kfold.eloss(LUAD.model);

## plot for every fold
if(plot_verbose){
  vioplot(LUAD.model$kfold$capri_bic$eloss,
          LUAD.model$kfold$aic$eloss,
          col = 'red',
          lty = 1, rectCol="gray",
          colMed = 'black',
          names = c('BIC', 'AIC'), 
          pchMed = 15, 
          horizontal = T)
  title(main = 'Entropy loss \n LUAD tumors');
  
  vioplot(LUAD.model$kfold$capri_bic$eloss,
          LUAD.model$kfold$capri_aic$eloss,
          col = 'red',
          lty = 1,
          rectCol="gray",
          colMed = 'black',
          names = c('BIC', 'AIC'),
          pchMed = 15, horizontal = T)
  title(main = 'Entropy loss \n LUAD tumors');
}

## k-fold cross validation, prediction error for each parent set X
LUAD.model <- tronco.kfold.prederr(LUAD.model);
kfold_pred <- as.kfold.prederr(LUAD.model);


## k-fold cross validation, posterior classification error for each edge
LUAD.model <- tronco.kfold.posterr(LUAD.model);
kfold_post <- as.kfold.posterr(LUAD.model);

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
             excel.file);