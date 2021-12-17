# LOAD DATA

## import mutex data
LUAD.mutex <- import.mutex.groups(file.mutex)

## clear mutex with selected genes
for(i in 1:length(LUAD.mutex)){
  LUAD.mutex[i][[1]] <- LUAD.mutex[i][[1]][LUAD.mutex[i][[1]] %in% pathway.genes]
}

## load MAF, reload if required
if(maf_reload){
  
  ## download MAF file from TCGA
  LUAD.maf <- GDCquery_Maf(tumor = "LUAD", 
                           pipelines = "mutect2")
  
  ## cool visualization of data in MAF
  ## TODO fix errors sometimes
  if(plot_verbose){
    maftools.input <- LUAD.maf %>% read.maf
    plotmafSummary(maf = maftools.input, 
                   rmOutlier = TRUE, 
                   addStat = 'median', 
                   dashboard = TRUE)
    LUAD.maf.genes <- LUAD.maf
    LUAD.maf.genes <- LUAD.maf.genes[LUAD.maf.genes$Hugo_Symbol %in% pathway.genes,]
    maftools.input <- LUAD.maf.genes %>% read.maf
    plotmafSummary(maf = maftools.input, 
                   rmOutlier = TRUE, 
                   addStat = 'median', 
                   dashboard = TRUE)
    
  }
  LUAD.mafdf <- as.data.frame(LUAD.maf)
  
  ## create TRONCO object
  if(all_mut){
    LUAD <- import.MAF(LUAD.mafdf,
                       is.TCGA = TRUE,
                       merge.mutation.types = FALSE,
                       filter.fun = function(x) {
                         return(x['Hugo_Symbol'] %in% pathway.genes)
                       })
  }else{
    LUAD <- import.MAF(LUAD.mafdf,
                       is.TCGA = TRUE,
                       merge.mutation.types = TRUE,
                       filter.fun = function(x) {
                         return(x['Hugo_Symbol'] %in% pathway.genes)
                       })
    ## change color due to labels
    LUAD <- change.color(LUAD, "Mutation", "antiquewhite")
  }
  LUAD <- annotate.description(LUAD, 
                               "Lung cancer data from GDC portal")
  save(LUAD, file = "input/luadDef.rda")
}else{
  LUAD <- loadRData("input/luadDef.rda")
}

## a first and brutal oncoprint
if(plot_verbose){
  oncoprint(LUAD)
}

##  oncoprint after deletions
if(plot_verbose){
  oncoprint(LUAD)
}



## print some informations about genes
if(verbose){
  print(paste("number of LUAD genes:",
              ngenes(LUAD)))    
  print("LUAD genes:")
  print(as.genes(LUAD))                          
  print("LUAD genes for type")
  for (type in as.types(LUAD)) {
    print(paste("gene for type:", 
                type))
    print(as.genes(LUAD,
                   types = type)) 
  }
}

## print some informations about events
if(verbose){
  print(paste("number of LUAD events:",
              nevents(LUAD)))    
  print("LUAD events:")
  print(as.events(LUAD))                       
  ## TODO add other similar events that can be queried together
  print("similar events:")
  print(as.events(LUAD, genes = c('KRAS', 
                                  'TP53')))
  print("LUAD alteration types compacted:")
  print(as.events(LUAD),
        keysToNames = TRUE)
}

## print some informations about types of alteration
if(verbose){
  print(paste("number of LUAD alteration types:", 
              ntypes(LUAD)))    
  print("LUAD alteration types:")
  print(as.types(LUAD)) 
  ## TODO add other similar events that can be queried together
  print("alteration types for similar events:")
  print(as.types(LUAD, genes = c('KRAS', '
                                 TP53')))
}


## print samples information
if(verbose){
  print(paste("number of LUAD samples:", 
              nsamples(LUAD)))    
  print("LUAD samples:")
  print(as.samples(LUAD)) 
  ## TODO add more of these
  if(all_mut){
    print(which.samples(LUAD, 
                        gene = 'KRAS', 
                        type = 'Missense_Mutation'))
  }else{
    print(which.samples(LUAD, 
                        gene = 'KRAS', 
                        type = 'Mutation'))
  }
}

## clinical data
if(clinic_reload){
  data <- getFirehoseData("LUAD")
  clinical <- RTCGAToolbox::getData(data, "clinical")
  df <- as.data.frame(clinical)
  for (i in 1:length(rownames(df))) {
    row.names(df)[i] <- toupper(gsub("\\.", 
                                     "-", 
                                     row.names(df)[i]))
    df$pathologic_stage[i] <- gsub("STAGE", 
                                   "Stage", 
                                   toupper(df$pathologic_stage[i]))
  }
  
  dtutils::write_tsv(df, 
                     file = file.clinical,
                     row_names = TRUE)
  
  smoker <- TCGA.map.clinical.data(file = file.clinical,
                                   column.samples = 'rn',
                                   column.map = 'number_pack_years_smoked')
  
  smoker$number_pack_years_smoked <- round_any(smoker$number_pack_years_smoked,
                                                 10,
                                                 f = ceiling)
  
  smoker <- smoker[order(-smoker$number_pack_years_smoked), ,
                   drop=FALSE]
  
  for(i in 1:length(smoker[[1]])){
    if(!is.na(smoker[[1]][i])){
      if(as.integer(smoker[[1]][i]) < 100){
        smoker[[1]][i]<-paste("0", smoker[[1]][i], sel="")
      }
    }
  }
  
  clinical.data <- TCGA.map.clinical.data(file = file.clinical,
                                          column.samples = 'rn',
                                          column.map = 'pathologic_stage')
  save(clinical.data, file = "input/luadClinical.rda")
}else{
  clinical.data <- loadRData("input/luadClinical.rda")
}

if(verbose){
  print(head(clinical.data))
  ## frequency of gender
  print(data.frame(table(clinical$gender)))
  ## frequency of histolocical type
  ## TODO check if they could be used for a study on subtypes
  ## for: acinar(18), ronchioloalveolar-nonmucinous (19), 
  ## papillary(23), mucinous-colloid (10)
  ## see page 2 marker paper
  print(data.frame(table(clinical$histological_type)))
  ## frequency of ages
  print(data.frame(table(clinical$years_to_birth)))
}

## histogram of age distribution and stages, race and ethnicity
## TODO change colors, lol
if(plot_verbose){
  hist(as.numeric(clinical$years_to_birth), 
       col = c("blue", 
               "green"),
       xlab = "Age",
       main = "Age histogram")
  barplot(table(clinical$pathologic_stage), 
       col = c("blue", 
               "green"),
       xlab = "Stages",
       main = "Stages histogram")
  barplot(table(clinical$race), 
          col = c("blue", 
                  "green"),
          xlab = "Races",
          main = "Races histogram")
  barplot(table(clinical$ethnicity), 
          col = c("blue", 
                  "green"),
          xlab = "Ethnicities",
          main = "Ethnicities histogram")
}


## match samples and stages
LUAD.smoke <- annotate.stages(LUAD, 
                              smoker, 
                              match.TCGA.patients = TRUE)
LUAD <- annotate.stages(LUAD, 
                        clinical.data, 
                        match.TCGA.patients = TRUE)


## clear multiple stages
LUAD <- TCGA.remove.multiple.samples(LUAD)
LUAD <- TCGA.shorten.barcodes(LUAD)
LUAD <- annotate.stages(LUAD, 
                        clinical.data)

## another brutal oncoprint
if(plot_verbose){
  oncoprint(LUAD)
}

## another brutal oncoprint with smoker
if(plot_verbose){
  oncoprint(LUAD.smoke)
}

## other clinicl data from cBIO
luad_cbio <- cbio.query(
  genes = pathway.genes,
  cbio.study = 'luad_tcga',
  cbio.dataset = 'luad_tcga_3way_complete',
  cbio.profile = 'luad_tcga_mutations')
clinical_bio <- luad_cbio$clinical

## here nothing useful
if(verbose){
  ## frequency of histolocical type
  print(data.frame(table(clinical_bio$CANCER_TYPE)))
  ## frequency of histolocical type detailed
  print(data.frame(table(clinical_bio$CANCER_TYPE_DETAILED)))
}

## load gistic 
if(gistic_reload){
  gistic.query <- GDCquery(project = "TCGA-LUAD",
                           data.category = "Copy Number Variation",
                           data.type = "Gene Level Copy Number Scores",
                           access = "open")
  GDCdownload(gistic.query)
  gistic <- GDCprepare(gistic.query)
  gist <- getGistic("LUAD", type = "thresholded")
  LUAD.gistic <- as.data.frame(gist)
  LUAD.gistic <- LUAD.gistic[LUAD.gistic$`Gene Symbol` %in% pathway.genes,]
  LUAD.gistic <- LUAD.gistic[ , ! names(LUAD.gistic) %in% c("Locus ID", 
                                                            "Cytoband")]
  LUAD.gistic <- t(LUAD.gistic)
  colnames(LUAD.gistic) <- lapply(LUAD.gistic[1, ], as.character)
  LUAD.gistic <- LUAD.gistic[-1,]
  LUADGistic <- import.GISTIC(LUAD.gistic,
                              trim = FALSE)
  LUADGistic <- annotate.description(LUADGistic,
                                     label = "LUAD GISTIC for driver genes")
  save(LUADGistic, 
       file = "input/luadDefGistic.rda")
  
}else{
  LUADGistic <- loadRData("input/luadDefGistic.rda")
}

## print types for GISTIC
if(verbose){
  print(as.types(LUADGistic))
}

## oncoprint for GISTIC
if(plot_verbose){
  oncoprint(LUADGistic)
}

## some changes in LUADGistic
LUADGistic <- delete.type(LUADGistic, 
                          'Heterozygous Loss')
LUADGistic <- delete.type(LUADGistic, 
                          'Low-level Gain')
LUADGistic <- rename.type(LUADGistic, 
                          'Homozygous Loss', 
                          'Deletion')
LUADGistic <- rename.type(LUADGistic, 
                          'High-level Gain', 
                          'Amplification')

## oncoprint for reduced GISTIC
if(plot_verbose){
  oncoprint(LUADGistic)
}

## stage annnotation in GISTIC
LUADGistic <- annotate.stages(LUADGistic, 
                              clinical.data,
                              match.TCGA.patients = TRUE)

## clear multiple stages
LUADGistic <- TCGA.remove.multiple.samples(LUADGistic)
LUADGistic <- TCGA.shorten.barcodes(LUADGistic)
LUADGistic <- annotate.stages(LUADGistic, 
                              clinical.data)

## oncoprint for GISTIC with stages
if(plot_verbose){
  oncoprint(LUADGistic)
}

## print informations for cleared GISTIC
if(verbose){
  print(paste("number of LUADGistic samples:", 
              nsamples(LUADGistic)))    
  print("LUADGistic samples:")
  print(as.samples(LUADGistic)) 
}

## Difference from MAF and GISTIC
if(verbose){
  print("LUAD MAF and GISTIC check:")
  print(as.samples(LUAD) %in% as.samples(LUADGistic)) 
}

## intersect MAF and GISTIC and save in LUAD object
if(intersect_reload){
  LUAD <- intersect.datasets(LUADGistic,
                             LUAD, 
                             intersect.genomes = FALSE)
  LUAD <- trim(LUAD)
  LUAD <- annotate.stages(LUAD, 
                          clinical.data)
  LUAD <- annotate.description(x = LUAD,
                               label = "LUAD MAF/CNA data for driver genes")
  save(LUAD, 
       file = "input/luadDefInt.rda")
  ## oncoprint of intersect
  if(plot_verbose){
    oncoprint(LUAD)
  }
  
}else{
  LUAD <- loadRData("input/luadDefInt.rda")
}

## editing types
## TODO check if it ha any sense
## join mutations
if(all_mut){
  LUAD <- join.types(LUAD,
                     "Frame_Shift_Del",
                     "In_Frame_Del",
                     "Frame_Shift_Ins",
                     "In_Frame_Ins",
                     new.type = "Del/Ins_Mutations")
  ## remove quasi-useless mutation type mutations
  # LUAD <- delete.type(LUAD,
  #                     type = "In_Frame_Ins")
  # LUAD <- delete.type(LUAD,
  #                     type = "In_Frame_Del")
  # LUAD <- delete.type(LUAD,
  #                     type = "Frame_Shift_Ins")
  # LUAD <- delete.type(LUAD,
  #                     type = "Frame_Shift_Del")
  LUAD <- delete.type(LUAD,
                      type = "3'UTR")
  LUAD <- delete.type(LUAD,
                      type = "5'UTR")
  LUAD <- delete.type(LUAD,
                      type = "Splice_Region")
  LUAD <- delete.type(LUAD,
                      type = "Splice_Site")
  LUAD <- delete.type(LUAD,
                      type = "Intron")
  LUAD <- delete.type(LUAD,
                      type = "Silent")
}
## select only event with a minfreq
LUAD <- events.selection(LUAD, 
                         filter.freq = min_freq)
LUAD.smoke <- annotate.stages(LUAD, 
                              smoker, 
                              match.TCGA.patients = TRUE)


## oncoprint of intersect with selection
if(plot_verbose){
  oncoprint(LUAD)
  ## TODO make it visualizable
  oncoprint(LUAD,
            legend.cex = .5,          
            cellwidth = 2,           
            cellheight = 2,
            gene.annot = pathway.list, 
            gene.annot.color = pathways.color)
}

## oncoprint of intersect with selection and smoker
if(plot_verbose){
  oncoprint(LUAD.smoke)
}


## other fancy plots
MAF.dataframe <- import.MAF(LUAD.mafdf,
                            is.TCGA = TRUE,
                            sep = ';',
                            to.TRONCO = FALSE,
                            #merge.mutation.types = TRUE,
                            filter.fun = function(x) {
                              return(x['Hugo_Symbol'] %in% pathway.genes)
                            }) 
MAF.dataframe <-
  MAF.dataframe[which(MAF.dataframe$Variant_Classification !=
                        'Splice_Region'),] 

mut.colors <- brewer.pal('Set3', n = 11)
names(mut.colors) <- unique(MAF.dataframe$Variant_Classification)

waterfall(MAF.dataframe,
          mainGrid = T,
          mainDropMut = T,
          mainPalette = mut.colors)
