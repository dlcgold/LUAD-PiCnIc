library("xlsx")
library(TRONCO)
library(RTCGAToolbox)
library(TCGAbiolinks)


loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}


#' input_genes_gdc
#'
#' function to load desired data, directly from TCGA using GDC, default genes 
#' are IMCDriver genes driver
#' @param interactive bool to display some options for loading data
#' @param genes_custom custom genes array to use, if not empty interactive is 
#' disabled
#'
#' @return LUAD dataset for TRONCO
#' @export
#'
#' @examples
#' # simple load without custom array
#' input_genes()
#' # no interactive mode
#' input_genes(interactive = FALSE)
#' # load with custom genes array
#' input_genes(gene_custom = c('TP53', 'KRAS', 'KEAP11', 'STK11'))
input_genes_gdc <- function(interactive = TRUE, genes_custom = c()) {
    
    ## use custom genes
    if(length(genes_custom) != 0){
        LUADmut <- GDCquery_Maf(tumor = "LUAD", pipelines = "mutect2");
        LUADmut <- as.data.frame(LUADmut);
        LUADclear <- subset(LUADmut, Hugo_Symbol %in% genes_custom);
        cols_exc <- names(LUADclear) %in% c("Hugo_Symbol", 
                                            "Tumor_Sample_Barcode",
                                            "Variant_Classification");
        LUADclear <- LUADclear[cols_exc];
        ##LUADclearr <- reduce_samples(LUADclear, 9);
        
        LUAD = import.MAF(LUADclear,
                          merge.mutation.types = FALSE, 
                          irregular = TRUE);
        
        LUAD <- annotate.description(LUAD, 
                                     "Lung cancer data from GDC portal");
        return(LUAD);
    }
    
    if(interactive){
        ## choose if use default drivers
        c <- readline(prompt = "Use default gene drivers (IMC)? [yes/no] ");
        c <- as.character(c);
        while(c != "yes" && c != "no"){
            c <- readline(prompt = "Use default gene drivers? [yes/no] ");
            c <- as.character(c);
        }
        if (c == "yes"){
            choice_bool <- TRUE;
        }else{
            choice_bool <- FALSE;
        }
        
        if(!choice_bool){
            ## if not default load xlsx file from suppl. 4
            file_drivers <- "input/gene_drivers.xlsx";
            wb <- loadWorkbook(file_drivers);
            sheets <- names(getSheets(wb));
            count <- 1;
            ## choose desired software
            print("Select gene driver software from");
            for (elem in sheets) {
                if(elem == "IMCDriver"){
                    cat(elem, "(suggested)", count, "\n");
                }else{
                    cat(elem, count, "\n");
                }
                count <- count + 1
            }
            
            choice_raw <- readline(prompt = "Enter your choice: ");
            choice <- as.integer(choice_raw);
            choice_str <- as.character(choice);
            
            while(choice_raw != choice_str ||
                  choice <= 0 ||
                  choice > length(sheets)){
                      choice_raw <- readline(prompt = "Enter your choice: "); 
                      choice <- as.integer(choice_raw);
                      choice_str <- as.character(choice);
                  }
            ## load genes
            genes <- read.xlsx(file_drivers, 
                               sheetIndex = choice, 
                               header = TRUE);
            
            LUADmut <- GDCquery_Maf(tumor = "LUAD", pipelines = "mutect2");
            LUADmut <- as.data.frame(LUADmut);
            genes <- read.xlsx("input/gene_drivers.xlsx", 
                               sheetIndex = 1, 
                               header = TRUE);
            LUADclear <- subset(LUADmut, Hugo_Symbol %in% genes$LUAD);
            ## cols_exc <- names(LUADclear) %in% c("Hugo_Symbol", 
            ##                                     "Tumor_Sample_Barcode",
            ##                                     "Variant_Classification");
            ## LUADclear <- LUADclear[cols_exc];
            ## LUADclearr <- reduce_samples(LUADclear, 9);
            
            LUAD = import.MAF(LUADclear,
                              merge.mutation.types = FALSE);
            LUAD <- annotate.description(LUAD, 
                                         "Lung cancer data from GDC portal");
            ##save(LUAD, file = "input/luadIMCgdc.rda");
            return(LUAD);
        }else{
            ## default genes already in a RDATA file
            LUAD <- loadRData("input/luadIMCgdc.rda");
            return (LUAD);
        }
    }else{
        ## default genes already in a RDATA file
        LUAD <- loadRData("input/luadIMCgdc.rda");
        return(LUAD);
    }
}


#' input_genes_cbio
#'
#' function to load desired data, using TRONCO cBIO interface, default genes are 
#' IMCDriver genes driver
#' @param interactive bool to display some options for loading data
#' @param genes_custom custom genes array to use, if not empty interactive is 
#' disabled
#'
#' @return LUAD dataset for TRONCO
#' @export
#'
#' @examples
#' # simple load without custom array
#' input_genes()
#' # no interactive mode
#' input_genes(interactive = FALSE)
#' # load with custom genes array
#' input_genes(gene_custom = c('TP53', 'KRAS', 'KEAP11', 'STK11'))
input_genes_cbio <- function(interactive = TRUE, genes_custom = c()) {

    ## use custom genes
    if(length(genes_custom) != 0){
        luad_data <- cbio.query(
            genes = genes_custom[!is.na(genes_custom)],
            cbio.study = 'luad_tcga_pub',
            cbio.dataset = 'luad_tcga_pub_cnaseq',
            cbio.profile = 'luad_tcga_pub_mutations');
        
        m <- as.matrix(luad_data$profile);
        m[is.na(m)] <- 0;
        m[m == 'NaN'] = 0;
        m[m != '0'] <- 1;
        LUAD <- import.genotypes(m, 
                                 event.type = 'Mutation', 
                                 color = 'brown3');
        
        LUAD <- annotate.description(LUAD, 
                                     "Lung cancer data from Cbio portal");
        return(LUAD);
    }
    
    if(interactive){
        ## choose if use default drivers
        c <- readline(prompt = "Use default gene drivers (IMC)? [yes/no] ");
        c <- as.character(c);
        while(c != "yes" && c != "no"){
            c <- readline(prompt = "Use default gene drivers? [yes/no] ");
            c <- as.character(c);
        }
        if (c == "yes"){
            choice_bool <- TRUE;
        }else{
            choice_bool <- FALSE;
        }
        
        if(!choice_bool){
            ## if not default load xlsx file from suppl. 4
            file_drivers <- "input/gene_drivers.xlsx";
            wb <- loadWorkbook(file_drivers);
            sheets <- names(getSheets(wb));
            count <- 1;
            ## choose desired software
            print("Select gene driver software from");
            for (elem in sheets) {
                if(elem == "IMCDriver"){
                    cat(elem, "(suggested)", count, "\n");
                }else{
                    cat(elem, count, "\n");
                }
                count <- count + 1
            }
            
            choice_raw <- readline(prompt = "Enter your choice: ");
            choice <- as.integer(choice_raw);
            choice_str <- as.character(choice);
            
            while(choice_raw != choice_str ||
                  choice <= 0 ||
                  choice > length(sheets)){
                      choice_raw <- readline(prompt = "Enter your choice: "); 
                      choice <- as.integer(choice_raw);
                      choice_str <- as.character(choice);
                  }
            ## load genes
            genes <- read.xlsx(file_drivers, 
                               sheetIndex = choice, 
                               header = TRUE);
            
            luad_data <- cbio.query(
                genes = genes$LUAD[!is.na(genes$LUAD)],
                cbio.study = 'luad_tcga',
                cbio.dataset = 'luad_tcga_3way_complete',
                cbio.profile = 'luad_tcga_mutations');
            
            ## prepare data for TRONCO usage
            m <- as.matrix(luad_data$profile);
            m[is.na(m)] <- 0;
            m[m == 'NaN'] = 0;
            m[m != '0'] <- 1
            LUAD_mut <- import.genotypes(m, 
                                         event.type = 'Mutation', 
                                         color = 'brown3');
            
            ## luad_data <- cbio.query(
            ##     genes = genes$LUAD[!is.na(genes$LUAD)],
            ##     cbio.study = 'luad_tcga',
            ##     cbio.dataset = 'luad_tcga_3way_complete',
            ##     cbio.profile = 'luad_tcga_gistic');
            
            ## m <- as.matrix(luad_data$profile);
            ## m[is.na(m)] <- 0;
            ## m[m == 'NaN'] = 0;
            ## m[m != '0'] <- 1
            
            ## LUAD_cna <- import.genotypes(m, 
            ##                              event.type = 'CNA', 
            ##                              color = 'blue3');
            ## luad_data <- cbio.query(
            ##     genes = genes$LUAD[!is.na(genes$LUAD)],
            ##     cbio.study = 'luad_tcga',
            ##     cbio.dataset = 'luad_tcga_3way_complete',
            ##     cbio.profile = 'luad_tcga_methylation_hm27');
            
            ## m <- as.matrix(luad_data$profile);
            ## m[is.na(m)] <- 0;
            ## m[m == 'NaN'] = 0;
            ## m[m != '0'] <- 1
            
            ## LUAD_hm27 <- import.genotypes(m, 
            ##                               event.type = 'hm27', 
            ##                               color = 'green3');
            ## LUAD <- intersect.datasets(LUAD_mut,
            ##                            LUAD_cna, 
            ##                            intersect.genomes = FALSE);
            ## LUAD <- intersect.datasets(LUAD,
            ##                            LUAD_hm27, 
            ##                            intersect.genomes = FALSE);
            LUAD <- annotate.description(LUAD, 
                                         "Lung cancer data from Cbio portal");
            ## save(LUAD, file = "input/luadIMCcbio.rda");
            return(LUAD);
        }else{
            ## default genes already in a RDATA file
            LUAD <- loadRData("input/luadIMCcbio.rda");
            return (LUAD);
        }
    }else{
        ## default genes already in a RDATA file
        LUAD <- loadRData("input/luadIMCcbio.rda");
        return(LUAD);
    }
}

input_gistic <- function(interactive = TRUE, genes_custom = c()) {
    
    ## use custom genes
    if(length(genes_custom) != 0){
        lastAnalyseDate <- getFirehoseAnalyzeDates(1);
        
        gistic <- getFirehoseData("LUAD",
                                  gistic2_Date = lastAnalyseDate, 
                                  GISTIC = TRUE);
        ##gistic <- as.data.frame(gistic);
        ## get GISTIC results
        LUADGistic <- getData(gistic,
                              type = "GISTIC", 
                              platform = "AllByGene");


        LUAD <- import.GISTIC(LUADGistic,
                              filter.genes = genes_custom);
        return(LUAD);
    }
    
    if(interactive){
        ## choose if use default drivers
        c <- readline(prompt = "Use default gene drivers (IMC)? [yes/no] ");
        c <- as.character(c);
        while(c != "yes" && c != "no"){
            c <- readline(prompt = "Use default gene drivers? [yes/no] ");
            c <- as.character(c);
        }
        if (c == "yes"){
            choice_bool <- TRUE;
        }else{
            choice_bool <- FALSE;
        }
        
        if(!choice_bool){
            ## if not default load xlsx file from suppl. 4
            file_drivers <- "input/gene_drivers.xlsx";
            wb <- loadWorkbook(file_drivers);
            sheets <- names(getSheets(wb));
            count <- 1;
            ## choose desired software
            print("Select gene driver software (the same used for MAF) from");
            for (elem in sheets) {
                if(elem == "IMCDriver"){
                    cat(elem, "(suggested)", count, "\n");
                }else{
                    cat(elem, count, "\n");
                }
                count <- count + 1
            }
            
            choice_raw <- readline(prompt = "Enter your choice: ");
            choice <- as.integer(choice_raw);
            choice_str <- as.character(choice);
            
            while(choice_raw != choice_str ||
                  choice <= 0 ||
                  choice > length(sheets)){
                      choice_raw <- readline(prompt = "Enter your choice: "); 
                      choice <- as.integer(choice_raw);
                      choice_str <- as.character(choice);
                  }
            ## load genes
            genes <- read.xlsx(file_drivers, 
                               sheetIndex = choice, 
                               header = TRUE);
            
            lastAnalyseDate <- getFirehoseAnalyzeDates(1);
            
            gistic <- getFirehoseData("LUAD",
                                      gistic2_Date = lastAnalyseDate, 
                                      GISTIC = TRUE);
            ##gistic <- as.data.frame(gistic);
            ## get GISTIC results
            LUADGistic <- getData(gistic,
                                  type = "GISTIC", 
                                  platform = "AllByGene");
            ## gistic.thresholedbygene <- getData(gistic,
            ##                                    type = "GISTIC", 
            ##                                    platform = "ThresholdedByGene");
            
            LUADGistic <- as.data.frame(LUADGistic);
            
            LUAD <- import.GISTIC(LUADGistic,
                                  filter.genes = genes$LUAD);
            
            save(LUAD, file = "input/luadGisticgdc.rda");
            return(LUAD);
        }else{
            ## default genes already in a RDATA file
            LUAD <- loadRData("input/luadGisticgdc.rda");
            return (LUAD);
        }
    }else{
        ## default genes already in a RDATA file
        LUAD <- loadRData("input/luadGisticgdc.rda");
        return(LUAD);
    }
}


reduce_samples <- function(LUAD, n_samples){
    uni_samples <- unique(LUAD$Tumor_Sample_Barcode);
    if(length(uni_samples) < n_samples){
        n_samples <- length(uni_samples);
    }
    
    samples <- uni_samples[1:n_samples];
    LUADreduced <- subset(LUAD, Tumor_Sample_Barcode %in% samples);
    return(LUADreduced);
}
