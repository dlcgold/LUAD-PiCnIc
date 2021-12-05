library("xlsx")
library(TRONCO)

loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}


#' input_genes
#'
#' function to load desired data, default are IMCDriver genes driver
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
input_genes <- function(interactive = TRUE, genes_custom = c()) {

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
            
            luad_data <- cbio.query(
                genes = genes$LUAD[!is.na(genes$LUAD)],
                cbio.study = 'luad_tcga',
                cbio.dataset = 'luad_tcga_3way_complete',
                cbio.profile = 'luad_tcga_gistic');
            
            m <- as.matrix(luad_data$profile);
            m[is.na(m)] <- 0;
            m[m == 'NaN'] = 0;
            m[m != '0'] <- 1
            
            LUAD_cna <- import.genotypes(m, 
                                         event.type = 'CNA', 
                                         color = 'blue3');
            luad_data <- cbio.query(
                genes = genes$LUAD[!is.na(genes$LUAD)],
                cbio.study = 'luad_tcga',
                cbio.dataset = 'luad_tcga_3way_complete',
                cbio.profile = 'luad_tcga_methylation_hm27');
            
            m <- as.matrix(luad_data$profile);
            m[is.na(m)] <- 0;
            m[m == 'NaN'] = 0;
            m[m != '0'] <- 1
            
            LUAD_hm27 <- import.genotypes(m, 
                                          event.type = 'hm27', 
                                          color = 'green3');
            LUAD <- intersect.datasets(LUAD_mut,
                                       LUAD_cna, 
                                       intersect.genomes = FALSE);
            LUAD <- intersect.datasets(LUAD,
                                       LUAD_hm27, 
                                       intersect.genomes = FALSE);
            LUAD <- annotate.description(LUAD, 
                                         "Lung cancer data from Cbio portal");
            #save(LUAD, file = "input/luadIMC2020full.rda");
            return(LUAD);
        }else{
            ## default genes already in a RDATA file
            LUAD <- loadRData("input/luadIMC2020full.rda");
            return (LUAD);
        }
    }else{
        ## default genes already in a RDATA file
        LUAD <- loadRData("input/luadIMC2020full.rda");
        return(LUAD);
    }
}
