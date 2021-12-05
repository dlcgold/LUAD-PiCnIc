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

  # use custom genes
  if(length(genes_custom) != 0){
    luad_data <- cbio.query(
      genes = genes_custom[!is.na(genes_custom)],
      cbio.study = 'luad_tcga_pub',
      cbio.dataset = 'luad_tcga_pub_cnaseq',
      cbio.profile = 'luad_tcga_pub_mutations');
    
    m <- as.matrix(luad_data$profile);
    m[is.na(m)] <- 0;
    m[m == 'NaN'] = 0;
    m[m != '0'] <- 1
    LUAD <- import.genotypes(m, 
                             event.type = 'Mutation', 
                             color = 'brown3');
    
    LUAD <- annotate.description(LUAD, 
                                 "Lung cancer data from Cbio portal");
    return(LUAD);
  }
  
  if(interactive){
    # choose if use default drivers
    choice_def <- readline(prompt = "Use default gene drivers (IMCDriver)? [yes/no] ");
    choice_def <- as.character(choice_def);
    while(choice_def != "yes" && choice_def != "no"){
      choice_def <- readline(prompt = "Use default gene drivers? [yes/no] ");
      choice_def <- as.character(choice_def);
    }
    if (choice_def == "yes"){
      choice_bool <- TRUE;
    }else{
      choice_bool <- FALSE;
    }
    
    if(!choice_bool){
      # if not default load xlsx file from suppl. 4
      file_drivers <- "input/gene_drivers.xlsx"
      wb <- loadWorkbook(file_drivers)
      sheets <- names(getSheets(wb))
      count <- 1
      # choose desired software
      print("Select gene driver software from")
      for (elem in sheets) {
        if(elem == "IMCDriver"){
          cat(elem, "(suggested)", count, "\n")
        }else{
          cat(elem, count, "\n")
        }
        count <- count + 1
      }
      
      choice_raw <- readline(prompt = "Enter your choice: ");
      choice <- as.integer(choice_raw);
      choice_str <- as.character(choice);
      
      while(choice_raw != choice_str || choice <= 0 || choice > length(sheets)){
        choice_raw <- readline(); 
        choice <- as.integer(choice_raw);
        choice_str <- as.character(choice);
      }
      # load genes
      genes <- read.xlsx(file_drivers, 
                         sheetIndex = choice, 
                         header = TRUE);
      # download from cBIO
      luad_data <- cbio.query(
        genes = genes$LUAD[!is.na(genes$LUAD)],
        cbio.study = 'luad_tcga_pub',
        cbio.dataset = 'luad_tcga_pub_cnaseq',
        cbio.profile = 'luad_tcga_pub_mutations');
      
      # prepare data for TRONCO usage
      m <- as.matrix(luad_data$profile);
      m[is.na(m)] <- 0;
      m[m == 'NaN'] = 0;
      m[m != '0'] <- 1
      LUAD <- import.genotypes(m, 
                               event.type = 'Mutation', 
                               color = 'brown3');
      
      LUAD <- annotate.description(LUAD, 
                                   "Lung cancer data from Cbio portal");
      return(LUAD);
    }else{
      # default genes already in a RDATA file
      LUAD <- loadRData("input/luadIMC.rda");
      return (LUAD);
    }
  }else{
    # default genes already in a RDATA file
    LUAD <- loadRData("input/luadIMC.rda");
    return(LUAD);
  }
}