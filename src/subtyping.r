# SUBTYPING

## TODO understand if we have cluster (type of tumor like MSI/MSS for CRC)
## and if we can extract them from clinical sheet
## maybe from cBIO
# luad_cbio <- cbio.query(
#   genes = pathway.genes,
#   cbio.study = 'luad_tcga',
#   cbio.dataset = 'luad_tcga_3way_complete',
#   cbio.profile = 'luad_tcga_mutations')

## reload clinical data

clinical_sub <- read_tsv(file.clinical)

## names of interesting subtypes
acinar <- "lung acinar adenocarcinoma"
nonmucinous <- "lung bronchioloalveolar carcinoma nonmucinous"
papillary <- "lung papillary adenocarcinoma"
mucinous <- "mucinous (colloid) carcinoma"

## detect samples
acinar_samples <- clinical_sub[clinical_sub$histological_type == acinar, ]
acinar_samples <- acinar_samples[acinar_samples$rn %in% as.samples(LUAD),]
acinar_samples <- acinar_samples$rn

nonmucinous_samples <- clinical_sub[clinical_sub$histological_type == nonmucinous, ]
nonmucinous_samples <- nonmucinous_samples[nonmucinous_samples$rn %in% as.samples(LUAD),]
nonmucinous_samples <- nonmucinous_samples$rn

## detect samples
papillary_samples <- clinical_sub[clinical_sub$histological_type == papillary, ]
papillary_samples <- papillary_samples[papillary_samples$rn %in% as.samples(LUAD),]
papillary_samples <- papillary_samples$rn

mucinous_samples <- clinical_sub[clinical_sub$histological_type == mucinous, ]
mucinous_samples <- mucinous_samples[mucinous_samples$rn %in% as.samples(LUAD),]
mucinous_samples <- mucinous_samples$rn


## subset of LUAD tronco object
LUAD.acinar <- trim(samples.selection(LUAD, acinar_samples))
LUAD.acinar <- annotate.description(LUAD.acinar, 
                                    "LUAD acinar subtype")

LUAD.nonmucinous <- trim(samples.selection(LUAD, nonmucinous_samples))
LUAD.nonmucinous <- annotate.description(LUAD.nonmucinous, 
                                         "LUAD nonmucinous subtype")

LUAD.papillary <- trim(samples.selection(LUAD, papillary_samples))
LUAD.papillary <- annotate.description(LUAD.papillary, 
                                       "LUAD papillary subtype")

LUAD.mucinous <- trim(samples.selection(LUAD, mucinous_samples))
LUAD.mucinous <- annotate.description(LUAD.mucinous, 
                                      "LUAD mucinous subtype")

old_events <- nevents(LUAD)
## select events for complete analysis with a min freq
LUAD <- events.selection(LUAD, 
                         filter.freq = min_freq)
new_events <- nevents(LUAD)
if(plot_verbose){dev.new()
  oncoprint(LUAD)
}

LUAD.acinar <- events.selection(LUAD.acinar, 
                                filter.freq = min_freq)
LUAD.nonmucinous <- events.selection(LUAD.nonmucinous, 
                                     filter.freq = min_freq)
LUAD.papillary <- events.selection(LUAD.papillary, 
                                   filter.freq = min_freq)
LUAD.mucinous <- events.selection(LUAD.mucinous, 
                                  filter.freq = min_freq)

## print number of samples and events
if(verbose){
  print(paste("Acinar subtype has",
              nsamples(LUAD.acinar), 
              "events and", 
              nevents(LUAD.acinar), 
              "events"))
  print(paste("Nonmucinous subtype has",
              nsamples(LUAD.nonmucinous), 
              "events and", 
              nevents(LUAD.nonmucinous), 
              "events"))
  print(paste("Papillary subtype has",
              nsamples(LUAD.papillary), 
              "events and", 
              nevents(LUAD.papillary), 
              "events"))
  print(paste("Mucinous subtype has",
              nsamples(LUAD.mucinous), 
              "events and", 
              nevents(LUAD.mucinous), 
              "events"))
}

## oncoprint of the subtypes
if(plot_verbose){dev.new()
  oncoprint(LUAD.acinar)
  oncoprint(LUAD.acinar,
            legend.cex = .5,          
            cellwidth = 3,           
            cellheight = 10,
            gene.annot = pathway.list, 
            gene.annot.color = pathways.color,
            title = "")
  
  oncoprint(LUAD.nonmucinous)
  oncoprint(LUAD.nonmucinous,
            legend.cex = .5,          
            cellwidth = 3,           
            cellheight = 10,
            gene.annot = pathway.list, 
            gene.annot.color = pathways.color,
            title = "")
  
  oncoprint(LUAD.papillary)
  oncoprint(LUAD.papillary,
            legend.cex = .5,          
            cellwidth = 3,           
            cellheight = 10,
            gene.annot = pathway.list, 
            gene.annot.color = pathways.color,
            title = "")
  
  # oncoprint(LUAD.mucinous)
  # oncoprint(LUAD.mucinous,
  #           legend.cex = .5,          
  #           cellwidth = 3,           
  #           cellheight = 10,
  #           gene.annot = pathway.list, 
  #           gene.annot.color = pathways.color,
  #           title = "")
}