# SUBTYPING

### Molecular Subtypes
subtypes <- TCGAquery_subtype(tumor = "luad")

subtypes.TRU <- subtypes[subtypes$expression_subtype %in% "TRU",]
subtypes.PP <- subtypes[subtypes$expression_subtype %in% "prox.-prolif.",]
subtypes.PI <- subtypes[subtypes$expression_subtype %in% "prox.-inflam",]

LUAD.TRU <- trim(samples.selection(LUAD, subtypes.TRU$patient))
LUAD.PP <- trim(samples.selection(LUAD, subtypes.PP$patient))
LUAD.PI <- trim(samples.selection(LUAD, subtypes.PI$patient))


if (histological_verbose) {
  ## Histologic Types
  ## reload clinical data
  clinical_sub <- read_tsv(file.clinical)

  ## names of interesting subtypes
  acinar <- "lung acinar adenocarcinoma"
  nonmucinous <- "lung bronchioloalveolar carcinoma nonmucinous"
  papillary <- "lung papillary adenocarcinoma"
  mucinous <- "mucinous (colloid) carcinoma"

  ## detect samples
  acinar_samples <-
    clinical_sub[clinical_sub$histological_type == acinar,]
  acinar_samples <-
    acinar_samples[acinar_samples$rn %in% as.samples(LUAD),]
  acinar_samples <- acinar_samples$rn

  nonmucinous_samples <-
    clinical_sub[clinical_sub$histological_type == nonmucinous,]
  nonmucinous_samples <-
    nonmucinous_samples[nonmucinous_samples$rn %in% as.samples(LUAD),]
  nonmucinous_samples <- nonmucinous_samples$rn

  ## detect samples
  papillary_samples <-
    clinical_sub[clinical_sub$histological_type == papillary,]
  papillary_samples <-
    papillary_samples[papillary_samples$rn %in% as.samples(LUAD),]
  papillary_samples <- papillary_samples$rn

  mucinous_samples <-
    clinical_sub[clinical_sub$histological_type == mucinous,]
  mucinous_samples <-
    mucinous_samples[mucinous_samples$rn %in% as.samples(LUAD),]
  mucinous_samples <- mucinous_samples$rn


  ## subset of LUAD tronco object for every subtype
  LUAD.acinar <- trim(samples.selection(LUAD, acinar_samples))
  LUAD.acinar <- annotate.description(LUAD.acinar,
                                      "LUAD acinar subtype")

  LUAD.nonmucinous <-
    trim(samples.selection(LUAD, nonmucinous_samples))
  LUAD.nonmucinous <- annotate.description(LUAD.nonmucinous,
                                           "LUAD nonmucinous subtype")

  LUAD.papillary <- trim(samples.selection(LUAD, papillary_samples))
  LUAD.papillary <- annotate.description(LUAD.papillary,
                                         "LUAD papillary subtype")

  LUAD.mucinous <- trim(samples.selection(LUAD, mucinous_samples))
  LUAD.mucinous <- annotate.description(LUAD.mucinous,
                                        "LUAD mucinous subtype")

}

old_events <- nevents(LUAD)

## select events for complete analysis (alla data) with a min freq
LUAD <- events.selection(LUAD,
                         filter.freq = min_freq)
new_events <- nevents(LUAD)

LUAD.TRU <- events.selection(LUAD.TRU,
                             filter.freq = min_freq)
LUAD.PP <- events.selection(LUAD.PP,
                            filter.freq = min_freq)
LUAD.PI <- events.selection(LUAD.PI,
                            filter.freq = min_freq)

## some prints about molecular subtypes
if (verbose) {
  print(paste(
    "TRU subtype has",
    nsamples(LUAD.TRU),
    "samples and",
    nevents(LUAD.TRU),
    "events"
  ))
  print(paste(
    "PP subtype has",
    nsamples(LUAD.PP),
    "samples and",
    nevents(LUAD.PP),
    "events"
  ))
  print(paste(
    "PI subtype has",
    nsamples(LUAD.PI),
    "samples and",
    nevents(LUAD.PI),
    "events"
  )) }

LUAD.TRU <- annotate.description(LUAD.TRU,
                                 "LUAD TRU subtype")
LUAD.PP <- annotate.description(LUAD.PP,
                                "LUAD PP subtype")
LUAD.PI <- annotate.description(LUAD.PI,
                                "LUAD PI subtype")
LUAD <- annotate.description(
  LUAD,
  "LUAD somatic mutations and CNA from GCD portal, final")


if (plot_verbose) {
  oncoprint(
    LUAD.TRU,
    gene.annot = pathway.list,
    gene.annot.color = pathways.color,
  )
}

if (plot_verbose) {
  oncoprint(
    LUAD.PP,
    gene.annot = pathway.list,
    gene.annot.color = pathways.color,
  )
}

if (plot_verbose) {
  oncoprint(
    LUAD.PI,
    gene.annot = pathway.list,
    gene.annot.color = pathways.color,
  )
}

if (plot_verbose) {
  oncoprint(
    LUAD,
    gene.annot = pathway.list,
    gene.annot.color = pathways.color,
  )
}

if (histological_verbose) {
  ## select events for every subtype and get the corresponding TRONCO objects
  LUAD.acinar <- events.selection(LUAD.acinar, filter.freq = min_freq)
  LUAD.nonmucinous <- events.selection(LUAD.nonmucinous, filter.freq = min_freq)
  LUAD.papillary <- events.selection(LUAD.papillary, filter.freq = min_freq)
  LUAD.mucinous <- events.selection(LUAD.mucinous, filter.freq = min_freq)

  ## print number of samples and events
  if (verbose) {
    print(paste(
      "Acinar subtype has",
      nsamples(LUAD.acinar),
      "events and",
      nevents(LUAD.acinar),
      "events"
    ))
    print(paste(
      "Nonmucinous subtype has",
      nsamples(LUAD.nonmucinous),
      "events and",
      nevents(LUAD.nonmucinous),
      "events"
    ))
    print(paste(
      "Papillary subtype has",
      nsamples(LUAD.papillary),
      "events and",
      nevents(LUAD.papillary),
      "events"
    ))
    print(paste(
      "Mucinous subtype has",
      nsamples(LUAD.mucinous),
      "events and",
      nevents(LUAD.mucinous),
      "events"
    ))
  }

  ## oncoprint of the subtypes
  if (plot_verbose) {
    oncoprint(LUAD.acinar)
  }
  if (plot_verbose) {
    oncoprint(
      LUAD.acinar,
      gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      title = ""
    )
  }
  if (plot_verbose) {
    oncoprint(LUAD.nonmucinous)
  }
  if (plot_verbose) {
    oncoprint(
      LUAD.nonmucinous,
      gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      title = ""
    )

  }
  if (plot_verbose) {
    oncoprint(LUAD.papillary)
  }
  if (plot_verbose) {
    oncoprint(
      LUAD.papillary,
      gene.annot = pathway.list,
      gene.annot.color = pathways.color,
      title = ""
    )

    # oncoprint(LUAD.mucinous)
    # oncoprint(LUAD.mucinous,
    #           legend.cex = .5,
    #           cellwidth = 3,
    #           cellheight = 10,
    #           gene.annot = pathway.list,
    #           gene.annot.color = pathways.color,
    #           title = "")
  }
}