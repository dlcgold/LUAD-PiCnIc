library(TRONCO)

stagesStudy <- function(LUAD, LUAD_samples, n_stages){
    if(n_stages >= LUAD_samples){
        stop('more stages than samples')
    }
    patient_for_stage = as.integer(LUAD_samples / n_stages);
    print(patient_for_stage);
    stages = c();
    for (i in 1:n_stages) {
        stage = paste("stage", as.character(i), collapse = "");
        print(stage);
        if(i != n_stages){
            stages <- c(stages, rep(stage, patient_for_stage));
        }else{
            fill = LUAD_samples - (patient_for_stage * n_stages);
            stages <- c(stages, rep(stage, patient_for_stage + fill));
        }
    }
    stages = as.matrix(stages);
    rownames(stages) = as.samples(LUAD);
    LUAD_grouped = annotate.stages(LUAD, stages = stages);
    print(head(as.stages(LUAD_grouped)));
    oncoprint(LUAD_grouped);
    oncoprint(LUAD_grouped, 
              group.samples = as.stages(LUAD_grouped));
    oncoprint(LUAD_grouped,
              samples.cluster = TRUE,
              genes.cluster = TRUE);
}
