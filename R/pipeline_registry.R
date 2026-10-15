# Describe stages, dependencies and expected principal outputs
# These declarations allow a dry run without loading the research data.
ccrvi_registry <- function(cfg) {
  cohort <- cfg$files$compiled_cohort
  original <- cfg$files$original_results
  no_aspirin <- cfg$files$no_aspirin_results
  candidate <- "outputs/models/lead_candidate_no_aspirin/05_no_aspirin_candidate_outputs.rds"
  core <- c("tidyverse", "caret", "glmnet", "pROC", "PRROC")
  stage <- function(script, inputs, outputs, packages, description) {
    list(script = script, inputs = inputs, outputs = outputs,
         packages = packages, description = description)
  }
  list(
    compile = stage("scripts/01_compile_cohort.R", unlist(cfg$files[c("clinical", "outcome", "workbook")]),
                    c(cohort, cfg$files$unfiltered_cohort), c("tidyverse", "readxl"), "Compile the complete-intersection cohort"),
    eda = stage("scripts/02_exploratory_analysis.R", cohort, "outputs/tables/eda_cohort_overview.csv",
                c("tidyverse", "janitor", "broom", "pROC", "pheatmap", "recipes", "factoextra", "naniar", "patchwork"),
                "Describe the cohort, missingness and molecular structure"),
    endpoint = stage("scripts/02a_endpoint_characterisation.R", c(cohort, cfg$files$clinical),
                     "outputs/models/mace_clinical_event_characterisation_outputs.rds", "tidyverse", "Characterise recorded MACE events"),
    original = stage("scripts/03_compare_modalities.R", cohort, original, core, "Original exhaustive elastic-net comparisons"),
    `no-aspirin` = stage("scripts/03b_compare_no_aspirin.R", c(cohort, original), no_aspirin, core,
                        "Compare aspirin-excluded expanded-clinical combinations"),
    clinical = stage("scripts/03d_clinical_interrogation.R", cohort, "outputs/models/clinical_interrogation_results.rds", core,
                     "Supporting clinical benchmark and ablation analyses"),
    medications = stage("scripts/03e_medication_investigation.R", cohort, "outputs/models/03E_medication_confounding_investigation_outputs.rds",
                        c("tidyverse", "broom"), "Supporting medication and phenotype investigation"),
    sensitivity = stage("scripts/03f_medication_sensitivity.R", cohort,
                        "outputs/models/medication_sensitivity_script03/03F_medication_sensitivity_script03_outputs.rds", core,
                        "Refit the focused medication sensitivity panel"),
    selection = stage("scripts/04_select_models.R", no_aspirin,
                      "outputs/models/model_selection_no_aspirin/04_no_aspirin_model_selection_outputs.rds",
                      c("tidyverse", "pROC", "PRROC"), "Rank models and calculate paired AUC comparisons"),
    rf = stage("scripts/04b_random_forest.R", cohort,
               "outputs/models/random_forest_complete_modality_no_aspirin_nested_cv_outputs.rds",
               c("tidyverse", "caret", "pROC", "PRROC", "ranger"), "Fit the optional complete-modality random forest"),
    candidate = stage("scripts/05_characterise_candidate.R", c(cohort, no_aspirin), candidate,
                      c("tidyverse", "glmnet", "pROC", "PRROC"), "Characterise the selected candidate; includes an interpretation refit"),
    risk = stage("scripts/06_exploratory_risk_classes.R", candidate,
                 "outputs/models/risk_score/06_mace_risk_score_outputs.rds", c("tidyverse", "pROC"), "Derive exploratory percentile scores and classes"),
    masld = stage("scripts/07_masld_analysis.R", c(cohort, candidate),
                  "outputs/models/masld_stratification/07_masld_stratification_outputs.rds", c(core, "broom"),
                  "Evaluate MASLD subgroups; exploratory retraining is controlled separately"),
    biomarkers = stage("scripts/08_biomarker_recurrence.R", c(original, no_aspirin),
                       "outputs/models/biomarker_recurrence/08_cross_model_biomarker_recurrence_outputs.rds", "tidyverse", "Summarise cross-model feature recurrence"),
    synthesis = stage("scripts/09_final_synthesis.R", c(cohort, original, no_aspirin),
                      c("outputs/models/final_thesis_synthesis/09_final_model_and_biomarker_synthesis_outputs.rds",
                        "outputs/tables/final_thesis_synthesis/09B_top5_auc_comparison_vs_baseline.csv"),
                      c("tidyverse", "pROC", "caret"), "Assemble final performance and biomarker summaries"),
    slides = stage("scripts/09b_slide_figure.R", "outputs/tables/final_thesis_synthesis/09B_top5_auc_comparison_vs_baseline.csv",
                   character(), "tidyverse", "Export the presentation figure; existing filenames may be suffixed")
  )
}

ccrvi_profiles <- list(
  main = c("compile", "original", "no-aspirin", "selection", "candidate", "risk", "masld", "biomarkers", "synthesis"),
  descriptive = c("eda", "endpoint"),
  sensitivity = c("clinical", "medications", "sensitivity"),
  rf = "rf",
  `masld-exploratory` = "masld",
  results = c("selection", "candidate", "risk", "masld", "biomarkers", "synthesis"),
  figures = c("synthesis", "slides"),
  all = c("compile", "eda", "endpoint", "original", "no-aspirin", "clinical", "medications", "sensitivity",
          "selection", "rf", "candidate", "risk", "masld", "biomarkers", "synthesis", "slides")
)
