# Repository stage: synthesis
# Source supplied as: 09_final_model_and_biomarker_synthesis_v3(1).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# 09 Final model performance and biomarker synthesis
#
# Project:
#   Cardiovascular Risk Stratification in Liver Cirrhosis
#
# Purpose:
#   Produce the final thesis-facing comparison of model discrimination and
#   biomarker stability without retraining any models.
#
#   This script combines the two completed nested cross-validation analyses:
#
#   1. The original complete multimodality analysis, which contains the compact
#      baseline clinical model and the six individual omics/pathology models.
#   2. The aspirin-excluded expanded-clinical analysis, which contains the
#      locked candidate, the all-modality model and the final integrated model
#      ranking used for model selection.
#
#   The principal clinical benchmark is the compact baseline clinical model.
#   All AUCs are calculated from pooled outer-fold out-of-sample predictions.
#   Paired DeLong tests are performed only after confirming that predictions
#   from the two Script 03 runs can be mapped back to the same held-out cohort
#   rows and validation folds.
#
# Main thesis outputs:
#   09A_prespecified_auc_comparison_vs_baseline.csv/png
#     Baseline clinical + all six single-modality models + locked candidate +
#     all-modality model.
#
#   09B_top5_auc_comparison_vs_baseline.csv/png
#     Top 5 aspirin-excluded integrated models, each compared with the compact
#     baseline clinical benchmark.
#
#   09C_monomodality_delong_vs_baseline.csv
#     Paired DeLong tests for each single-modality model versus baseline.
#
#   09D_top5_delong_vs_baseline.csv
#     Paired DeLong tests for the top 5 final models versus baseline.
#
#   09E_top5_recurrent_stable_biomarkers.csv
#   09F_candidate_stable_biomarkers.csv
#     Biomarker recurrence across the empirically top-ranked final models and
#     stable predictors in the locked candidate.
#
#   09G_cross_model_recurrent_biomarker_summary.csv
#     Optional broader recurrence summary extracted from Script 08 when available.
#
# Statistical interpretation:
#   DeLong testing evaluates whether two correlated ROC AUCs differ when they
#   are estimated in the same participants. Raw and Benjamini-Hochberg adjusted
#   p-values are both retained. The six single-modality comparisons are a
#   prespecified comparison family. DeLong tests applied to the empirically top
#   10 models are explicitly exploratory because those same out-of-fold AUCs
#   were used to select the top 5. These tests do not replace assessment of
#   effect size, confidence intervals, calibration, or external validation.

library(tidyverse)
library(pROC)
library(caret)

set.seed(123)

# Output folders
output_table_dir <- "outputs/tables/final_thesis_synthesis"
output_figure_dir <- "outputs/figures/final_thesis_synthesis"
output_model_dir <- "outputs/models/final_thesis_synthesis"

dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_model_dir, recursive = TRUE, showWarnings = FALSE)

# Core settings
baseline_model <- "clinical_baseline"

single_modality_models <- c(
  "lipoproteins",
  "lipids",
  "small_molecules",
  "cytokines",
  "bile_acids",
  "pathology"
)

candidate_model <- paste0(
  "clinical_expanded_no_aspirin_",
  "lipoproteins_lipids_small_molecules"
)

all_modalities_model <- paste0(
  "clinical_expanded_no_aspirin_",
  "lipoproteins_lipids_small_molecules_",
  "cytokines_bile_acids_pathology"
)

top_n_models <- 5
stable_selection_cutoff <- 0.60
expanded_clinical_model <- "clinical_expanded"

# Thesis figure colours
# The palette deliberately uses the orange/blue visual language used in the
# current thesis figures. Shades of blue distinguish model classes while the
# clinical baseline remains visually distinct as the reference.
clinical_orange <- "#E69F00"
omics_blue <- "#56B4E9"
candidate_blue <- "#0072B2"
all_modalities_blue <- "#004C6D"

model_group_colours <- c(
  "Clinical baseline" = clinical_orange,
  "Single modality" = omics_blue,
  "Locked candidate" = candidate_blue,
  "All modalities" = all_modalities_blue,
  "Expanded clinical" = "#F3B34C",
  "Top 5 model" = candidate_blue
)

# Input files
original_output_candidates <- c(ccrvi_config$files$original_results)

no_aspirin_output_candidates <- c(ccrvi_config$files$no_aspirin_results)

biomarker_output_candidates <- c(
  "outputs/models/biomarker_recurrence/08_cross_model_biomarker_recurrence_outputs.rds"
)

compiled_cohort_file <- ccrvi_config$files$compiled_cohort

# Save helpers
save_table <- function(x, filename) {
  write.csv(
    x,
    file.path(output_table_dir, filename),
    row.names = FALSE
  )
}

save_plot <- function(plot_object, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(output_figure_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

# File and label helpers
find_first_existing_file <- function(candidates, description) {
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    stop(
      description,
      " was not found. Expected one of:\n",
      paste(candidates, collapse = "\n")
    )
  }

  existing[1]
}

make_model_label <- function(model_name) {
  case_when(
    model_name == baseline_model ~ "Baseline clinical",
    model_name == expanded_clinical_model ~ "Expanded clinical",
    model_name == "lipoproteins" ~ "Lipoproteins only",
    model_name == "lipids" ~ "Lipids only",
    model_name == "small_molecules" ~ "Small molecules only",
    model_name == "cytokines" ~ "Cytokines only",
    model_name == "bile_acids" ~ "Bile acids only",
    model_name == "pathology" ~ "Pathology only",
    model_name == candidate_model ~ "Locked candidate",
    model_name == all_modalities_model ~ "All modalities",
    TRUE ~ model_name %>%
      str_replace(
        "^clinical_expanded_no_aspirin_",
        "Expanded clinical (aspirin excluded) + "
      ) %>%
      str_replace(
        "^clinical_expanded_no_aspirin$",
        "Expanded clinical (aspirin excluded)"
      ) %>%
      str_replace("^clinical_expanded_", "Expanded clinical + ") %>%
      str_replace("^clinical_baseline_", "Baseline clinical + ") %>%
      str_replace_all("lipoproteins", "Lipoproteins") %>%
      str_replace_all("small_molecules", "Small molecules") %>%
      str_replace_all("bile_acids", "Bile acids") %>%
      str_replace_all("lipids", "Lipids") %>%
      str_replace_all("cytokines", "Cytokines") %>%
      str_replace_all("pathology", "Pathology") %>%
      str_replace_all("_", " + ") %>%
      str_squish()
  )
}

format_p_value <- function(p) {
  case_when(
    is.na(p) ~ "p = NA",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", formatC(p, format = "f", digits = 3))
  )
}

format_delta_auc <- function(x) {
  ifelse(
    is.na(x),
    "",
    paste0(ifelse(x >= 0, "+", ""), formatC(x, format = "f", digits = 3))
  )
}

classify_feature_modality <- function(feature_name) {
  case_when(
    str_starts(feature_name, "lipo_") ~ "Lipoproteins",
    str_starts(feature_name, "lipid_") ~ "Lipids",
    str_starts(feature_name, "sm_") ~ "Small molecules",
    str_starts(feature_name, "cyto_") ~ "Cytokines",
    str_starts(feature_name, "ba_") ~ "Bile acids",
    str_starts(feature_name, "path_") ~ "Pathology",
    TRUE ~ "Clinical"
  )
}

clean_feature_label <- function(feature_name) {
  feature_name %>%
    str_remove("^(lipo_|lipid_|sm_|cyto_|ba_|path_)") %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

# Consistent thesis plot theme
thesis_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, colour = "#1F2A44"),
    plot.subtitle = element_text(size = 11.5, colour = "#2F3A4A"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12.5, colour = "#1F2A44"),
    axis.text = element_text(colour = "#4A4A4A"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#E6E6E6", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10.5),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(14, 72, 14, 18)
  )

# Load completed model outputs
original_output_file <- find_first_existing_file(
  original_output_candidates,
  "Original Script 03 output"
)

no_aspirin_output_file <- find_first_existing_file(
  no_aspirin_output_candidates,
  "Aspirin-excluded Script 03 output"
)

if (!file.exists(compiled_cohort_file)) {
  stop(
    "Compiled cohort file was not found: ",
    compiled_cohort_file,
    "\nThis file is required to reconstruct the outer-fold participant mapping ",
    "for cross-run paired DeLong comparisons."
  )
}

original_outputs <- readRDS(original_output_file)
no_aspirin_outputs <- readRDS(no_aspirin_output_file)

required_output_objects <- c(
  "modality_predictions",
  "pooled_modality_performance",
  "modality_feature_votes"
)

for (object_name in required_output_objects) {
  if (!object_name %in% names(original_outputs)) {
    stop("Original Script 03 output is missing: ", object_name)
  }

  if (!object_name %in% names(no_aspirin_outputs)) {
    stop("Aspirin-excluded Script 03 output is missing: ", object_name)
  }
}

# Reconstruct the exact global outer-fold participant mapping
#
# Script 03 does not save participant IDs in the prediction table. Both Script 03
# versions, however, create the same five global validation folds from the same
# complete cohort with set.seed(123) before any model fitting. Recreating those
# folds here allows each out-of-fold prediction to be mapped back to its original
# cohort row. This provides a defensible participant-level key for paired DeLong
# tests across the two separately saved model runs.
cohort_for_pairing <- read.csv(compiled_cohort_file) %>%
  mutate(mace = as.numeric(mace)) %>%
  filter(!is.na(mace))

set.seed(123)
recreated_validation_folds <- caret::createFolds(
  factor(cohort_for_pairing$mace),
  k = 5,
  returnTrain = FALSE
)

# Retain a participant identifier in the audit table when one is available.
possible_id_columns <- c(
  "participant_id",
  "participantid",
  "study_id",
  "studyid",
  "subject_id",
  "subjectid",
  "record_id",
  "recordid",
  "id"
)

participant_id_column <- possible_id_columns[
  possible_id_columns %in% names(cohort_for_pairing)
][1]

fold_key <- map2_dfr(
  recreated_validation_folds,
  seq_along(recreated_validation_folds),
  function(cohort_rows, validation_fold) {
    output <- tibble(
      validation_fold = as.integer(validation_fold),
      row_in_fold = seq_along(cohort_rows),
      cohort_row = cohort_rows,
      expected_outcome = cohort_for_pairing$mace[cohort_rows]
    )

    if (!is.na(participant_id_column)) {
      output <- output %>%
        mutate(
          participant_id = as.character(
            cohort_for_pairing[[participant_id_column]][cohort_rows]
          )
        )
    }

    output
  }
)

# Attach the reconstructed participant key to each saved prediction table.
prepare_prediction_table <- function(prediction_table, source_name) {
  required_columns <- c("model", "validation_fold", "outcome", "pred")
  missing_columns <- setdiff(required_columns, names(prediction_table))

  if (length(missing_columns) > 0) {
    stop(
      "Prediction table for ", source_name, " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  prepared <- prediction_table %>%
    transmute(
      source = source_name,
      model = as.character(model),
      validation_fold = as.integer(validation_fold),
      outcome = as.numeric(outcome),
      pred = as.numeric(pred)
    ) %>%
    group_by(source, model, validation_fold) %>%
    mutate(row_in_fold = row_number()) %>%
    ungroup() %>%
    left_join(
      fold_key,
      by = c("validation_fold", "row_in_fold")
    )

  unmatched_rows <- prepared %>%
    filter(is.na(cohort_row))

  if (nrow(unmatched_rows) > 0) {
    stop(
      "Could not map all ", source_name,
      " predictions to the reconstructed validation folds."
    )
  }

  outcome_mismatches <- prepared %>%
    filter(outcome != expected_outcome)

  if (nrow(outcome_mismatches) > 0) {
    stop(
      "Outcome mismatch detected while reconstructing participant pairing for ",
      source_name,
      ". Cross-run paired DeLong testing has been stopped rather than assuming ",
      "that prediction rows are aligned."
    )
  }

  prepared
}

original_predictions <- prepare_prediction_table(
  original_outputs$modality_predictions,
  "original"
)

no_aspirin_predictions <- prepare_prediction_table(
  no_aspirin_outputs$modality_predictions,
  "no_aspirin"
)

all_predictions <- bind_rows(
  original_predictions,
  no_aspirin_predictions
)

# Save a compact participant-pairing audit for reproducibility.
pairing_audit <- all_predictions %>%
  group_by(source, model) %>%
  summarise(
    n_predictions = n(),
    n_unique_cohort_rows = n_distinct(cohort_row),
    n_events = sum(outcome == 1, na.rm = TRUE),
    n_nonevents = sum(outcome == 0, na.rm = TRUE),
    outcome_mapping_matches = all(outcome == expected_outcome),
    .groups = "drop"
  )

save_table(pairing_audit, "09_pairing_audit.csv")

# Calculate pooled AUC and 95% DeLong confidence intervals directly from the
# saved out-of-fold predictions. Recalculation here ensures that all models in
# the final synthesis use exactly the same AUC definition and direction.
calculate_auc_table <- function(prediction_table) {
  prediction_table %>%
    group_by(source, model) %>%
    group_modify(function(model_df, model_key) {
      roc_object <- pROC::roc(
        response = model_df$outcome,
        predictor = model_df$pred,
        direction = "<",
        quiet = TRUE
      )

      auc_ci <- pROC::ci.auc(
        roc_object,
        method = "delong"
      )

      tibble(
        pooled_auc = as.numeric(pROC::auc(roc_object)),
        auc_ci_low = as.numeric(auc_ci[1]),
        auc_ci_high = as.numeric(auc_ci[3]),
        n = nrow(model_df),
        n_events = sum(model_df$outcome == 1, na.rm = TRUE),
        n_nonevents = sum(model_df$outcome == 0, na.rm = TRUE)
      )
    }) %>%
    ungroup()
}

all_auc_performance <- calculate_auc_table(all_predictions)

# Confirm all models required for the final thesis comparison are available.
required_original_models <- c(
  baseline_model,
  expanded_clinical_model,
  single_modality_models
)

required_no_aspirin_models <- c(
  candidate_model,
  all_modalities_model
)

missing_original_models <- setdiff(
  required_original_models,
  original_predictions$model
)

missing_no_aspirin_models <- setdiff(
  required_no_aspirin_models,
  no_aspirin_predictions$model
)

if (length(missing_original_models) > 0) {
  stop(
    "Required model(s) missing from original Script 03 output: ",
    paste(missing_original_models, collapse = ", ")
  )
}

if (length(missing_no_aspirin_models) > 0) {
  stop(
    "Required model(s) missing from aspirin-excluded Script 03 output: ",
    paste(missing_no_aspirin_models, collapse = ", ")
  )
}

# Paired DeLong comparison against the compact baseline clinical model
#
# Comparisons may span the two Script 03 runs. Predictions are joined by the
# reconstructed original cohort row rather than by row position alone.
paired_delong_vs_baseline <- function(comparator_source, comparator_model) {
  benchmark_df <- all_predictions %>%
    filter(
      source == "original",
      model == baseline_model
    ) %>%
    select(
      cohort_row,
      outcome_baseline = outcome,
      pred_baseline = pred
    )

  comparator_df <- all_predictions %>%
    filter(
      source == comparator_source,
      model == comparator_model
    ) %>%
    select(
      cohort_row,
      outcome_comparator = outcome,
      pred_comparator = pred
    )

  paired_df <- inner_join(
    comparator_df,
    benchmark_df,
    by = "cohort_row"
  ) %>%
    filter(
      !is.na(pred_comparator),
      !is.na(pred_baseline)
    )

  if (nrow(paired_df) != nrow(comparator_df) ||
      nrow(paired_df) != nrow(benchmark_df)) {
    stop(
      "Incomplete participant pairing for model: ",
      comparator_model,
      ". Paired DeLong testing has been stopped."
    )
  }

  if (!all(paired_df$outcome_comparator == paired_df$outcome_baseline)) {
    stop(
      "Outcome mismatch between comparator and baseline for model: ",
      comparator_model
    )
  }

  roc_comparator <- pROC::roc(
    response = paired_df$outcome_comparator,
    predictor = paired_df$pred_comparator,
    direction = "<",
    quiet = TRUE
  )

  roc_baseline <- pROC::roc(
    response = paired_df$outcome_baseline,
    predictor = paired_df$pred_baseline,
    direction = "<",
    quiet = TRUE
  )

  test_result <- tryCatch(
    pROC::roc.test(
      roc_comparator,
      roc_baseline,
      method = "delong",
      paired = TRUE,
      conf.level = 0.95
    ),
    error = function(e) NULL
  )

  auc_comparator <- as.numeric(pROC::auc(roc_comparator))
  auc_baseline <- as.numeric(pROC::auc(roc_baseline))

  delta_ci <- if (!is.null(test_result) && !is.null(test_result$conf.int)) {
    as.numeric(test_result$conf.int)
  } else {
    c(NA_real_, NA_real_)
  }

  tibble(
    source = comparator_source,
    model = comparator_model,
    model_label = make_model_label(comparator_model),
    n_paired = nrow(paired_df),
    auc_model = auc_comparator,
    auc_baseline = auc_baseline,
    delta_auc_vs_baseline = auc_comparator - auc_baseline,
    delta_auc_ci_low = delta_ci[1],
    delta_auc_ci_high = delta_ci[2],
    delong_p_value = ifelse(
      is.null(test_result),
      NA_real_,
      as.numeric(test_result$p.value)
    )
  )
}

# Baseline AUC is retained as the common reference for all final figures.
baseline_performance <- all_auc_performance %>%
  filter(
    source == "original",
    model == baseline_model
  ) %>%
  slice(1)

if (nrow(baseline_performance) != 1) {
  stop("Baseline clinical performance could not be uniquely identified.")
}

baseline_auc <- baseline_performance$pooled_auc

# Final top-5 model ranking
#
# The final model-development pathway uses the aspirin-excluded expanded
# clinical analysis. The top 5 are therefore ranked within that completed run,
# not by mixing model specifications from different sensitivity analyses.
top_models <- all_auc_performance %>%
  filter(source == "no_aspirin") %>%
  arrange(desc(pooled_auc)) %>%
  slice_head(n = top_n_models) %>%
  mutate(
    rank = row_number(),
    model_label = map_chr(model, make_model_label)
  )

# Biomarker recurrence across the empirically top-ranked final models
#
# This analysis is deliberately descriptive. It asks which predictors repeatedly
# survive penalisation across the top 5 aspirin-excluded models. A biomarker is
# considered stable within a model when selected in at least three of five outer
# folds. Cross-model recurrence is then summarised without creating a composite
# importance score.
no_aspirin_feature_votes <- no_aspirin_outputs$modality_feature_votes %>%
  mutate(
    model = as.character(model),
    feature = as.character(feature),
    feature_modality = map_chr(feature, classify_feature_modality),
    feature_label = map_chr(feature, clean_feature_label),
    stable_3_of_5 = selection_frequency >= stable_selection_cutoff
  )

top5_model_lookup <- top_models %>%
  select(model, rank, model_label)

top5_feature_votes <- no_aspirin_feature_votes %>%
  inner_join(
    top5_model_lookup,
    by = "model"
  )

top5_recurrent_stable_biomarkers <- top5_feature_votes %>%
  group_by(feature, feature_label, feature_modality) %>%
  summarise(
    top_models_selected_n = n_distinct(model),
    top_models_stable_n = n_distinct(model[stable_3_of_5]),
    total_outer_fold_selections = sum(selected_n, na.rm = TRUE),
    mean_model_selection_frequency = mean(selection_frequency, na.rm = TRUE),
    max_model_selection_frequency = max(selection_frequency, na.rm = TRUE),
    mean_coefficient_across_models = mean(mean_coefficient, na.rm = TRUE),
    mean_abs_coefficient_across_models = mean(mean_abs_coefficient, na.rm = TRUE),
    positive_models_n = sum(mean_coefficient > 0, na.rm = TRUE),
    negative_models_n = sum(mean_coefficient < 0, na.rm = TRUE),
    direction_agreement = if_else(
      positive_models_n + negative_models_n > 0,
      pmax(positive_models_n, negative_models_n) /
        (positive_models_n + negative_models_n),
      NA_real_
    ),
    recurrent_direction = case_when(
      positive_models_n > negative_models_n ~ "Increase risk",
      negative_models_n > positive_models_n ~ "Decrease risk",
      TRUE ~ "Mixed or neutral"
    ),
    stable_model_labels = paste(
      model_label[stable_3_of_5][order(rank[stable_3_of_5])],
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  filter(top_models_stable_n >= 2) %>%
  arrange(
    desc(top_models_stable_n),
    desc(top_models_selected_n),
    desc(mean_model_selection_frequency),
    desc(mean_abs_coefficient_across_models)
  ) %>%
  mutate(top5_recurrence_rank = row_number()) %>%
  select(top5_recurrence_rank, everything())

save_table(
  top5_recurrent_stable_biomarkers,
  "09E_top5_recurrent_stable_biomarkers.csv"
)

# Stable biomarkers within the locked candidate are retained as a separate table
# because this is the model carried forward to risk scoring and MASLD subgroup
# analysis. Cross-model top-5 recurrence is attached where available.
candidate_stable_biomarkers_direct <- no_aspirin_feature_votes %>%
  filter(
    model == candidate_model,
    stable_3_of_5
  ) %>%
  left_join(
    top5_recurrent_stable_biomarkers %>%
      select(
        feature,
        top5_recurrence_rank,
        top_models_selected_n,
        top_models_stable_n,
        direction_agreement_top5 = direction_agreement,
        recurrent_direction_top5 = recurrent_direction
      ),
    by = "feature"
  ) %>%
  arrange(
    desc(selection_frequency),
    top5_recurrence_rank,
    desc(mean_abs_coefficient)
  ) %>%
  select(
    feature,
    feature_label,
    feature_modality,
    selected_n,
    selection_frequency,
    mean_coefficient,
    mean_abs_coefficient,
    direction,
    top5_recurrence_rank,
    top_models_selected_n,
    top_models_stable_n,
    direction_agreement_top5,
    recurrent_direction_top5
  )

save_table(
  candidate_stable_biomarkers_direct,
  "09F_candidate_stable_biomarkers.csv"
)

# DeLong tests: all six single-modality models versus baseline clinical.
monomodality_delong <- map_dfr(
  single_modality_models,
  ~ paired_delong_vs_baseline(
    comparator_source = "original",
    comparator_model = .x
  )
) %>%
  mutate(
    comparison_family = "Prespecified single-modality comparison",
    delong_p_adjust_bh = p.adjust(delong_p_value, method = "BH")
  ) %>%
  arrange(desc(delta_auc_vs_baseline))

save_table(
  monomodality_delong,
  "09C_monomodality_delong_vs_baseline.csv"
)

expanded_clinical_delong <- paired_delong_vs_baseline(
  comparator_source = "original",
  comparator_model = expanded_clinical_model
) %>%
  mutate(
    comparison_family = "Prespecified expanded-clinical comparison",
    delong_p_adjust_bh = delong_p_value
  )


# DeLong tests: top 5 final aspirin-excluded models versus baseline clinical.
top5_delong <- top_models %>%
  select(source, model, rank) %>%
  pmap_dfr(function(source, model, rank) {
    paired_delong_vs_baseline(
      comparator_source = source,
      comparator_model = model
    ) %>%
      mutate(rank = rank)
  }) %>%
  mutate(
    comparison_family = "Post-selection exploratory top-5 comparison",
    delong_p_adjust_bh = p.adjust(delong_p_value, method = "BH")
  ) %>%
  arrange(rank)

save_table(
  top5_delong,
  "09D_top5_delong_vs_baseline.csv"
)

# Additional prespecified DeLong comparisons for the locked candidate and the
# maximal all-modality model. These are included because they are central thesis
# models even if either falls outside the empirical top 5.
central_integrated_delong <- bind_rows(
  paired_delong_vs_baseline("no_aspirin", candidate_model),
  paired_delong_vs_baseline("no_aspirin", all_modalities_model)
) %>%
  mutate(
    comparison_family = "Prespecified central integrated-model comparison",
    delong_p_adjust_bh = p.adjust(delong_p_value, method = "BH")
  )

# Final prespecified AUC comparison table
prespecified_model_registry <- tribble(
  ~source, ~model, ~model_group, ~display_order,
  "original", baseline_model, "Clinical baseline", 1,
  "original", expanded_clinical_model, "Expanded clinical", 2,
  "original", "lipoproteins", "Single modality", 3,
  "original", "lipids", "Single modality", 4,
  "original", "small_molecules", "Single modality", 5,
  "original", "cytokines", "Single modality", 6,
  "original", "bile_acids", "Single modality", 7,
  "original", "pathology", "Single modality", 8,
  "no_aspirin", candidate_model, "Locked candidate", 9,
  "no_aspirin", all_modalities_model, "All modalities", 10
)

prespecified_delong <- bind_rows(
  monomodality_delong,
  expanded_clinical_delong,
  central_integrated_delong
) %>%
  select(
    source,
    model,
    delta_auc_vs_baseline,
    delta_auc_ci_low,
    delta_auc_ci_high,
    delong_p_value,
    delong_p_adjust_bh
  )

prespecified_auc_comparison <- prespecified_model_registry %>%
  left_join(
    all_auc_performance,
    by = c("source", "model")
  ) %>%
  left_join(
    prespecified_delong,
    by = c("source", "model")
  ) %>%
  mutate(
    model_label = map_chr(model, make_model_label),
    delta_auc_vs_baseline = if_else(
      model == baseline_model,
      0,
      delta_auc_vs_baseline
    ),
    comparison_status = if_else(
      model == baseline_model,
      "Reference",
      "Paired DeLong versus baseline clinical"
    )
  ) %>%
  arrange(display_order)

save_table(
  prespecified_auc_comparison,
  "09A_prespecified_auc_comparison_vs_baseline.csv"
)

# Final top-5 AUC table with the compact clinical benchmark included.
top5_auc_comparison <- top_models %>%
  select(
    source,
    model,
    rank,
    pooled_auc,
    auc_ci_low,
    auc_ci_high,
    n,
    n_events,
    n_nonevents
  ) %>%
  left_join(
    top5_delong %>%
      select(
        source,
        model,
        delta_auc_vs_baseline,
        delta_auc_ci_low,
        delta_auc_ci_high,
        delong_p_value,
        delong_p_adjust_bh
      ),
    by = c("source", "model")
  ) %>%
  mutate(
    model_label = map_chr(model, make_model_label),
    model_group = "Top 5 model"
  ) %>%
  select(
    rank,
    source,
    model,
    model_label,
    model_group,
    everything()
  )

baseline_top5_row <- baseline_performance %>%
  transmute(
    rank = NA_integer_,
    source,
    model,
    model_label = "Baseline clinical",
    model_group = "Clinical baseline",
    pooled_auc,
    auc_ci_low,
    auc_ci_high,
    n,
    n_events,
    n_nonevents,
    delta_auc_vs_baseline = 0,
    delta_auc_ci_low = NA_real_,
    delta_auc_ci_high = NA_real_,
    delong_p_value = NA_real_,
    delong_p_adjust_bh = NA_real_
  )

expanded_clinical_top5_row <- all_auc_performance %>%
  filter(source == "original", model == expanded_clinical_model) %>%
  left_join(
    expanded_clinical_delong %>%
      select(source, model, delta_auc_vs_baseline, delta_auc_ci_low, delta_auc_ci_high, delong_p_value, delong_p_adjust_bh),
    by = c("source", "model")
  ) %>%
  transmute(
    rank = NA_integer_,
    source,
    model,
    model_label = "Expanded clinical",
    model_group = "Expanded clinical",
    pooled_auc,
    auc_ci_low,
    auc_ci_high,
    n,
    n_events,
    n_nonevents,
    delta_auc_vs_baseline,
    delta_auc_ci_low,
    delta_auc_ci_high,
    delong_p_value,
    delong_p_adjust_bh
  )

top5_auc_comparison_with_baseline <- bind_rows(
  top5_auc_comparison,
  expanded_clinical_top5_row,
  baseline_top5_row
)

save_table(
  top5_auc_comparison_with_baseline,
  "09B_top5_auc_comparison_vs_baseline.csv"
)

# Figure 09A: prespecified thesis model panel versus baseline clinical
#
# Horizontal confidence intervals show each model's pooled out-of-fold AUC 95%
# CI. The vertical dashed line is the baseline clinical AUC. The right-hand text
# reports delta AUC and the raw paired DeLong p-value for direct thesis reading.
prespecified_plot_data <- prespecified_auc_comparison %>%
  mutate(
    model_label_wrapped = str_wrap(model_label, width = 36),
    model_label_wrapped = factor(
      model_label_wrapped,
      levels = rev(model_label_wrapped)
    ),
    comparison_text = case_when(
      model == baseline_model ~ paste0(
        "AUC ", formatC(pooled_auc, format = "f", digits = 3),
        " | reference"
      ),
      TRUE ~ paste0(
        "Δ ", format_delta_auc(delta_auc_vs_baseline),
        " | ", format_p_value(delong_p_value)
      )
    )
  )

prespecified_label_x <- max(
  prespecified_plot_data$auc_ci_high,
  na.rm = TRUE
) + 0.012

prespecified_x_min <- max(
  0.40,
  min(c(prespecified_plot_data$auc_ci_low, prespecified_plot_data$pooled_auc), na.rm = TRUE) - 0.03
)

prespecified_x_max <- min(
  0.99,
  max(c(prespecified_plot_data$auc_ci_high, prespecified_plot_data$pooled_auc), na.rm = TRUE) + 0.22
)

prespecified_auc_plot <- ggplot(
  prespecified_plot_data,
  aes(
    x = pooled_auc,
    y = model_label_wrapped
  )
) +
  geom_vline(
    xintercept = baseline_auc,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey45"
  ) +
  geom_segment(
    aes(
      x = auc_ci_low,
      xend = auc_ci_high,
      yend = model_label_wrapped
    ),
    linewidth = 1.0,
    colour = "#5B5B5B"
  ) +
  geom_point(
    aes(fill = model_group),
    shape = 21,
    size = 4.6,
    stroke = 0.7,
    colour = "white"
  ) +
  geom_text(
    aes(
      x = prespecified_label_x,
      label = comparison_text
    ),
    hjust = 0,
    size = 3.45,
    colour = "#2F3A4A"
  ) +
  scale_fill_manual(values = model_group_colours) +
  scale_x_continuous(
    breaks = seq(0.40, 0.95, by = 0.05)
  ) +
  labs(
    title = "Final model discrimination against the clinical benchmark",
    subtitle = paste0(
      "Pooled outer-fold AUC with 95% DeLong confidence intervals. ",
      "All p-values are from paired DeLong tests versus the baseline clinical model."
    ),
    x = "Pooled out-of-fold AUC",
    fill = NULL
  ) +
  coord_cartesian(xlim = c(prespecified_x_min, prespecified_x_max), clip = "off") +
  thesis_theme

save_plot(
  prespecified_auc_plot,
  "09A_prespecified_auc_comparison_vs_baseline.png",
  width = 12.8,
  height = 7.8
)

# Figure 09B: top 5 final models versus baseline clinical
#
# The top 5 are ordered by AUC rank. The clinical baseline is added as an
# explicit reference rather than being included in the top-5 ranking itself.
top5_plot_data <- top5_auc_comparison_with_baseline %>%
  mutate(
    plot_order = case_when(
      !is.na(rank) ~ rank,
      model == expanded_clinical_model ~ top_n_models + 1,
      model == baseline_model ~ top_n_models + 2,
      TRUE ~ top_n_models + 3
    )
  ) %>%
  arrange(plot_order) %>%
  mutate(
    display_label = case_when(
      !is.na(rank) ~ paste0(rank, ". ", model_label),
      model == expanded_clinical_model ~ "Expanded clinical",
      model == baseline_model ~ "Baseline clinical",
      TRUE ~ model_label
    ),
    display_label = str_wrap(display_label, width = 48),
    display_label = factor(
      display_label,
      levels = rev(display_label)
    ),
    comparison_text = case_when(
      model == baseline_model ~ paste0(
        "AUC ", formatC(pooled_auc, format = "f", digits = 3),
        " | reference"
      ),
      TRUE ~ paste0(
        "Δ ", format_delta_auc(delta_auc_vs_baseline),
        " | ", format_p_value(delong_p_value)
      )
    )
  )

top5_label_x <- max(
  top5_plot_data$auc_ci_high,
  na.rm = TRUE
) + 0.012

top5_x_min <- max(
  0.55,
  min(c(top5_plot_data$auc_ci_low, top5_plot_data$pooled_auc), na.rm = TRUE) - 0.03
)

top5_x_max <- min(
  0.99,
  max(c(top5_plot_data$auc_ci_high, top5_plot_data$pooled_auc), na.rm = TRUE) + 0.22
)

top5_auc_plot <- ggplot(
  top5_plot_data,
  aes(
    x = pooled_auc,
    y = display_label
  )
) +
  geom_vline(
    xintercept = baseline_auc,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey45"
  ) +
  geom_segment(
    aes(
      x = auc_ci_low,
      xend = auc_ci_high,
      yend = display_label
    ),
    linewidth = 1.0,
    colour = "#5B5B5B"
  ) +
  geom_point(
    aes(fill = model_group),
    shape = 21,
    size = 4.6,
    stroke = 0.7,
    colour = "white"
  ) +
  geom_text(
    aes(
      x = top5_label_x,
      label = comparison_text
    ),
    hjust = 0,
    size = 3.30,
    colour = "#2F3A4A"
  ) +
  scale_fill_manual(values = model_group_colours) +
  scale_x_continuous(
    breaks = seq(0.55, 0.95, by = 0.05)
  ) +
  labs(
    title = "Highest-performing integrated models compared with the clinical benchmarks",
    subtitle = paste0(
      "The highest-performing aspirin-excluded models are shown alongside expanded and baseline clinical models.\n",
      "Error bars show 95% DeLong AUC confidence intervals."
    ),
    x = "Pooled out-of-fold AUC",
    fill = NULL
  ) +
  coord_cartesian(xlim = c(top5_x_min, top5_x_max), clip = "off") +
  thesis_theme

save_plot(
  top5_auc_plot,
  "09B_top5_auc_comparison_vs_baseline.png",
  width = 14.5,
  height = 9.3
)

# Optional final biomarker synthesis from Script 08
#
# Script 08 already performs the detailed cross-model feature recurrence work.
# Rather than duplicate that analysis here, this section extracts the broader
# recurrent-biomarker summary when the saved Script 08 RDS is available.
biomarker_output_file <- biomarker_output_candidates[
  file.exists(biomarker_output_candidates)
][1]

cross_model_recurrent_biomarkers <- tibble()

if (!is.na(biomarker_output_file)) {
  biomarker_outputs <- readRDS(biomarker_output_file)

  if ("recurrent_summary" %in% names(biomarker_outputs)) {
    cross_model_recurrent_biomarkers <- biomarker_outputs$recurrent_summary %>%
      filter(models_stable_n >= 2) %>%
      arrange(
        desc(models_stable_n),
        desc(models_selected_n),
        desc(mean_model_selection_frequency),
        desc(mean_abs_coefficient)
      ) %>%
      select(
        recurrence_rank,
        feature,
        feature_label,
        feature_modality,
        models_selected_n,
        models_stable_n,
        pooled_selection_frequency,
        mean_model_selection_frequency,
        mean_coefficient,
        mean_abs_coefficient,
        direction_agreement,
        recurrent_direction,
        model_labels
      )

    save_table(
      cross_model_recurrent_biomarkers,
      "09G_cross_model_recurrent_biomarker_summary.csv"
    )
  }


} else {
  message(
    "Script 08 biomarker output was not found. Model/AUC synthesis was completed, ",
    "but optional 09G cross-model biomarker table was skipped."
  )
}

# Save compact RDS containing all final synthesis objects.
final_synthesis_outputs <- list(
  settings = list(
    original_output_file = original_output_file,
    no_aspirin_output_file = no_aspirin_output_file,
    biomarker_output_file = biomarker_output_file,
    baseline_model = baseline_model,
    single_modality_models = single_modality_models,
    candidate_model = candidate_model,
    all_modalities_model = all_modalities_model,
    top_n_models = top_n_models,
    expanded_clinical_model = expanded_clinical_model,
    stable_selection_cutoff = stable_selection_cutoff,
    participant_id_column = participant_id_column
  ),
  pairing_audit = pairing_audit,
  all_auc_performance = all_auc_performance,
  prespecified_auc_comparison = prespecified_auc_comparison,
  top5_auc_comparison = top5_auc_comparison_with_baseline,
  monomodality_delong = monomodality_delong,
  top5_delong = top5_delong,
  central_integrated_delong = central_integrated_delong,
  top5_recurrent_stable_biomarkers = top5_recurrent_stable_biomarkers,
  candidate_stable_biomarkers = candidate_stable_biomarkers_direct,
  cross_model_recurrent_biomarkers = cross_model_recurrent_biomarkers
)

saveRDS(
  final_synthesis_outputs,
  file.path(
    output_model_dir,
    "09_final_model_and_biomarker_synthesis_outputs.rds"
  )
)

# Console summary
message("Script 09 final thesis synthesis complete.")
message("Baseline clinical pooled AUC: ", round(baseline_auc, 3))
message("Top final model: ", top_models$model[1])
message("Top final model pooled AUC: ", round(top_models$pooled_auc[1], 3))
message("Locked candidate pooled AUC: ", round(
  prespecified_auc_comparison %>%
    filter(model == candidate_model) %>%
    pull(pooled_auc),
  3
))
message("All-modality pooled AUC: ", round(
  prespecified_auc_comparison %>%
    filter(model == all_modalities_model) %>%
    pull(pooled_auc),
  3
))
message("Main figures saved to: ", output_figure_dir)
message("Main tables saved to: ", output_table_dir)

print(
  prespecified_auc_comparison %>%
    select(
      display_order,
      model_label,
      pooled_auc,
      auc_ci_low,
      auc_ci_high,
      delta_auc_vs_baseline,
      delong_p_value,
      delong_p_adjust_bh
    )
)
