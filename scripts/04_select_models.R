# Repository stage: selection
# Source supplied as: 04_model_selection_and_delong_comparison_no_aspirin(6).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# 04 Model selection and DeLong comparison: expanded clinical models excluding aspirin
#
# Project:
#   Cardiovascular Risk Stratification in Liver Cirrhosis
#
# Purpose:
#   Simplified model-selection script.
#
#   Script 03 generated the full 191-model nested cross-validation comparison.
#   This script does not train new models. It uses the saved out-of-fold
#   predictions and pooled performance tables from Script 03 to:
#
#   1. Rank all models by pooled out-of-fold AUC.
#   2. Compare each model against the expanded clinical benchmark.
#   3. Perform paired DeLong tests using matched out-of-fold predictions.
#   4. Save a simple top-10 model panel for Script 05.
#   5. Produce a small set of model-selection figures.
#
#   Detailed clinical characterisation, confusion matrices, sensitivity/
#   specificity review, top analytes, SHAP, and biological interpretation are
#   intentionally left for Script 05.

library(tidyverse)
library(pROC)
library(PRROC)

set.seed(123)

# Output folders
dir.create("outputs/tables/model_selection_no_aspirin", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures/model_selection_no_aspirin", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/models/model_selection_no_aspirin", recursive = TRUE, showWarnings = FALSE)

# Settings
top_n_panel <- 10
top_n_plot <- 20
top_n_roc <- 5

baseline_benchmark_model <- "clinical_baseline"
expanded_benchmark_model <- "clinical_expanded_no_aspirin"

script03_output_candidates <- c(ccrvi_config$files$no_aspirin_results)

script03_output_file <- script03_output_candidates[
  file.exists(script03_output_candidates)
][1]

if (is.na(script03_output_file)) {
  stop(
    "No no-aspirin Script 03 RDS output was found. Expected one of:\n",
    paste(script03_output_candidates, collapse = "\n")
  )
}

# Helper functions
save_table <- function(x, filename) {
  write.csv(
    x,
    file.path("outputs/tables/model_selection_no_aspirin", filename),
    row.names = FALSE
  )
}

save_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path("outputs/figures/model_selection_no_aspirin", filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

normalise_omics_label <- function(x) {
  if_else(is.na(x) | x == "", "clinical_only", x)
}

make_model_label <- function(model_name) {
  model_name %>%
    str_replace(
      "^clinical_expanded_no_aspirin",
      "Expanded clinical (aspirin excluded)"
    ) %>%
    str_replace("^clinical_expanded", "Expanded clinical") %>%
    str_replace("^clinical_baseline", "Baseline clinical") %>%
    str_replace_all("lipoproteins", "Lipoproteins") %>%
    str_replace_all("lipids", "Lipids") %>%
    str_replace_all("small_molecules", "Small molecules") %>%
    str_replace_all("cytokines", "Cytokines") %>%
    str_replace_all("bile_acids", "Bile acids") %>%
    str_replace_all("pathology", "Pathology") %>%
    str_replace_all("_", " + ") %>%
    str_squish()
}

calc_auprc <- function(y_true, pred) {
  if (length(unique(y_true)) < 2) {
    return(NA_real_)
  }

  tryCatch(
    PRROC::pr.curve(
      scores.class0 = pred[y_true == 1],
      scores.class1 = pred[y_true == 0],
      curve = FALSE
    )$auc.integral,
    error = function(e) NA_real_
  )
}

calc_threshold_metrics <- function(y_true, pred) {
  roc_obj <- pROC::roc(y_true, pred, direction = "<", quiet = TRUE)

  best <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
    transpose = FALSE
  )

  predicted_class <- ifelse(pred >= best$threshold, 1, 0)

  tibble(
    threshold = as.numeric(best$threshold),
    sensitivity = as.numeric(best$sensitivity),
    specificity = as.numeric(best$specificity),
    ppv = as.numeric(best$ppv),
    npv = as.numeric(best$npv),
    accuracy = mean(predicted_class == y_true, na.rm = TRUE),
    balanced_accuracy = mean(c(best$sensitivity, best$specificity), na.rm = TRUE)
  )
}

safe_delong <- function(prediction_df, model_name, benchmark_name) {
  pair_df <- prediction_df %>%
    filter(.data$model %in% c(model_name, benchmark_name)) %>%
    select(model, outcome, pred) %>%
    group_by(model) %>%
    mutate(row_id = row_number()) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = row_id,
      names_from = model,
      values_from = pred
    )

  # Add the outcome vector once, using the prediction order from the model being
  # tested. Script 03 creates matched out-of-fold predictions in the same row
  # order for every model, so row_id is a valid pairing key.
  outcome_df <- prediction_df %>%
    filter(.data$model == model_name) %>%
    select(outcome) %>%
    mutate(row_id = row_number())

  pair_df <- pair_df %>%
    left_join(outcome_df, by = "row_id")

  if (!all(c(model_name, benchmark_name, "outcome") %in% names(pair_df))) {
    return(
      tibble(
        model = model_name,
        benchmark_model = benchmark_name,
        auc_model = NA_real_,
        auc_benchmark = NA_real_,
        delta_auc = NA_real_,
        delong_p_value = NA_real_
      )
    )
  }

  pair_df <- pair_df %>%
    filter(
      !is.na(.data$outcome),
      !is.na(.data[[model_name]]),
      !is.na(.data[[benchmark_name]])
    )

  if (nrow(pair_df) < 10 || n_distinct(pair_df$outcome) != 2) {
    return(
      tibble(
        model = model_name,
        benchmark_model = benchmark_name,
        auc_model = NA_real_,
        auc_benchmark = NA_real_,
        delta_auc = NA_real_,
        delong_p_value = NA_real_
      )
    )
  }

  roc_model <- pROC::roc(
    pair_df$outcome,
    pair_df[[model_name]],
    direction = "<",
    quiet = TRUE
  )

  roc_benchmark <- pROC::roc(
    pair_df$outcome,
    pair_df[[benchmark_name]],
    direction = "<",
    quiet = TRUE
  )

  test_result <- tryCatch(
    pROC::roc.test(
      roc_model,
      roc_benchmark,
      method = "delong",
      paired = TRUE
    ),
    error = function(e) NULL
  )

  tibble(
    model = model_name,
    benchmark_model = benchmark_name,
    auc_model = as.numeric(pROC::auc(roc_model)),
    auc_benchmark = as.numeric(pROC::auc(roc_benchmark)),
    delta_auc = auc_model - auc_benchmark,
    delong_p_value = ifelse(is.null(test_result), NA_real_, test_result$p.value)
  )
}


# Load Script 03 outputs
if (!file.exists(script03_output_file)) {
  stop("Script 03 output file not found: ", script03_output_file)
}

script03_outputs <- readRDS(script03_output_file)

required_objects <- c(
  "modality_predictions",
  "pooled_modality_performance",
  "modality_fold_performance"
)

missing_objects <- setdiff(required_objects, names(script03_outputs))

if (length(missing_objects) > 0) {
  stop(
    "Script 03 output object is missing required object(s): ",
    paste(missing_objects, collapse = ", ")
  )
}

modality_predictions <- script03_outputs$modality_predictions %>%
  mutate(
    model = as.character(model),
    outcome = as.numeric(outcome),
    pred = as.numeric(pred)
  )

pooled_modality_performance <- script03_outputs$pooled_modality_performance %>%
  mutate(
    model = as.character(model),
    omics_modalities = normalise_omics_label(omics_modalities)
  )

modality_fold_performance <- script03_outputs$modality_fold_performance %>%
  mutate(
    model = as.character(model),
    omics_modalities = normalise_omics_label(omics_modalities)
  )

# Basic prediction coverage audit
prediction_coverage <- modality_predictions %>%
  group_by(model) %>%
  summarise(
    n_predictions = n(),
    n_events = sum(outcome == 1, na.rm = TRUE),
    n_nonevents = sum(outcome == 0, na.rm = TRUE),
    event_rate = mean(outcome == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(model)

save_table(prediction_coverage, "04_no_aspirin_prediction_coverage_audit.csv")

# Build all-model ranking table
fold_summary <- modality_fold_performance %>%
  group_by(model) %>%
  summarise(
    mean_validation_auc = mean(validation_auc, na.rm = TRUE),
    sd_validation_auc = sd(validation_auc, na.rm = TRUE),
    min_validation_auc = min(validation_auc, na.rm = TRUE),
    max_validation_auc = max(validation_auc, na.rm = TRUE),
    mean_validation_auprc = mean(validation_auprc, na.rm = TRUE),
    sd_validation_auprc = sd(validation_auprc, na.rm = TRUE),
    mean_selected_features = mean(n_model_features_selected, na.rm = TRUE),
    .groups = "drop"
  )

expanded_reference <- pooled_modality_performance %>%
  filter(model == expanded_benchmark_model) %>%
  slice(1)

baseline_reference <- pooled_modality_performance %>%
  filter(model == baseline_benchmark_model) %>%
  slice(1)

if (nrow(expanded_reference) == 0) {
  stop("Expanded clinical benchmark not found: ", expanded_benchmark_model)
}

all_model_ranking <- pooled_modality_performance %>%
  left_join(fold_summary, by = "model") %>%
  mutate(
    model_label = make_model_label(model),
    model_rank_auc = rank(-pooled_auc, ties.method = "min"),
    complexity_class = case_when(
      n_omics_modalities == 0 ~ "clinical_only",
      n_omics_modalities == 1 ~ "single_modality",
      n_omics_modalities <= 3 ~ "limited_multimodal",
      TRUE ~ "high_dimensional_multimodal"
    ),
    delta_auc_vs_expanded = pooled_auc - expanded_reference$pooled_auc,
    delta_auprc_vs_expanded = pooled_auprc - expanded_reference$pooled_auprc,
    delta_brier_vs_expanded = brier_score - expanded_reference$brier_score,
    delta_auc_vs_baseline = if (nrow(baseline_reference) == 1) {
      pooled_auc - baseline_reference$pooled_auc
    } else {
      NA_real_
    }
  ) %>%
  arrange(desc(pooled_auc))

save_table(all_model_ranking, "04_no_aspirinA_all_models_ranked_by_auc.csv")

# DeLong tests against expanded clinical and baseline clinical
delong_vs_expanded <- all_model_ranking %>%
  filter(model != expanded_benchmark_model) %>%
  pull(model) %>%
  map_dfr(~ safe_delong(
    prediction_df = modality_predictions,
    model_name = .x,
    benchmark_name = expanded_benchmark_model
  )) %>%
  mutate(
    comparison = "vs expanded clinical (aspirin excluded)",
    delong_p_adjust_bh = p.adjust(delong_p_value, method = "BH")
  )

delong_vs_baseline <- if (nrow(baseline_reference) == 1) {
  all_model_ranking %>%
    filter(model != baseline_benchmark_model) %>%
    pull(model) %>%
    map_dfr(~ safe_delong(
      prediction_df = modality_predictions,
      model_name = .x,
      benchmark_name = baseline_benchmark_model
    )) %>%
    mutate(
      comparison = "vs baseline clinical",
      delong_p_adjust_bh = p.adjust(delong_p_value, method = "BH")
    )
} else {
  tibble()
}

delong_results <- bind_rows(
  delong_vs_expanded,
  delong_vs_baseline
) %>%
  arrange(comparison, desc(delta_auc))

save_table(delong_results, "04_no_aspirinB_delong_tests_against_clinical_benchmarks.csv")

# Add DeLong results and threshold metrics to the main model-selection table.
threshold_metric_table <- modality_predictions %>%
  group_split(model) %>%
  map_dfr(function(model_df) {
    metrics <- calc_threshold_metrics(
      y_true = model_df$outcome,
      pred = model_df$pred
    )

    metrics %>%
      mutate(model = model_df$model[1], .before = 1)
  })

required_threshold_columns <- c(
  "model",
  "threshold",
  "sensitivity",
  "specificity",
  "ppv",
  "npv",
  "accuracy",
  "balanced_accuracy"
)

missing_threshold_columns <- setdiff(
  required_threshold_columns,
  names(threshold_metric_table)
)

if (length(missing_threshold_columns) > 0) {
  stop(
    "Threshold metric table is missing: ",
    paste(missing_threshold_columns, collapse = ", ")
  )
}

# The pooled performance table from Script 03 may already contain threshold
# metrics. Joining the freshly calculated threshold table can therefore create
# .x/.y column suffixes. To avoid downstream ambiguity, remove any existing
# threshold columns before joining the fresh threshold metrics.
threshold_columns <- c(
  "threshold",
  "sensitivity",
  "specificity",
  "ppv",
  "npv",
  "accuracy",
  "balanced_accuracy"
)

model_selection_table <- all_model_ranking %>%
  select(-any_of(threshold_columns)) %>%
  left_join(
    delong_vs_expanded %>%
      select(
        model,
        delong_delta_auc_vs_expanded = delta_auc,
        delong_p_vs_expanded = delong_p_value,
        delong_p_adjust_bh_vs_expanded = delong_p_adjust_bh
      ),
    by = "model"
  ) %>%
  left_join(threshold_metric_table, by = "model") %>%
  mutate(
    delong_delta_auc_vs_expanded = if_else(
      model == expanded_benchmark_model,
      0,
      delong_delta_auc_vs_expanded
    ),
    delong_p_vs_expanded = if_else(
      model == expanded_benchmark_model,
      NA_real_,
      delong_p_vs_expanded
    )
  ) %>%
  arrange(desc(pooled_auc))

missing_model_selection_columns <- setdiff(
  c(
    "model",
    "pooled_auc",
    "sensitivity",
    "specificity",
    "ppv",
    "npv"
  ),
  names(model_selection_table)
)

if (length(missing_model_selection_columns) > 0) {
  stop(
    "Model selection table is missing: ",
    paste(missing_model_selection_columns, collapse = ", ")
  )
}

save_table(model_selection_table, "04_no_aspirinC_model_selection_table_with_delong.csv")

# Top-10 model panel for Script 05
top10_model_panel <- model_selection_table %>%
  slice_head(n = top_n_panel) %>%
  select(
    rank = model_rank_auc,
    model,
    model_label,
    model_family,
    omics_modalities,
    n_omics_modalities,
    complexity_class,
    pooled_auc,
    auc_ci_low,
    auc_ci_high,
    pooled_auprc,
    delta_auc_vs_expanded,
    delong_p_vs_expanded,
    brier_score,
    calibration_intercept,
    calibration_slope,
    sensitivity,
    specificity,
    ppv,
    npv,
    mean_selected_features
  )

# Include benchmarks underneath the top 10 if they are not already in the top 10.
benchmark_panel <- model_selection_table %>%
  filter(model %in% c(baseline_benchmark_model, expanded_benchmark_model)) %>%
  filter(!model %in% top10_model_panel$model) %>%
  select(
    rank = model_rank_auc,
    model,
    model_label,
    model_family,
    omics_modalities,
    n_omics_modalities,
    complexity_class,
    pooled_auc,
    auc_ci_low,
    auc_ci_high,
    pooled_auprc,
    delta_auc_vs_expanded,
    delong_p_vs_expanded,
    brier_score,
    calibration_intercept,
    calibration_slope,
    sensitivity,
    specificity,
    ppv,
    npv,
    mean_selected_features
  )

script05_selected_models <- bind_rows(top10_model_panel, benchmark_panel) %>%
  arrange(rank)

save_table(top10_model_panel, "04_no_aspirinD_top10_model_panel.csv")
save_table(script05_selected_models, "04_no_aspirinD_script05_selected_models.csv")

# ROC and PR curve data for top models and clinical benchmarks
curve_model_names <- unique(c(
  baseline_benchmark_model,
  expanded_benchmark_model,
  top10_model_panel %>%
    slice_head(n = top_n_roc) %>%
    pull(model)
))

curve_predictions <- modality_predictions %>%
  filter(model %in% curve_model_names)

roc_curve_data <- curve_predictions %>%
  group_by(model) %>%
  group_modify(~ {
    roc_obj <- pROC::roc(.x$outcome, .x$pred, direction = "<", quiet = TRUE)

    tibble(
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      false_positive_rate = 1 - specificity
    )
  }) %>%
  ungroup() %>%
  left_join(
    model_selection_table %>%
      select(model, model_label, pooled_auc, auc_ci_low, auc_ci_high),
    by = "model"
  )

pr_curve_data <- curve_predictions %>%
  group_by(model) %>%
  group_modify(~ {
    thresholds <- sort(unique(.x$pred), decreasing = TRUE)

    map_dfr(thresholds, function(th) {
      predicted_class <- ifelse(.x$pred >= th, 1, 0)

      tp <- sum(predicted_class == 1 & .x$outcome == 1)
      fp <- sum(predicted_class == 1 & .x$outcome == 0)
      fn <- sum(predicted_class == 0 & .x$outcome == 1)

      tibble(
        threshold = th,
        precision = ifelse(tp + fp == 0, NA_real_, tp / (tp + fp)),
        recall = ifelse(tp + fn == 0, NA_real_, tp / (tp + fn))
      )
    })
  }) %>%
  ungroup() %>%
  filter(!is.na(precision), !is.na(recall)) %>%
  left_join(
    model_selection_table %>%
      select(model, model_label, pooled_auprc),
    by = "model"
  )

save_table(roc_curve_data, "04_no_aspirinE_roc_curve_data_top_models.csv")
save_table(pr_curve_data, "04_no_aspirinE_precision_recall_curve_data_top_models.csv")

# Figures

# 1. Top 20 models ranked by pooled AUC.
top20_auc_plot <- model_selection_table %>%
  slice_head(n = top_n_plot) %>%
  mutate(
    model_label_wrapped = str_wrap(model_label, width = 42),
    model_label_wrapped = factor(model_label_wrapped, levels = rev(model_label_wrapped))
  ) %>%
  ggplot(
    aes(
      x = model_label_wrapped,
      y = pooled_auc,
      fill = complexity_class
    )
  ) +
  geom_col() +
  geom_errorbar(
    aes(ymin = auc_ci_low, ymax = auc_ci_high),
    width = 0.25
  ) +
  coord_flip() +
  labs(
    title = "Top models by pooled out-of-fold AUC",
    subtitle = "Error bars show 95% AUC confidence intervals",
    x = NULL,
    y = "Pooled AUC",
    fill = "Model complexity"
  ) +
  theme_minimal(base_size = 12)

save_plot(top20_auc_plot, "04_no_aspirinA_top20_models_by_auc.png", width = 10, height = 8)

# 2. Delta AUC versus expanded clinical (aspirin excluded).
delta_auc_plot <- model_selection_table %>%
  filter(model != expanded_benchmark_model) %>%
  slice_head(n = top_n_plot) %>%
  mutate(
    model_label_wrapped = str_wrap(model_label, width = 42),
    model_label_wrapped = factor(model_label_wrapped, levels = rev(model_label_wrapped)),
    significance_label = case_when(
      delong_p_vs_expanded < 0.05 ~ "DeLong p < 0.05",
      is.na(delong_p_vs_expanded) ~ "Not tested",
      TRUE ~ "DeLong p >= 0.05"
    )
  ) %>%
  ggplot(
    aes(
      x = model_label_wrapped,
      y = delta_auc_vs_expanded,
      fill = significance_label
    )
  ) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Model improvement versus expanded clinical (aspirin excluded)",
    subtitle = "Bars show delta pooled AUC; colour indicates paired DeLong result",
    x = NULL,
    y = "Delta pooled AUC versus expanded clinical (aspirin excluded)",
    fill = "Comparison"
  ) +
  theme_minimal(base_size = 12)

save_plot(delta_auc_plot, "04_no_aspirinB_delta_auc_vs_expanded_clinical.png", width = 10, height = 8)

# 3. Complexity versus AUC.
complexity_auc_plot <- model_selection_table %>%
  ggplot(
    aes(
      x = mean_selected_features,
      y = pooled_auc,
      colour = complexity_class
    )
  ) +
  geom_point(alpha = 0.75, size = 2.5) +
  geom_point(
    data = top10_model_panel,
    aes(x = mean_selected_features, y = pooled_auc),
    inherit.aes = FALSE,
    shape = 21,
    size = 4,
    stroke = 1.1,
    fill = "white"
  ) +
  labs(
    title = "Model complexity versus pooled AUC",
    subtitle = "Open circles indicate top-10 AUC models",
    x = "Mean selected features across validation folds",
    y = "Pooled AUC",
    colour = "Model complexity"
  ) +
  theme_minimal(base_size = 12)

save_plot(complexity_auc_plot, "04_no_aspirinC_complexity_vs_auc.png", width = 8, height = 6)

# 4. ROC curves for top 5 plus clinical benchmarks.
roc_label_table <- model_selection_table %>%
  filter(model %in% curve_model_names) %>%
  mutate(
    curve_label = paste0(
      str_wrap(model_label, width = 32),
      " | AUC ",
      round(pooled_auc, 3)
    )
  ) %>%
  select(model, curve_label)

roc_curve_plot <- roc_curve_data %>%
  left_join(roc_label_table, by = "model") %>%
  ggplot(
    aes(
      x = false_positive_rate,
      y = sensitivity,
      colour = curve_label
    )
  ) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "ROC curves for top candidate models",
    subtitle = "Top 5 by pooled AUC plus clinical benchmarks",
    x = "False positive rate",
    y = "Sensitivity",
    colour = "Model"
  ) +
  theme_minimal(base_size = 11)

save_plot(roc_curve_plot, "04_no_aspirinD_roc_curves_top5_plus_benchmarks.png", width = 9, height = 7)

# 5. Precision-recall curves for top 5 plus clinical benchmarks.
pr_label_table <- model_selection_table %>%
  filter(model %in% curve_model_names) %>%
  mutate(
    curve_label = paste0(
      str_wrap(model_label, width = 32),
      " | AUPRC ",
      round(pooled_auprc, 3)
    )
  ) %>%
  select(model, curve_label)

pr_curve_plot <- pr_curve_data %>%
  left_join(pr_label_table, by = "model") %>%
  ggplot(
    aes(
      x = recall,
      y = precision,
      colour = curve_label
    )
  ) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  labs(
    title = "Precision-recall curves for top candidate models",
    subtitle = "Top 5 by pooled AUC plus clinical benchmarks",
    x = "Recall",
    y = "Precision",
    colour = "Model"
  ) +
  theme_minimal(base_size = 11)

save_plot(pr_curve_plot, "04_no_aspirinE_precision_recall_curves_top5_plus_benchmarks.png", width = 9, height = 7)

# Save Script 04 objects
script04_outputs <- list(
  settings = list(
    top_n_panel = top_n_panel,
    top_n_plot = top_n_plot,
    top_n_roc = top_n_roc,
    baseline_benchmark_model = baseline_benchmark_model,
    expanded_benchmark_model = expanded_benchmark_model
  ),
  prediction_coverage = prediction_coverage,
  all_model_ranking = all_model_ranking,
  delong_results = delong_results,
  model_selection_table = model_selection_table,
  top10_model_panel = top10_model_panel,
  script05_selected_models = script05_selected_models,
  roc_curve_data = roc_curve_data,
  pr_curve_data = pr_curve_data
)

saveRDS(
  script04_outputs,
  "outputs/models/model_selection_no_aspirin/04_no_aspirin_model_selection_outputs.rds"
)

message("Script 04 no-aspirin model selection complete.")
message("Models ranked: ", nrow(model_selection_table))
message("Top model: ", top10_model_panel$model[1])
message("Top-10 panel saved to outputs/tables/model_selection_no_aspirin/04_no_aspirinD_top10_model_panel.csv")
message("Script 05 selected models saved to outputs/tables/model_selection_no_aspirin/04_no_aspirinD_script05_selected_models.csv")
print(top10_model_panel)
