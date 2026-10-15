# Repository stage: masld
# Source supplied as: 07_masld_stratification_candidate_model(3).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# 07 MASLD stratification of the selected MACE candidate model
#
# Project:
#   Cardiovascular Risk Stratification in Liver Cirrhosis
#
# Locked candidate:
#   Expanded clinical (aspirin excluded) +
#   Lipoproteins + Lipids + Small molecules
#
# Purpose:
#   Evaluate whether the selected candidate model behaves differently in
#   participants with and without MASLD, while keeping two distinct questions
#   methodologically separate:
#
#   Part A: Locked-model subgroup evaluation
#     Uses the original pooled outer-fold out-of-sample predictions from Script 03.
#     The model is not retrained, the overall threshold is not reselected, and
#     MASLD status is attached only after reconstructing the original validation
#     fold membership.
#
#   Part B: Exploratory MASLD-only model development
#     Repeats the candidate modelling strategy within the MASLD-only cohort using
#     nested cross-validation. This is intended for exploratory biological and
#     feature-stability analysis, not as validation of the locked candidate.
#
# Main outputs:
#   outputs/tables/masld_stratification/
#   outputs/figures/masld_stratification/
#   outputs/models/masld_stratification/
#
# Required inputs:
#   combined_multimodal_complete_cohort.csv
#   outputs/models/lead_candidate_no_aspirin/
#     05_no_aspirin_candidate_outputs.rds
#
# Optional input:
#   outputs/models/risk_score/
#     06_mace_risk_score_outputs.rds
#
# Important interpretation:
#   Subgroup estimates may be imprecise because the MASLD and non-MASLD groups
#   contain fewer patients and fewer MACE events than the complete cohort.
#   Confidence intervals and event counts should therefore be reported alongside
#   all subgroup performance estimates.

library(tidyverse)
library(caret)
library(glmnet)
library(pROC)
library(PRROC)
library(broom)

set.seed(123)

# Settings

compiled_cohort_file <- ccrvi_config$files$compiled_cohort

candidate_output_file <- paste0(
  "outputs/models/lead_candidate_no_aspirin/",
  "05_no_aspirin_candidate_outputs.rds"
)

risk_score_output_file <- paste0(
  "outputs/models/risk_score/",
  "06_mace_risk_score_outputs.rds"
)

output_table_dir <- "outputs/tables/masld_stratification"
output_figure_dir <- "outputs/figures/masld_stratification"
output_model_dir <- "outputs/models/masld_stratification"

candidate_model <- paste0(
  "clinical_expanded_no_aspirin_",
  "lipoproteins_",
  "lipids_",
  "small_molecules"
)

candidate_label <- paste(
  "Expanded clinical (aspirin excluded) + Lipoproteins +",
  "Lipids + Small molecules"
)

# This must match the seed and outer-fold count used in Script 03.
original_outer_seed <- 123
original_outer_folds <- 5

# Calibration groups are reduced automatically when a subgroup is small.
maximum_calibration_groups <- 5
minimum_patients_per_calibration_group <- 8

# Run the exploratory MASLD-only nested model.
# Set to FALSE when only the locked subgroup evaluation is required.
run_exploratory_masld_model <- ccrvi_config$run_exploratory_masld_model

# MASLD-only nested modelling settings.
masld_outer_folds <- 5
masld_inner_folds <- 5
alpha_grid <- ccrvi_config$alpha_grid
top_k_grid <- ccrvi_config$top_k_grid

dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_model_dir, recursive = TRUE, showWarnings = FALSE)

# Save helpers

save_table <- function(x, filename) {
  write.csv(
    x,
    file.path(output_table_dir, filename),
    row.names = FALSE
  )
}

save_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(output_figure_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

# General helpers

normalise_binary <- function(x) {
  x_character <- str_to_lower(str_trim(as.character(x)))

  case_when(
    x_character %in% c("1", "yes", "y", "true", "positive", "present") ~ 1,
    x_character %in% c("0", "no", "n", "false", "negative", "absent") ~ 0,
    TRUE ~ suppressWarnings(as.numeric(x_character))
  )
}

make_weights <- function(y) {
  n_positive <- sum(y == 1, na.rm = TRUE)
  n_negative <- sum(y == 0, na.rm = TRUE)

  if (n_positive == 0 || n_negative == 0) {
    return(rep(1, length(y)))
  }

  ifelse(y == 1, n_negative / n_positive, 1)
}

calc_auprc <- function(y_true, pred) {
  if (length(unique(y_true[!is.na(y_true)])) < 2) {
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

calc_calibration <- function(y_true, pred) {
  valid <- !is.na(y_true) & !is.na(pred)
  y_true <- y_true[valid]
  pred <- pred[valid]

  if (length(unique(y_true)) < 2) {
    return(
      tibble(
        calibration_intercept = NA_real_,
        calibration_slope = NA_real_,
        brier_score = mean((y_true - pred)^2)
      )
    )
  }

  pred_clipped <- pmin(pmax(pred, 1e-6), 1 - 1e-6)
  logit_pred <- qlogis(pred_clipped)

  calibration_fit <- tryCatch(
    glm(y_true ~ logit_pred, family = binomial),
    error = function(e) NULL
  )

  if (is.null(calibration_fit)) {
    return(
      tibble(
        calibration_intercept = NA_real_,
        calibration_slope = NA_real_,
        brier_score = mean((y_true - pred)^2)
      )
    )
  }

  tibble(
    calibration_intercept = unname(coef(calibration_fit)[1]),
    calibration_slope = unname(coef(calibration_fit)[2]),
    brier_score = mean((y_true - pred)^2)
  )
}

calc_metrics_at_threshold <- function(y_true, pred, threshold) {
  predicted_class <- ifelse(pred >= threshold, 1, 0)

  true_positive <- sum(predicted_class == 1 & y_true == 1, na.rm = TRUE)
  false_positive <- sum(predicted_class == 1 & y_true == 0, na.rm = TRUE)
  true_negative <- sum(predicted_class == 0 & y_true == 0, na.rm = TRUE)
  false_negative <- sum(predicted_class == 0 & y_true == 1, na.rm = TRUE)

  sensitivity <- ifelse(
    true_positive + false_negative == 0,
    NA_real_,
    true_positive / (true_positive + false_negative)
  )

  specificity <- ifelse(
    true_negative + false_positive == 0,
    NA_real_,
    true_negative / (true_negative + false_positive)
  )

  ppv <- ifelse(
    true_positive + false_positive == 0,
    NA_real_,
    true_positive / (true_positive + false_positive)
  )

  npv <- ifelse(
    true_negative + false_negative == 0,
    NA_real_,
    true_negative / (true_negative + false_negative)
  )

  tibble(
    threshold = threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    accuracy = mean(predicted_class == y_true, na.rm = TRUE),
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    true_positive = true_positive,
    false_positive = false_positive,
    true_negative = true_negative,
    false_negative = false_negative
  )
}

calc_youden_metrics <- function(y_true, pred) {
  if (length(unique(y_true[!is.na(y_true)])) < 2) {
    return(
      tibble(
        threshold = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        ppv = NA_real_,
        npv = NA_real_,
        accuracy = NA_real_,
        balanced_accuracy = NA_real_,
        true_positive = NA_integer_,
        false_positive = NA_integer_,
        true_negative = NA_integer_,
        false_negative = NA_integer_
      )
    )
  }

  roc_object <- pROC::roc(
    response = y_true,
    predictor = pred,
    direction = "<",
    quiet = TRUE
  )

  best <- pROC::coords(
    roc_object,
    x = "best",
    best.method = "youden",
    ret = "threshold",
    transpose = FALSE
  )

  calc_metrics_at_threshold(
    y_true = y_true,
    pred = pred,
    threshold = as.numeric(best$threshold)
  )
}

calc_performance <- function(data, group_label, fixed_threshold) {
  y_true <- data$outcome
  pred <- data$pred

  if (length(unique(y_true[!is.na(y_true)])) < 2) {
    return(
      tibble(
        subgroup = group_label,
        n = nrow(data),
        n_mace = sum(y_true == 1, na.rm = TRUE),
        n_no_mace = sum(y_true == 0, na.rm = TRUE),
        mace_prevalence = mean(y_true == 1, na.rm = TRUE),
        auc = NA_real_,
        auc_ci_low = NA_real_,
        auc_ci_high = NA_real_,
        auprc = NA_real_
      ) %>%
        bind_cols(calc_calibration(y_true, pred)) %>%
        bind_cols(
          calc_metrics_at_threshold(
            y_true,
            pred,
            fixed_threshold
          ) %>%
            rename_with(~ paste0("fixed_", .x))
        ) %>%
        bind_cols(
          calc_youden_metrics(y_true, pred) %>%
            rename_with(~ paste0("subgroup_youden_", .x))
        )
    )
  }

  roc_object <- pROC::roc(
    response = y_true,
    predictor = pred,
    direction = "<",
    quiet = TRUE
  )

  auc_ci <- as.numeric(pROC::ci.auc(roc_object))

  tibble(
    subgroup = group_label,
    n = nrow(data),
    n_mace = sum(y_true == 1, na.rm = TRUE),
    n_no_mace = sum(y_true == 0, na.rm = TRUE),
    mace_prevalence = mean(y_true == 1, na.rm = TRUE),
    auc = as.numeric(pROC::auc(roc_object)),
    auc_ci_low = auc_ci[1],
    auc_ci_high = auc_ci[3],
    auprc = calc_auprc(y_true, pred)
  ) %>%
    bind_cols(calc_calibration(y_true, pred)) %>%
    bind_cols(
      calc_metrics_at_threshold(
        y_true,
        pred,
        fixed_threshold
      ) %>%
        rename_with(~ paste0("fixed_", .x))
    ) %>%
    bind_cols(
      calc_youden_metrics(y_true, pred) %>%
        rename_with(~ paste0("subgroup_youden_", .x))
    )
}

exact_binomial_ci <- function(events, total, conf_level = 0.95) {
  if (total == 0) {
    return(tibble(rate = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }

  result <- binom.test(events, total, conf.level = conf_level)

  tibble(
    rate = events / total,
    ci_low = unname(result$conf.int[1]),
    ci_high = unname(result$conf.int[2])
  )
}

make_global_percentile_score <- function(x) {
  if (sum(!is.na(x)) <= 1) {
    return(rep(50, length(x)))
  }

  100 * (
    rank(x, ties.method = "average", na.last = "keep") - 1
  ) / (
    sum(!is.na(x)) - 1
  )
}

classify_feature_modality <- function(feature_name) {
  case_when(
    str_starts(feature_name, "lipo_") ~ "Lipoproteins",
    str_starts(feature_name, "lipid_") ~ "Lipids",
    str_starts(feature_name, "sm_") ~ "Small molecules",
    TRUE ~ "Expanded clinical"
  )
}

clean_feature_label <- function(feature_name) {
  feature_name %>%
    str_remove("^lipo_") %>%
    str_remove("^lipid_") %>%
    str_remove("^sm_") %>%
    str_replace_all("\\.", " ") %>%
    str_squish()
}

# Load data and locked candidate outputs

if (!file.exists(compiled_cohort_file)) {
  stop("Compiled cohort file not found: ", compiled_cohort_file)
}

if (!file.exists(candidate_output_file)) {
  stop("Script 05 candidate output file not found: ", candidate_output_file)
}

model_source <- read.csv(
  compiled_cohort_file,
  check.names = FALSE
) %>%
  as_tibble() %>%
  mutate(
    mace = normalise_binary(mace),
    masld_binary = normalise_binary(masld)
  ) %>%
  filter(!is.na(mace))

if (any(!model_source$masld_binary %in% c(0, 1, NA))) {
  stop("MASLD could not be consistently converted to binary 0/1 coding.")
}

candidate_outputs <- readRDS(candidate_output_file)

required_candidate_objects <- c(
  "predictions",
  "performance_summary",
  "settings",
  "feature_votes"
)

missing_candidate_objects <- setdiff(
  required_candidate_objects,
  names(candidate_outputs)
)

if (length(missing_candidate_objects) > 0) {
  stop(
    "Script 05 output is missing: ",
    paste(missing_candidate_objects, collapse = ", ")
  )
}

candidate_predictions <- candidate_outputs$predictions %>%
  as_tibble() %>%
  mutate(
    validation_fold = as.integer(validation_fold),
    outcome = as.numeric(outcome),
    pred = as.numeric(pred)
  )

candidate_threshold <- as.numeric(
  candidate_outputs$performance_summary$threshold[1]
)

if (!is.finite(candidate_threshold)) {
  stop("A valid overall candidate threshold was not found in Script 05 outputs.")
}

# Reconstruct original patient indices for the pooled outer-fold predictions
#
# Script 03 stores validation predictions in fold order but does not retain the
# original row index. Recreating the folds with the original seed gives the
# original indices. The outcome vector is then audited before MASLD is attached.

set.seed(original_outer_seed)

reconstructed_validation_folds <- caret::createFolds(
  factor(model_source$mace),
  k = original_outer_folds,
  returnTrain = FALSE
)

prediction_index_map <- map_dfr(
  seq_along(reconstructed_validation_folds),
  function(fold_i) {
    tibble(
      validation_fold = fold_i,
      fold_row_order = seq_along(reconstructed_validation_folds[[fold_i]]),
      cohort_row_index = reconstructed_validation_folds[[fold_i]]
    )
  }
)

candidate_predictions_indexed <- candidate_predictions %>%
  group_by(validation_fold) %>%
  mutate(fold_row_order = row_number()) %>%
  ungroup() %>%
  left_join(
    prediction_index_map,
    by = c("validation_fold", "fold_row_order")
  )

if (any(is.na(candidate_predictions_indexed$cohort_row_index))) {
  stop(
    "At least one prediction could not be mapped to the reconstructed ",
    "validation folds."
  )
}

candidate_predictions_indexed <- candidate_predictions_indexed %>%
  mutate(
    reconstructed_outcome =
      model_source$mace[cohort_row_index],
    masld_binary =
      model_source$masld_binary[cohort_row_index],
    masld_group = case_when(
      masld_binary == 1 ~ "MASLD",
      masld_binary == 0 ~ "Non-MASLD",
      TRUE ~ "MASLD missing"
    )
  )

outcome_mismatch_n <- sum(
  candidate_predictions_indexed$outcome !=
    candidate_predictions_indexed$reconstructed_outcome,
  na.rm = TRUE
)

fold_reconstruction_audit <- candidate_predictions_indexed %>%
  group_by(validation_fold) %>%
  summarise(
    predictions_n = n(),
    prediction_mace_n = sum(outcome == 1),
    reconstructed_mace_n = sum(reconstructed_outcome == 1),
    outcome_mismatch_n = sum(outcome != reconstructed_outcome),
    .groups = "drop"
  )

save_table(
  fold_reconstruction_audit,
  "07A_fold_reconstruction_audit.csv"
)

if (outcome_mismatch_n > 0) {
  stop(
    "Fold reconstruction failed: ",
    outcome_mismatch_n,
    " prediction outcomes did not match the compiled cohort. ",
    "Do not continue until the Script 03 seed and input row order are confirmed."
  )
}

analysis_predictions <- candidate_predictions_indexed %>%
  filter(masld_group != "MASLD missing") %>%
  mutate(
    subgroup = factor(
      masld_group,
      levels = c("Non-MASLD", "MASLD")
    )
  )

missing_masld_summary <- tibble(
  total_predictions = nrow(candidate_predictions_indexed),
  masld_missing_n = sum(
    candidate_predictions_indexed$masld_group == "MASLD missing"
  ),
  masld_complete_n = nrow(analysis_predictions)
)

save_table(
  missing_masld_summary,
  "07A_masld_prediction_coverage.csv"
)

# Part A: locked-model subgroup evaluation

overall_performance <- calc_performance(
  data = analysis_predictions,
  group_label = "Overall MASLD-complete cohort",
  fixed_threshold = candidate_threshold
)

masld_performance <- calc_performance(
  data = analysis_predictions %>% filter(subgroup == "MASLD"),
  group_label = "MASLD",
  fixed_threshold = candidate_threshold
)

non_masld_performance <- calc_performance(
  data = analysis_predictions %>% filter(subgroup == "Non-MASLD"),
  group_label = "Non-MASLD",
  fixed_threshold = candidate_threshold
)

subgroup_performance <- bind_rows(
  overall_performance,
  masld_performance,
  non_masld_performance
)

save_table(
  subgroup_performance,
  "07B_locked_candidate_subgroup_performance.csv"
)

subgroup_performance_presentation <- subgroup_performance %>%
  transmute(
    Cohort = subgroup,
    `Patients, n` = n,
    `MACE events, n` = n_mace,
    `MACE prevalence, %` = round(100 * mace_prevalence, 1),
    `AUC` = round(auc, 3),
    `AUC 95% CI` = paste0(
      round(auc_ci_low, 3),
      "–",
      round(auc_ci_high, 3)
    ),
    `AUPRC` = round(auprc, 3),
    `Brier score` = round(brier_score, 3),
    `Calibration intercept` = round(calibration_intercept, 3),
    `Calibration slope` = round(calibration_slope, 3),
    `Overall threshold` = round(fixed_threshold, 3),
    `Sensitivity at overall threshold` = round(fixed_sensitivity, 3),
    `Specificity at overall threshold` = round(fixed_specificity, 3),
    `PPV at overall threshold` = round(fixed_ppv, 3),
    `NPV at overall threshold` = round(fixed_npv, 3),
    `Subgroup Youden threshold` = round(subgroup_youden_threshold, 3)
  )

save_table(
  subgroup_performance_presentation,
  "07B_locked_candidate_subgroup_performance_presentation.csv"
)

# Compare AUCs between independent MASLD and non-MASLD groups
#
# This is an exploratory test of whether discrimination differs between groups.
# It is not a test of whether either subgroup AUC is clinically acceptable.

masld_data <- analysis_predictions %>%
  filter(subgroup == "MASLD")

non_masld_data <- analysis_predictions %>%
  filter(subgroup == "Non-MASLD")

auc_difference_test <- tryCatch(
  {
    roc_masld <- pROC::roc(
      masld_data$outcome,
      masld_data$pred,
      direction = "<",
      quiet = TRUE
    )

    roc_non_masld <- pROC::roc(
      non_masld_data$outcome,
      non_masld_data$pred,
      direction = "<",
      quiet = TRUE
    )

    test <- pROC::roc.test(
      roc_masld,
      roc_non_masld,
      method = "delong",
      paired = FALSE
    )

    tibble(
      auc_masld = as.numeric(pROC::auc(roc_masld)),
      auc_non_masld = as.numeric(pROC::auc(roc_non_masld)),
      delta_auc_masld_minus_non_masld =
        auc_masld - auc_non_masld,
      p_value = test$p.value,
      interpretation =
        "Exploratory unpaired DeLong comparison between subgroups"
    )
  },
  error = function(e) {
    tibble(
      auc_masld = NA_real_,
      auc_non_masld = NA_real_,
      delta_auc_masld_minus_non_masld = NA_real_,
      p_value = NA_real_,
      interpretation = paste(
        "AUC comparison unavailable:",
        conditionMessage(e)
      )
    )
  }
)

save_table(
  auc_difference_test,
  "07C_masld_vs_non_masld_auc_difference.csv"
)

# Model-score interaction with MASLD status
#
# The interaction tests whether the association between the locked model score
# and MACE differs by MASLD status. It is complementary to the subgroup AUC
# comparison and should be interpreted cautiously in this sample size.

interaction_data <- analysis_predictions %>%
  mutate(
    pred_clipped = pmin(pmax(pred, 1e-6), 1 - 1e-6),
    logit_pred = qlogis(pred_clipped),
    masld_binary = as.numeric(subgroup == "MASLD")
  )

score_interaction_fit <- tryCatch(
  glm(
    outcome ~ logit_pred * masld_binary,
    data = interaction_data,
    family = binomial
  ),
  error = function(e) NULL
)

score_interaction_summary <- if (is.null(score_interaction_fit)) {
  tibble(
    term = "logit_pred:masld_binary",
    estimate = NA_real_,
    std.error = NA_real_,
    statistic = NA_real_,
    p.value = NA_real_,
    odds_ratio = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_
  )
} else {
  broom::tidy(
    score_interaction_fit,
    conf.int = TRUE
  ) %>%
    filter(term == "logit_pred:masld_binary") %>%
    mutate(
      odds_ratio = exp(estimate),
      conf.low = exp(conf.low),
      conf.high = exp(conf.high)
    )
}

save_table(
  score_interaction_summary,
  "07C_model_score_by_masld_interaction.csv"
)

# Calibration groups by MASLD status

make_calibration_groups <- function(data) {
  number_of_groups <- min(
    maximum_calibration_groups,
    max(
      2,
      floor(nrow(data) / minimum_patients_per_calibration_group)
    )
  )

  data %>%
    mutate(
      calibration_group = ntile(pred, number_of_groups)
    ) %>%
    group_by(subgroup, calibration_group) %>%
    summarise(
      n = n(),
      mean_predicted_score = mean(pred),
      median_predicted_score = median(pred),
      observed_mace_n = sum(outcome == 1),
      observed_mace_rate = mean(outcome == 1),
      .groups = "drop"
    ) %>%
    mutate(
      ci_low = case_when(
        observed_mace_n == 0 ~ 0,
        TRUE ~ qbeta(
          0.025,
          observed_mace_n,
          n - observed_mace_n + 1
        )
      ),
      ci_high = case_when(
        observed_mace_n == n ~ 1,
        TRUE ~ qbeta(
          0.975,
          observed_mace_n + 1,
          n - observed_mace_n
        )
      )
    )
}

calibration_by_subgroup <- analysis_predictions %>%
  group_split(subgroup) %>%
  map_dfr(make_calibration_groups)

save_table(
  calibration_by_subgroup,
  "07D_calibration_groups_by_masld.csv"
)

calibration_plot <- ggplot(
  calibration_by_subgroup,
  aes(
    x = mean_predicted_score,
    y = observed_mace_rate,
    shape = subgroup
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    width = 0.015
  ) +
  geom_point(size = 3) +
  facet_wrap(~ subgroup) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  labs(
    title = "Calibration of the locked candidate by MASLD status",
    subtitle = "Pooled outer-fold out-of-sample predictions",
    x = "Mean predicted score",
    y = "Observed MACE proportion",
    shape = NULL,
    caption = paste0(
      "Vertical bars are exact binomial 95% confidence intervals. ",
      "Sparse subgroup estimates should be interpreted cautiously."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none"
  )

save_plot(
  calibration_plot,
  "07D_calibration_by_masld.png",
  width = 10,
  height = 5
)

# ROC curves by MASLD status

roc_curve_data <- analysis_predictions %>%
  group_by(subgroup) %>%
  group_modify(
    ~ {
      if (n_distinct(.x$outcome) < 2) {
        return(tibble())
      }

      roc_object <- pROC::roc(
        .x$outcome,
        .x$pred,
        direction = "<",
        quiet = TRUE
      )

      tibble(
        specificity = roc_object$specificities,
        sensitivity = roc_object$sensitivities,
        false_positive_rate = 1 - specificity
      )
    }
  ) %>%
  ungroup() %>%
  left_join(
    subgroup_performance %>%
      filter(subgroup %in% c("MASLD", "Non-MASLD")) %>%
      select(subgroup, auc),
    by = "subgroup"
  ) %>%
  mutate(
    curve_label = paste0(
      subgroup,
      " | AUC ",
      round(auc, 3)
    )
  )

save_table(
  roc_curve_data,
  "07E_roc_curve_data_by_masld.csv"
)

roc_plot <- ggplot(
  roc_curve_data,
  aes(
    x = false_positive_rate,
    y = sensitivity,
    linetype = curve_label
  )
) +
  geom_line(linewidth = 1) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dotted"
  ) +
  coord_equal() +
  labs(
    title = "Candidate-model discrimination by MASLD status",
    subtitle = "Locked pooled outer-fold predictions",
    x = "False positive rate",
    y = "Sensitivity",
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

save_plot(
  roc_plot,
  "07E_roc_curves_by_masld.png",
  width = 8,
  height = 6
)

# Predicted-score distribution by MASLD status and outcome

score_distribution_plot <- ggplot(
  analysis_predictions,
  aes(
    x = subgroup,
    y = pred,
    shape = factor(
      outcome,
      levels = c(0, 1),
      labels = c("No MACE", "MACE")
    )
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.25
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    alpha = 0.45
  ) +
  geom_jitter(
    width = 0.08,
    alpha = 0.55,
    size = 1.6
  ) +
  geom_hline(
    yintercept = candidate_threshold,
    linetype = "dashed"
  ) +
  labs(
    title = "Out-of-sample candidate scores by MASLD status",
    subtitle = paste0(
      "Dashed line: overall pooled Youden threshold = ",
      round(candidate_threshold, 3)
    ),
    x = NULL,
    y = "Out-of-fold predicted score",
    shape = "Observed outcome"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

save_plot(
  score_distribution_plot,
  "07F_candidate_score_distribution_by_masld.png",
  width = 8,
  height = 6
)

# Risk-score and risk-class stratification
#
# Prefer Script 06 patient scores because they preserve the exact overall
# percentile-score framework. If unavailable, reconstruct the same overall
# percentile score from the locked candidate predictions.

if (file.exists(risk_score_output_file)) {
  risk_score_outputs <- readRDS(risk_score_output_file)

  risk_score_data <- risk_score_outputs$patient_scores %>%
    as_tibble()

  if (nrow(risk_score_data) != nrow(candidate_predictions)) {
    stop(
      "Script 06 patient-score rows do not match Script 05 predictions."
    )
  }

  risk_score_indexed <- risk_score_data %>%
    select(
      validation_fold,
      outcome,
      pred,
      risk_score_100_raw,
      risk_score_100,
      risk_class,
      above_youden
    ) %>%
    group_by(validation_fold) %>%
    mutate(fold_row_order = row_number()) %>%
    ungroup() %>%
    left_join(
      prediction_index_map,
      by = c("validation_fold", "fold_row_order")
    ) %>%
    mutate(
      masld_binary =
        model_source$masld_binary[cohort_row_index],
      subgroup = case_when(
        masld_binary == 1 ~ "MASLD",
        masld_binary == 0 ~ "Non-MASLD",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(subgroup))
} else {
  warning(
    "Script 06 output was not found. Reconstructing the overall ",
    "percentile-score framework directly from Script 05 predictions."
  )

  risk_score_indexed <- analysis_predictions %>%
    mutate(
      risk_score_100_raw =
        make_global_percentile_score(pred),
      risk_score_100 =
        round(risk_score_100_raw),
      above_youden =
        pred >= candidate_threshold,
      risk_class = case_when(
        above_youden ~ "Higher predicted risk",
        risk_score_100_raw <= 33 ~ "Lower predicted risk",
        TRUE ~ "Intermediate predicted risk"
      )
    )
}

risk_score_indexed <- risk_score_indexed %>%
  mutate(
    risk_class = factor(
      risk_class,
      levels = c(
        "Lower predicted risk",
        "Intermediate predicted risk",
        "Higher predicted risk"
      )
    )
  )

risk_class_by_masld <- risk_score_indexed %>%
  group_by(subgroup, risk_class) %>%
  summarise(
    n = n(),
    observed_mace_n = sum(outcome == 1),
    observed_mace_rate = mean(outcome == 1),
    median_score = median(risk_score_100_raw),
    .groups = "drop"
  ) %>%
  group_by(subgroup) %>%
  mutate(
    subgroup_n = sum(n),
    subgroup_proportion = n / subgroup_n
  ) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    ci = list(
      exact_binomial_ci(
        observed_mace_n,
        n
      )
    )
  ) %>%
  unnest(ci) %>%
  ungroup() %>%
  select(-rate)

save_table(
  risk_class_by_masld,
  "07G_risk_classes_by_masld.csv"
)

risk_class_composition_plot <- ggplot(
  risk_class_by_masld,
  aes(
    x = subgroup,
    y = 100 * subgroup_proportion,
    fill = risk_class
  )
) +
  geom_col(position = "stack") +
  geom_text(
    aes(
      label = paste0(
        n,
        "\n",
        round(100 * subgroup_proportion, 1),
        "%"
      )
    ),
    position = position_stack(vjust = 0.5),
    size = 3.5
  ) +
  labs(
    title = "Predicted-risk class composition by MASLD status",
    subtitle = "Classes use the overall cohort score and overall Youden boundary",
    x = NULL,
    y = "Participants (%)",
    fill = "Predicted-risk class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

save_plot(
  risk_class_composition_plot,
  "07G_risk_class_composition_by_masld.png",
  width = 8,
  height = 6
)

risk_class_event_plot <- ggplot(
  risk_class_by_masld,
  aes(
    x = risk_class,
    y = 100 * observed_mace_rate,
    shape = subgroup,
    group = subgroup
  )
) +
  geom_errorbar(
    aes(
      ymin = 100 * ci_low,
      ymax = 100 * ci_high
    ),
    position = position_dodge(width = 0.35),
    width = 0.15
  ) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 3
  ) +
  labs(
    title = "Observed MACE rate within risk classes",
    subtitle = "Comparison of MASLD and non-MASLD participants",
    x = NULL,
    y = "Observed MACE rate (%)",
    shape = NULL,
    caption = "Exact binomial 95% confidence intervals"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

save_plot(
  risk_class_event_plot,
  "07G_observed_mace_by_risk_class_and_masld.png",
  width = 9,
  height = 6
)

# Part B: exploratory MASLD-only nested model development

omics_prefixes <- c(
  "lipo_",
  "lipid_",
  "sm_"
)

metadata_pattern <- paste0(
  "Visit|Timepoint|Label|sampleID|plateID|tubeLabel|",
  "sourceID|Batch|Location|Sample|Study\\.ID|Study.ID"
)

clinical_expanded_vars <- c(
  "age",
  "gender",
  "ethnicity",
  "bmi",
  "smoking_bin",
  "alcohol_frequency",
  "alcohol_drinks",
  "alcohol_dependency",
  "diabetes_bin",
  "hypertension",
  "dyslipidemia",
  "family_history_cardiac",
  "portal_hypertension",
  "cirrhosis_stage",
  "thrombocytopenia",
  "tips",
  "decomp_ascites",
  "decomp_encephalopathy",
  "decomp_variceal_bleed",
  "decomp_child_pugh7",
  "primary_aetiology",
  "med_metformin",
  "med_statin",
  "med_ppi",
  "med_glp1_gip",
  "med_none"
) %>%
  intersect(names(model_source))

lipoprotein_vars <- names(model_source) %>%
  str_subset("^lipo_") %>%
  str_subset(metadata_pattern, negate = TRUE)

lipid_vars <- names(model_source) %>%
  str_subset("^lipid_") %>%
  str_subset(metadata_pattern, negate = TRUE)

small_molecule_vars <- names(model_source) %>%
  str_subset("^sm_") %>%
  str_subset(metadata_pattern, negate = TRUE)

candidate_omics_vars <- unique(
  c(
    lipoprotein_vars,
    lipid_vars,
    small_molecule_vars
  )
)

get_omics_features <- function(x_train) {
  colnames(x_train)[
    str_detect(
      colnames(x_train),
      paste0(
        "^(",
        paste(omics_prefixes, collapse = "|"),
        ")"
      )
    )
  ]
}

get_clinical_features <- function(x_train) {
  setdiff(colnames(x_train), get_omics_features(x_train))
}

select_auc_features <- function(x_train, y_train, top_k) {
  omics_columns <- get_omics_features(x_train)

  if (length(omics_columns) == 0 || is.na(top_k)) {
    return(character())
  }

  auc_table <- map_dfr(
    omics_columns,
    function(feature_name) {
      feature_values <- x_train[, feature_name]

      if (
        length(
          unique(feature_values[!is.na(feature_values)])
        ) < 2
      ) {
        return(
          tibble(
            feature = feature_name,
            auc_strength = NA_real_
          )
        )
      }

      auc_value <- tryCatch(
        as.numeric(
          pROC::auc(
            pROC::roc(
              y_train,
              feature_values,
              quiet = TRUE
            )
          )
        ),
        error = function(e) NA_real_
      )

      tibble(
        feature = feature_name,
        auc_strength = ifelse(
          is.na(auc_value),
          NA_real_,
          max(auc_value, 1 - auc_value)
        )
      )
    }
  )

  auc_table %>%
    filter(!is.na(auc_strength)) %>%
    arrange(desc(auc_strength)) %>%
    slice_head(n = min(top_k, nrow(.))) %>%
    pull(feature)
}

prep_train_test <- function(train_df, test_df) {
  y_train <- train_df$outcome
  y_test <- test_df$outcome

  x_train_raw <- train_df %>% select(-outcome)
  x_test_raw <- test_df %>% select(-outcome)

  numeric_variables <- names(x_train_raw)[
    map_lgl(x_train_raw, is.numeric)
  ]

  categorical_variables <- setdiff(
    names(x_train_raw),
    numeric_variables
  )

  # Numeric preprocessing
  if (length(numeric_variables) > 0) {
    x_train_numeric <- x_train_raw %>%
      select(all_of(numeric_variables))

    x_test_numeric <- x_test_raw %>%
      select(all_of(numeric_variables))

    training_medians <- apply(
      x_train_numeric,
      2,
      median,
      na.rm = TRUE
    )

    training_medians[is.na(training_medians)] <- 0

    for (variable_name in numeric_variables) {
      x_train_numeric[[variable_name]][
        is.na(x_train_numeric[[variable_name]])
      ] <- training_medians[[variable_name]]

      x_test_numeric[[variable_name]][
        is.na(x_test_numeric[[variable_name]])
      ] <- training_medians[[variable_name]]
    }
  } else {
    x_train_numeric <- tibble(.rows = nrow(x_train_raw))
    x_test_numeric <- tibble(.rows = nrow(x_test_raw))
  }

  # Categorical preprocessing
  if (length(categorical_variables) > 0) {
    x_train_categorical <- x_train_raw %>%
      select(all_of(categorical_variables)) %>%
      mutate(
        across(
          everything(),
          ~ replace_na(as.character(.x), "Missing")
        )
      )

    x_test_categorical <- x_test_raw %>%
      select(all_of(categorical_variables)) %>%
      mutate(
        across(
          everything(),
          ~ replace_na(as.character(.x), "Missing")
        )
      )

    for (variable_name in categorical_variables) {
      training_levels <- unique(
        c(
          x_train_categorical[[variable_name]],
          "Other"
        )
      )

      x_test_categorical[[variable_name]][
        !x_test_categorical[[variable_name]] %in%
          training_levels
      ] <- "Other"

      x_train_categorical[[variable_name]] <- factor(
        x_train_categorical[[variable_name]],
        levels = training_levels
      )

      x_test_categorical[[variable_name]] <- factor(
        x_test_categorical[[variable_name]],
        levels = training_levels
      )
    }

    x_train_categorical <- model.matrix(
      ~ . - 1,
      data = x_train_categorical
    ) %>%
      as.data.frame(check.names = FALSE)

    x_test_categorical <- model.matrix(
      ~ . - 1,
      data = x_test_categorical
    ) %>%
      as.data.frame(check.names = FALSE)

    missing_test_columns <- setdiff(
      names(x_train_categorical),
      names(x_test_categorical)
    )

    extra_test_columns <- setdiff(
      names(x_test_categorical),
      names(x_train_categorical)
    )

    for (missing_column in missing_test_columns) {
      x_test_categorical[[missing_column]] <- 0
    }

    if (length(extra_test_columns) > 0) {
      x_test_categorical <- x_test_categorical %>%
        select(-all_of(extra_test_columns))
    }

    x_test_categorical <- x_test_categorical[
      ,
      names(x_train_categorical),
      drop = FALSE
    ]
  } else {
    x_train_categorical <- tibble(.rows = nrow(x_train_raw))
    x_test_categorical <- tibble(.rows = nrow(x_test_raw))
  }

  x_train <- bind_cols(
    x_train_numeric,
    x_train_categorical
  )

  x_test <- bind_cols(
    x_test_numeric,
    x_test_categorical
  )

  training_standard_deviations <- apply(
    x_train,
    2,
    sd,
    na.rm = TRUE
  )

  retained_variables <- names(
    training_standard_deviations
  )[
    is.finite(training_standard_deviations) &
      training_standard_deviations > 0
  ]

  x_train <- x_train[
    ,
    retained_variables,
    drop = FALSE
  ]

  x_test <- x_test[
    ,
    retained_variables,
    drop = FALSE
  ]

  training_means <- colMeans(x_train)
  training_standard_deviations <- apply(
    x_train,
    2,
    sd
  )

  list(
    x_train = as.matrix(
      scale(
        x_train,
        center = training_means,
        scale = training_standard_deviations
      )
    ),
    x_test = as.matrix(
      scale(
        x_test,
        center = training_means,
        scale = training_standard_deviations
      )
    ),
    y_train = y_train,
    y_test = y_test
  )
}

fit_predict_enet <- function(
  x_train,
  y_train,
  x_test,
  alpha_value,
  internal_folds
) {
  weights <- make_weights(y_train)

  cv_fit <- cv.glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = alpha_value,
    weights = weights,
    nfolds = internal_folds,
    type.measure = "auc",
    maxit = 200000
  )

  final_fit <- glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = alpha_value,
    weights = weights,
    lambda = cv_fit$lambda.min,
    maxit = 200000
  )

  predictions <- as.numeric(
    predict(
      final_fit,
      newx = x_test,
      type = "response"
    )
  )

  list(
    fit = final_fit,
    pred = predictions,
    lambda_min = cv_fit$lambda.min
  )
}

safe_fold_count <- function(y, requested_folds) {
  smallest_class <- min(table(y))

  max(
    2,
    min(
      requested_folds,
      as.integer(smallest_class)
    )
  )
}

run_masld_only_nested_model <- function(data) {
  masld_model_data <- data %>%
    filter(masld_binary == 1) %>%
    select(
      mace,
      all_of(clinical_expanded_vars),
      all_of(candidate_omics_vars)
    ) %>%
    mutate(outcome = as.numeric(mace)) %>%
    select(-mace)

  if (n_distinct(masld_model_data$outcome) < 2) {
    stop("The MASLD-only cohort does not contain both MACE outcome classes.")
  }

  outer_k <- safe_fold_count(
    masld_model_data$outcome,
    masld_outer_folds
  )

  set.seed(123)

  outer_folds <- caret::createFolds(
    factor(masld_model_data$outcome),
    k = outer_k,
    returnTrain = FALSE
  )

  fold_results <- map_dfr(
    seq_along(outer_folds),
    function(outer_i) {
      message(
        "Exploratory MASLD-only model: outer fold ",
        outer_i,
        "/",
        length(outer_folds)
      )

      validation_indices <- outer_folds[[outer_i]]

      development_data <- masld_model_data[
        -validation_indices,
        ,
        drop = FALSE
      ]

      validation_data <- masld_model_data[
        validation_indices,
        ,
        drop = FALSE
      ]

      inner_k <- safe_fold_count(
        development_data$outcome,
        masld_inner_folds
      )

      inner_folds <- caret::createFolds(
        factor(development_data$outcome),
        k = inner_k,
        returnTrain = FALSE
      )

      tuning_results <- crossing(
        top_k = top_k_grid,
        alpha = alpha_grid
      ) %>%
        mutate(
          mean_inner_auc = map2_dbl(
            top_k,
            alpha,
            function(current_top_k, current_alpha) {
              inner_auc_values <- map_dbl(
                seq_along(inner_folds),
                function(inner_i) {
                  inner_validation_indices <-
                    inner_folds[[inner_i]]

                  inner_training_data <- development_data[
                    -inner_validation_indices,
                    ,
                    drop = FALSE
                  ]

                  inner_validation_data <- development_data[
                    inner_validation_indices,
                    ,
                    drop = FALSE
                  ]

                  processed_inner <- prep_train_test(
                    inner_training_data,
                    inner_validation_data
                  )

                  clinical_features <- get_clinical_features(
                    processed_inner$x_train
                  )

                  omics_features <- select_auc_features(
                    processed_inner$x_train,
                    processed_inner$y_train,
                    current_top_k
                  )

                  selected_features <- unique(
                    c(
                      clinical_features,
                      omics_features
                    )
                  )

                  if (
                    length(selected_features) < 2 ||
                      n_distinct(processed_inner$y_test) < 2
                  ) {
                    return(NA_real_)
                  }

                  fit_result <- tryCatch(
                    fit_predict_enet(
                      x_train =
                        processed_inner$x_train[
                          ,
                          selected_features,
                          drop = FALSE
                        ],
                      y_train =
                        processed_inner$y_train,
                      x_test =
                        processed_inner$x_test[
                          ,
                          selected_features,
                          drop = FALSE
                        ],
                      alpha_value =
                        current_alpha,
                      internal_folds =
                        safe_fold_count(
                          processed_inner$y_train,
                          inner_k
                        )
                    ),
                    error = function(e) NULL
                  )

                  if (is.null(fit_result)) {
                    return(NA_real_)
                  }

                  as.numeric(
                    pROC::auc(
                      pROC::roc(
                        processed_inner$y_test,
                        fit_result$pred,
                        direction = "<",
                        quiet = TRUE
                      )
                    )
                  )
                }
              )

              if (all(is.na(inner_auc_values))) {
                return(NA_real_)
              }

              mean(inner_auc_values, na.rm = TRUE)
            }
          )
        ) %>%
        filter(!is.na(mean_inner_auc)) %>%
        arrange(
          desc(mean_inner_auc),
          top_k,
          desc(alpha)
        )

      if (nrow(tuning_results) == 0) {
        stop(
          "No valid hyperparameter combination in MASLD outer fold ",
          outer_i
        )
      }

      best_settings <- tuning_results %>%
        slice(1)

      processed_final <- prep_train_test(
        development_data,
        validation_data
      )

      final_clinical_features <- get_clinical_features(
        processed_final$x_train
      )

      final_omics_features <- select_auc_features(
        processed_final$x_train,
        processed_final$y_train,
        best_settings$top_k
      )

      final_selected_features <- unique(
        c(
          final_clinical_features,
          final_omics_features
        )
      )

      final_fit <- fit_predict_enet(
        x_train =
          processed_final$x_train[
            ,
            final_selected_features,
            drop = FALSE
          ],
        y_train =
          processed_final$y_train,
        x_test =
          processed_final$x_test[
            ,
            final_selected_features,
            drop = FALSE
          ],
        alpha_value =
          best_settings$alpha,
        internal_folds =
          safe_fold_count(
            processed_final$y_train,
            inner_k
          )
      )

      coefficient_matrix <- as.matrix(
        coef(final_fit$fit)
      )

      selected_feature_table <- tibble(
        feature = rownames(coefficient_matrix),
        coefficient = as.numeric(
          coefficient_matrix[, 1]
        )
      ) %>%
        filter(
          feature != "(Intercept)",
          coefficient != 0
        ) %>%
        mutate(
          validation_fold = outer_i,
          modality = classify_feature_modality(feature),
          direction = if_else(
            coefficient > 0,
            "Increase risk",
            "Decrease risk"
          )
        )

      validation_outcome <- processed_final$y_test
      validation_prediction <- final_fit$pred

      fold_auc <- if (
        n_distinct(validation_outcome) == 2
      ) {
        as.numeric(
          pROC::auc(
            pROC::roc(
              validation_outcome,
              validation_prediction,
              direction = "<",
              quiet = TRUE
            )
          )
        )
      } else {
        NA_real_
      }

      tibble(
        validation_fold = outer_i,
        final_top_k = best_settings$top_k,
        final_alpha = best_settings$alpha,
        mean_inner_auc = best_settings$mean_inner_auc,
        validation_auc = fold_auc,
        validation_auprc = calc_auprc(
          validation_outcome,
          validation_prediction
        ),
        n_development = nrow(development_data),
        n_validation = nrow(validation_data),
        n_validation_mace =
          sum(validation_outcome == 1),
        n_selected_features =
          nrow(selected_feature_table),
        predictions = list(
          tibble(
            validation_fold = outer_i,
            outcome = validation_outcome,
            pred = validation_prediction
          )
        ),
        selected_features = list(
          selected_feature_table
        )
      )
    }
  )

  predictions <- fold_results %>%
    select(predictions) %>%
    unnest(predictions)

  selected_features <- fold_results %>%
    select(selected_features) %>%
    unnest(selected_features)

  feature_votes <- selected_features %>%
    group_by(feature, modality) %>%
    summarise(
      selected_n = n(),
      selection_frequency =
        selected_n / n_distinct(fold_results$validation_fold),
      mean_coefficient = mean(coefficient),
      median_coefficient = median(coefficient),
      mean_abs_coefficient = mean(abs(coefficient)),
      direction = case_when(
        mean_coefficient > 0 ~ "Increase risk",
        mean_coefficient < 0 ~ "Decrease risk",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      feature_label = clean_feature_label(feature)
    ) %>%
    arrange(
      desc(selection_frequency),
      desc(mean_abs_coefficient)
    )

  pooled_performance <- calc_performance(
    data = predictions,
    group_label = "Exploratory MASLD-only nested model",
    fixed_threshold = candidate_threshold
  ) %>%
    mutate(
      mean_fold_auc =
        mean(fold_results$validation_auc, na.rm = TRUE),
      sd_fold_auc =
        sd(fold_results$validation_auc, na.rm = TRUE),
      mean_fold_auprc =
        mean(fold_results$validation_auprc, na.rm = TRUE)
    )

  list(
    fold_results = fold_results %>%
      select(-predictions, -selected_features),
    predictions = predictions,
    selected_features = selected_features,
    feature_votes = feature_votes,
    pooled_performance = pooled_performance
  )
}

masld_exploratory_results <- NULL

if (run_exploratory_masld_model) {
  masld_exploratory_results <- run_masld_only_nested_model(
    model_source
  )

  save_table(
    masld_exploratory_results$fold_results,
    "07H_masld_only_nested_fold_performance.csv"
  )

  save_table(
    masld_exploratory_results$pooled_performance,
    "07H_masld_only_nested_pooled_performance.csv"
  )

  save_table(
    masld_exploratory_results$feature_votes,
    "07I_masld_only_feature_stability.csv"
  )

  # Compare exploratory MASLD-only feature selection with the original candidate.
  overall_feature_votes <- candidate_outputs$feature_votes %>%
    as_tibble() %>%
    select(
      feature,
      overall_selection_frequency = selection_frequency,
      overall_mean_coefficient = mean_coefficient,
      overall_mean_abs_coefficient = mean_abs_coefficient,
      overall_direction = direction
    )

  feature_stability_comparison <-
    masld_exploratory_results$feature_votes %>%
    select(
      feature,
      modality,
      feature_label,
      masld_selection_frequency = selection_frequency,
      masld_mean_coefficient = mean_coefficient,
      masld_mean_abs_coefficient = mean_abs_coefficient,
      masld_direction = direction
    ) %>%
    full_join(
      overall_feature_votes,
      by = "feature"
    ) %>%
    mutate(
      modality = coalesce(
        modality,
        classify_feature_modality(feature)
      ),
      feature_label = coalesce(
        feature_label,
        clean_feature_label(feature)
      ),
      masld_selection_frequency =
        replace_na(masld_selection_frequency, 0),
      overall_selection_frequency =
        replace_na(overall_selection_frequency, 0),
      selection_frequency_difference =
        masld_selection_frequency -
        overall_selection_frequency
    ) %>%
    arrange(
      desc(
        pmax(
          masld_selection_frequency,
          overall_selection_frequency
        )
      ),
      desc(abs(selection_frequency_difference))
    )

  save_table(
    feature_stability_comparison,
    "07I_masld_vs_overall_feature_stability.csv"
  )

  top_feature_comparison <- feature_stability_comparison %>%
    slice_head(n = 25) %>%
    select(
      feature,
      feature_label,
      modality,
      overall_selection_frequency,
      masld_selection_frequency
    ) %>%
    pivot_longer(
      cols = c(
        overall_selection_frequency,
        masld_selection_frequency
      ),
      names_to = "analysis",
      values_to = "selection_frequency"
    ) %>%
    mutate(
      analysis = recode(
        analysis,
        overall_selection_frequency =
          "Original full-cohort candidate",
        masld_selection_frequency =
          "Exploratory MASLD-only model"
      ),
      feature_label = fct_reorder(
        feature_label,
        selection_frequency,
        .fun = max
      )
    )

  feature_stability_plot <- ggplot(
    top_feature_comparison,
    aes(
      x = selection_frequency,
      y = feature_label,
      shape = analysis
    )
  ) +
    geom_point(size = 3) +
    facet_wrap(
      ~ modality,
      scales = "free_y"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    ) +
    labs(
      title = "Feature-selection stability: overall versus MASLD-only",
      subtitle = paste0(
        "MASLD-only model is exploratory and independently retuned"
      ),
      x = "Outer-fold selection frequency",
      y = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom"
    )

  save_plot(
    feature_stability_plot,
    "07I_feature_stability_overall_vs_masld.png",
    width = 11,
    height = 9
  )
}

# Save combined output object

masld_stratification_outputs <- list(
  settings = list(
    candidate_model = candidate_model,
    candidate_label = candidate_label,
    candidate_threshold = candidate_threshold,
    original_outer_seed = original_outer_seed,
    original_outer_folds = original_outer_folds,
    run_exploratory_masld_model =
      run_exploratory_masld_model,
    subgroup_evaluation =
      "Locked pooled outer-fold predictions; no subgroup retraining",
    exploratory_model_interpretation =
      "MASLD-only nested model; exploratory and not external validation"
  ),
  fold_reconstruction_audit =
    fold_reconstruction_audit,
  indexed_predictions =
    analysis_predictions,
  subgroup_performance =
    subgroup_performance,
  auc_difference_test =
    auc_difference_test,
  score_interaction_summary =
    score_interaction_summary,
  calibration_by_subgroup =
    calibration_by_subgroup,
  risk_class_by_masld =
    risk_class_by_masld,
  exploratory_masld_model =
    masld_exploratory_results
)

saveRDS(
  masld_stratification_outputs,
  file.path(
    output_model_dir,
    "07_masld_stratification_outputs.rds"
  )
)

message("Script 07 MASLD stratification complete.")
message("Locked candidate: ", candidate_label)
message("Overall candidate threshold: ", round(candidate_threshold, 4))
message("Outcome mismatches after fold reconstruction: ", outcome_mismatch_n)
message("Subgroup performance:")
print(subgroup_performance_presentation)

if (run_exploratory_masld_model) {
  message("Exploratory MASLD-only nested-model performance:")
  print(masld_exploratory_results$pooled_performance)
}
