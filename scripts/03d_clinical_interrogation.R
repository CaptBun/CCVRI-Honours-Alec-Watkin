# Repository stage: clinical
# Source supplied as: 03D_clinical_model_interrogation_and_biomarker_identification(5).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# 03D Clinical model interrogation and clinical biomarker identification
#
# Project:
#   Cardiovascular Risk Stratification in Liver Cirrhosis
#
# Purpose:
#   This script investigates why the expanded clinical model performs better than
#   the compact clinical baseline and identifies candidate clinical biomarkers
#   associated with MACE risk.
#
#   It is deliberately clinical-only. No omics variables are used.
#
# Main questions:
#   1. Does the expanded clinical model outperform the baseline clinical model?
#   2. Which clinical domains contribute most to predictive performance?
#   3. Are medications or liver-severity variables driving the expanded model?
#   4. Which individual clinical variables show independent or univariate
#      discriminatory signal?
#   5. Which variables should be treated as candidate clinical biomarkers for
#      later CCRVI development?
#
# Key outputs:
#   outputs/tables/clinical_interrogation/
#   outputs/figures/clinical_interrogation/
#
# Notes:
#   - Univariate AUC is used for single-variable clinical biomarker screening,
#     matching the feature-screening logic used elsewhere in this project.
#   - Elastic net is used only for multivariable benchmark and ablation models.

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(pROC)
library(PRROC)

set.seed(123)

# Output folders
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables/clinical_interrogation", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures/clinical_interrogation", recursive = TRUE, showWarnings = FALSE)

# Load compiled modelling cohort
model_source <- read.csv(ccrvi_config$files$compiled_cohort) %>%
  mutate(mace = as.numeric(mace)) %>%
  filter(!is.na(mace))

# Helper functions
save_table <- function(x, filename) {
  write.csv(
    x,
    file = file.path("outputs/tables/clinical_interrogation", filename),
    row.names = FALSE
  )
}

save_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path("outputs/figures/clinical_interrogation", filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

make_weights <- function(y) {
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)

  if (n_pos == 0 || n_neg == 0) {
    return(rep(1, length(y)))
  }

  ifelse(y == 1, n_neg / n_pos, 1)
}

calc_auprc <- function(y_true, pred) {
  if (length(unique(y_true)) < 2) {
    return(NA_real_)
  }

  pr_val <- PRROC::pr.curve(
    scores.class0 = pred[y_true == 1],
    scores.class1 = pred[y_true == 0],
    curve = FALSE
  )

  pr_val$auc.integral
}

calc_calibration <- function(y_true, pred) {
  pred_clipped <- pmin(pmax(pred, 1e-6), 1 - 1e-6)
  logit_pred <- qlogis(pred_clipped)

  cal_fit <- tryCatch(
    glm(y_true ~ logit_pred, family = binomial),
    error = function(e) NULL
  )

  if (is.null(cal_fit)) {
    return(
      tibble(
        calibration_intercept = NA_real_,
        calibration_slope = NA_real_,
        brier_score = mean((y_true - pred)^2)
      )
    )
  }

  tibble(
    calibration_intercept = as.numeric(coef(cal_fit)[1]),
    calibration_slope = as.numeric(coef(cal_fit)[2]),
    brier_score = mean((y_true - pred)^2)
  )
}

calc_threshold_metrics <- function(y_true, pred) {
  roc_obj <- pROC::roc(y_true, pred, quiet = TRUE)

  best <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
    transpose = FALSE
  )

  pred_class <- ifelse(pred >= best$threshold, 1, 0)

  tibble(
    threshold = best$threshold,
    sensitivity = best$sensitivity,
    specificity = best$specificity,
    ppv = best$ppv,
    npv = best$npv,
    youden_j = best$sensitivity + best$specificity - 1,
    accuracy = mean(pred_class == y_true),
    balanced_accuracy = mean(c(best$sensitivity, best$specificity))
  )
}

# Preprocess training and validation folds
#
# Numeric variables:
#   Median imputation learned from training data only.
#
# Categorical variables:
#   Missing values are encoded as "Missing".
#   Levels are learned from training data.
#   Unseen validation levels are mapped to "Other".
#   Dummy variables are created with model.matrix.
prep_train_test <- function(train_df, test_df) {
  y_train <- train_df$outcome
  y_test <- test_df$outcome

  x_train_raw <- train_df %>% select(-outcome)
  x_test_raw <- test_df %>% select(-outcome)

  numeric_vars <- names(x_train_raw)[map_lgl(x_train_raw, is.numeric)]
  categorical_vars <- setdiff(names(x_train_raw), numeric_vars)

  if (length(numeric_vars) > 0) {
    x_train_num <- x_train_raw %>% select(all_of(numeric_vars))
    x_test_num <- x_test_raw %>% select(all_of(numeric_vars))

    train_medians <- apply(x_train_num, 2, median, na.rm = TRUE)
    train_medians[is.na(train_medians)] <- 0

    for (nm in numeric_vars) {
      x_train_num[[nm]][is.na(x_train_num[[nm]])] <- train_medians[[nm]]
      x_test_num[[nm]][is.na(x_test_num[[nm]])] <- train_medians[[nm]]
    }
  } else {
    x_train_num <- tibble(.rows = nrow(x_train_raw))
    x_test_num <- tibble(.rows = nrow(x_test_raw))
  }

  if (length(categorical_vars) > 0) {
    x_train_cat <- x_train_raw %>%
      select(all_of(categorical_vars)) %>%
      mutate(across(everything(), ~ replace_na(as.character(.x), "Missing")))

    x_test_cat <- x_test_raw %>%
      select(all_of(categorical_vars)) %>%
      mutate(across(everything(), ~ replace_na(as.character(.x), "Missing")))

    for (nm in categorical_vars) {
      train_levels <- unique(c(x_train_cat[[nm]], "Other"))
      x_test_cat[[nm]][!x_test_cat[[nm]] %in% train_levels] <- "Other"

      x_train_cat[[nm]] <- factor(x_train_cat[[nm]], levels = train_levels)
      x_test_cat[[nm]] <- factor(x_test_cat[[nm]], levels = train_levels)
    }

    x_train_cat <- model.matrix(~ . - 1, data = x_train_cat) %>%
      as.data.frame(check.names = FALSE)

    x_test_cat <- model.matrix(~ . - 1, data = x_test_cat) %>%
      as.data.frame(check.names = FALSE)

    missing_test_cols <- setdiff(names(x_train_cat), names(x_test_cat))
    extra_test_cols <- setdiff(names(x_test_cat), names(x_train_cat))

    for (missing_col in missing_test_cols) {
      x_test_cat[[missing_col]] <- 0
    }

    if (length(extra_test_cols) > 0) {
      x_test_cat <- x_test_cat %>%
        select(-all_of(extra_test_cols))
    }

    x_test_cat <- x_test_cat[, names(x_train_cat), drop = FALSE]
  } else {
    x_train_cat <- tibble(.rows = nrow(x_train_raw))
    x_test_cat <- tibble(.rows = nrow(x_test_raw))
  }

  x_train <- bind_cols(x_train_num, x_train_cat)
  x_test <- bind_cols(x_test_num, x_test_cat)

  train_sds <- apply(x_train, 2, sd, na.rm = TRUE)
  keep_vars <- names(train_sds)[is.finite(train_sds) & train_sds > 0]

  x_train <- x_train[, keep_vars, drop = FALSE]
  x_test <- x_test[, keep_vars, drop = FALSE]

  train_means <- colMeans(x_train)
  train_sds <- apply(x_train, 2, sd)

  list(
    x_train = as.matrix(scale(x_train, center = train_means, scale = train_sds)),
    x_test = as.matrix(scale(x_test, center = train_means, scale = train_sds)),
    y_train = y_train,
    y_test = y_test,
    feature_names = keep_vars
  )
}

fit_predict_enet <- function(x_train, y_train, x_test, alpha_val) {
  w_train <- make_weights(y_train)

  cv_fit <- cv.glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = alpha_val,
    weights = w_train,
    nfolds = 5,
    type.measure = "auc",
    parallel = FALSE,
    maxit = 200000
  )

  final_fit <- glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = alpha_val,
    weights = w_train,
    lambda = cv_fit$lambda.min,
    maxit = 200000
  )

  pred <- as.numeric(
    predict(
      final_fit,
      newx = x_test,
      type = "response"
    )
  )

  list(
    pred = pred,
    fit = final_fit,
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se
  )
}

# Clinical variable definitions
#
# mets_count is intentionally excluded from all analyses.
demographic_vars <- c(
  "age",
  "gender",
  "ethnicity"
)

lifestyle_vars <- c(
  "bmi",
  "smoking_bin",
  "alcohol_frequency",
  "alcohol_drinks",
  "alcohol_dependency"
)

cardiometabolic_vars <- c(
  "masld",
  "diabetes_bin",
  "hypertension",
  "dyslipidemia",
  "family_history_cardiac"
)

liver_severity_vars <- c(
  "portal_hypertension",
  "cirrhosis_stage",
  "thrombocytopenia",
  "tips",
  "decomp_ascites",
  "decomp_encephalopathy",
  "decomp_variceal_bleed",
  "decomp_child_pugh7",
  "primary_aetiology"
)

medication_vars <- c(
  "med_metformin",
  "med_statin",
  "med_aspirin",
  "med_ppi",
  "med_glp1_gip",
  "med_none"
)

clinical_baseline_vars <- c(
  "age",
  "gender",
  "masld",
  "diabetes_bin",
  "hypertension",
  "dyslipidemia",
  "bmi",
  "smoking_bin"
)

clinical_block_list <- list(
  demographic = demographic_vars,
  lifestyle = lifestyle_vars,
  cardiometabolic = cardiometabolic_vars,
  liver_severity = liver_severity_vars,
  medication = medication_vars
)

clinical_block_list <- map(
  clinical_block_list,
  ~ intersect(.x, names(model_source))
)

clinical_baseline_vars <- intersect(clinical_baseline_vars, names(model_source))

clinical_expanded_vars <- clinical_block_list %>%
  unlist(use.names = FALSE) %>%
  unique()

clinical_variable_universe <- unique(c(
  clinical_baseline_vars,
  clinical_expanded_vars
))

clinical_variable_metadata <- imap_dfr(
  clinical_block_list,
  ~ tibble(
    clinical_block = .y,
    variable = .x
  )
) %>%
  mutate(
    is_baseline_variable = variable %in% clinical_baseline_vars,
    is_expanded_variable = variable %in% clinical_expanded_vars
  )

save_table(
  clinical_variable_metadata,
  "03D_clinical_variable_metadata.csv"
)

# Define multivariable interrogation models
#
# Models fall into four families:
#   1. benchmark
#   2. block_ablation
#   3. variable_ablation
#   4. progressive_build
clinical_model_sets <- list(
  clinical_baseline = clinical_baseline_vars,
  clinical_expanded_full = clinical_expanded_vars
)

model_family_lookup <- tibble(
  model = c("clinical_baseline", "clinical_expanded_full"),
  interrogation_family = c("benchmark", "benchmark"),
  ablated_block = NA_character_,
  ablated_variable = NA_character_,
  progressive_step = NA_character_
)

# Block ablation
for (block_name in names(clinical_block_list)) {
  block_vars <- clinical_block_list[[block_name]]
  model_name <- paste0("expanded_minus_", block_name)

  clinical_model_sets[[model_name]] <- setdiff(
    clinical_expanded_vars,
    block_vars
  )

  model_family_lookup <- bind_rows(
    model_family_lookup,
    tibble(
      model = model_name,
      interrogation_family = "block_ablation",
      ablated_block = block_name,
      ablated_variable = NA_character_,
      progressive_step = NA_character_
    )
  )
}

# Individual medication and liver-severity variable ablation
variables_for_individual_ablation <- unique(c(
  intersect(medication_vars, names(model_source)),
  intersect(liver_severity_vars, names(model_source))
))

for (var_name in variables_for_individual_ablation) {
  model_name <- paste0("expanded_minus_", var_name)

  clinical_model_sets[[model_name]] <- setdiff(
    clinical_expanded_vars,
    var_name
  )

  var_block <- clinical_variable_metadata %>%
    filter(variable == var_name) %>%
    pull(clinical_block) %>%
    first()

  model_family_lookup <- bind_rows(
    model_family_lookup,
    tibble(
      model = model_name,
      interrogation_family = "variable_ablation",
      ablated_block = var_block,
      ablated_variable = var_name,
      progressive_step = NA_character_
    )
  )
}

# Progressive build-up
progressive_definitions <- list(
  progressive_baseline_plus_liver = unique(c(
    clinical_baseline_vars,
    clinical_block_list$liver_severity
  )),
  progressive_baseline_plus_medication = unique(c(
    clinical_baseline_vars,
    clinical_block_list$medication
  )),
  progressive_baseline_plus_liver_medication = unique(c(
    clinical_baseline_vars,
    clinical_block_list$liver_severity,
    clinical_block_list$medication
  )),
  progressive_baseline_plus_liver_medication_cardiometabolic = unique(c(
    clinical_baseline_vars,
    clinical_block_list$liver_severity,
    clinical_block_list$medication,
    clinical_block_list$cardiometabolic
  )),
  progressive_full_expanded = clinical_expanded_vars
)

for (model_name in names(progressive_definitions)) {
  clinical_model_sets[[model_name]] <- progressive_definitions[[model_name]]

  model_family_lookup <- bind_rows(
    model_family_lookup,
    tibble(
      model = model_name,
      interrogation_family = "progressive_build",
      ablated_block = NA_character_,
      ablated_variable = NA_character_,
      progressive_step = model_name
    )
  )
}

# Remove duplicate model definitions but keep the first named definition.
model_signature <- map_chr(
  clinical_model_sets,
  ~ paste(sort(.x), collapse = "|")
)

clinical_model_sets <- clinical_model_sets[!duplicated(model_signature)]

model_family_lookup <- model_family_lookup %>%
  filter(model %in% names(clinical_model_sets)) %>%
  distinct(model, .keep_all = TRUE)

clinical_model_definition_table <- imap_dfr(
  clinical_model_sets,
  ~ tibble(
    model = .y,
    n_raw_variables = length(.x),
    variables = paste(.x, collapse = "; ")
  )
) %>%
  left_join(model_family_lookup, by = "model")

save_table(
  clinical_model_definition_table,
  "03D_clinical_interrogation_model_definitions.csv"
)

# Shared validation folds
validation_folds_global <- caret::createFolds(
  factor(model_source$mace),
  k = 5,
  returnTrain = FALSE
)

alpha_grid <- ccrvi_config$alpha_grid

# Run one clinical elastic-net model
run_clinical_model <- function(model_name, predictor_vars) {
  message("Running clinical interrogation model: ", model_name)

  model_data <- model_source %>%
    select(mace, all_of(predictor_vars)) %>%
    mutate(outcome = as.numeric(mace)) %>%
    select(-mace)

  map_dfr(seq_along(validation_folds_global), function(validation_i) {
    message("  Validation fold ", validation_i)

    validation_idx <- validation_folds_global[[validation_i]]

    development_df <- model_data[-validation_idx, ]
    validation_df <- model_data[validation_idx, ]

    dev_outer_folds <- caret::createFolds(
      factor(development_df$outcome),
      k = 5,
      returnTrain = FALSE
    )

    nested_results <- map_dfr(seq_along(dev_outer_folds), function(outer_i) {
      outer_test_idx <- dev_outer_folds[[outer_i]]

      outer_train_df <- development_df[-outer_test_idx, ]
      outer_test_df <- development_df[outer_test_idx, ]

      processed_outer <- prep_train_test(
        train_df = outer_train_df,
        test_df = outer_test_df
      )

      if (ncol(processed_outer$x_train) < 2) {
        return(
          tibble(
            validation_fold = validation_i,
            nested_outer_fold = outer_i,
            model = model_name,
            alpha = NA_real_,
            nested_auc = NA_real_
          )
        )
      }

      tuning_results <- tibble(alpha = alpha_grid) %>%
        mutate(
          nested_auc = map_dbl(alpha, function(alpha_value) {
            fit_out <- tryCatch(
              fit_predict_enet(
                x_train = processed_outer$x_train,
                y_train = processed_outer$y_train,
                x_test = processed_outer$x_test,
                alpha_val = alpha_value
              ),
              error = function(e) NULL
            )

            if (is.null(fit_out)) {
              return(NA_real_)
            }

            as.numeric(
              pROC::auc(
                pROC::roc(
                  response = processed_outer$y_test,
                  predictor = fit_out$pred,
                  quiet = TRUE
                )
              )
            )
          })
        )

      best <- tuning_results %>%
        filter(!is.na(nested_auc)) %>%
        arrange(desc(nested_auc)) %>%
        slice(1)

      if (nrow(best) == 0) {
        return(
          tibble(
            validation_fold = validation_i,
            nested_outer_fold = outer_i,
            model = model_name,
            alpha = NA_real_,
            nested_auc = NA_real_
          )
        )
      }

      tibble(
        validation_fold = validation_i,
        nested_outer_fold = outer_i,
        model = model_name,
        alpha = best$alpha,
        nested_auc = best$nested_auc
      )
    })

    final_params <- nested_results %>%
      filter(!is.na(nested_auc), !is.na(alpha)) %>%
      group_by(alpha) %>%
      summarise(
        mean_nested_auc = mean(nested_auc, na.rm = TRUE),
        sd_nested_auc = sd(nested_auc, na.rm = TRUE),
        n_selected = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_nested_auc), sd_nested_auc, desc(n_selected)) %>%
      slice(1)

    if (nrow(final_params) == 0) {
      stop("No valid alpha selected for ", model_name, " validation fold ", validation_i)
    }

    final_alpha <- final_params$alpha

    processed_final <- prep_train_test(
      train_df = development_df,
      test_df = validation_df
    )

    if (ncol(processed_final$x_train) < 2) {
      stop("Fewer than two model features after preprocessing for ", model_name)
    }

    final_fit <- fit_predict_enet(
      x_train = processed_final$x_train,
      y_train = processed_final$y_train,
      x_test = processed_final$x_test,
      alpha_val = final_alpha
    )

    y_val <- processed_final$y_test
    pred_val <- final_fit$pred

    coefs <- coef(final_fit$fit)

    selected_model_features <- tibble(
      feature = rownames(coefs),
      coefficient = as.numeric(coefs[, 1])
    ) %>%
      filter(feature != "(Intercept)", coefficient != 0) %>%
      mutate(
        direction = case_when(
          coefficient > 0 ~ "Increase risk",
          coefficient < 0 ~ "Decrease risk",
          TRUE ~ NA_character_
        )
      )

    tibble(
      model = model_name,
      validation_fold = validation_i,
      n_raw_variables = length(predictor_vars),
      n_features_available = ncol(processed_final$x_train),
      final_alpha = final_alpha,
      final_lambda_min = final_fit$lambda_min,
      final_lambda_1se = final_fit$lambda_1se,
      n_model_features_selected = nrow(selected_model_features),
      validation_auc = as.numeric(pROC::auc(pROC::roc(y_val, pred_val, quiet = TRUE))),
      validation_auprc = calc_auprc(y_val, pred_val),
      predictions = list(
        tibble(
          model = model_name,
          validation_fold = validation_i,
          outcome = y_val,
          pred = pred_val
        )
      ),
      selected_features = list(
        selected_model_features %>%
          mutate(
            model = model_name,
            validation_fold = validation_i
          )
      ),
      nested_results = list(nested_results)
    )
  })
}

# Run multivariable benchmark, ablation, and progressive models
clinical_interrogation_results <- imap_dfr(
  clinical_model_sets,
  ~ run_clinical_model(
    model_name = .y,
    predictor_vars = .x
  )
)

saveRDS(
  clinical_interrogation_results,
  "outputs/models/clinical_interrogation_results.rds"
)

# Flatten multivariable model outputs
clinical_predictions <- clinical_interrogation_results %>%
  select(predictions) %>%
  tidyr::unnest(cols = c(predictions))

clinical_selected_features <- clinical_interrogation_results %>%
  select(selected_features) %>%
  tidyr::unnest(cols = c(selected_features))

clinical_nested_results <- clinical_interrogation_results %>%
  select(nested_results) %>%
  tidyr::unnest(cols = c(nested_results))

clinical_fold_performance <- clinical_interrogation_results %>%
  select(
    model,
    validation_fold,
    n_raw_variables,
    n_features_available,
    final_alpha,
    final_lambda_min,
    final_lambda_1se,
    n_model_features_selected,
    validation_auc,
    validation_auprc
  ) %>%
  left_join(model_family_lookup, by = "model")

clinical_fold_summary <- clinical_fold_performance %>%
  group_by(model, interrogation_family, ablated_block, ablated_variable, progressive_step) %>%
  summarise(
    mean_auc = mean(validation_auc, na.rm = TRUE),
    sd_auc = sd(validation_auc, na.rm = TRUE),
    min_auc = min(validation_auc, na.rm = TRUE),
    max_auc = max(validation_auc, na.rm = TRUE),
    mean_auprc = mean(validation_auprc, na.rm = TRUE),
    sd_auprc = sd(validation_auprc, na.rm = TRUE),
    mean_features_available = mean(n_features_available, na.rm = TRUE),
    mean_features_selected = mean(n_model_features_selected, na.rm = TRUE),
    .groups = "drop"
  )

pooled_clinical_performance <- clinical_predictions %>%
  group_by(model) %>%
  group_modify(~ {
    roc_obj <- pROC::roc(.x$outcome, .x$pred, quiet = TRUE)
    roc_ci <- pROC::ci.auc(roc_obj)

    pooled_auprc <- calc_auprc(.x$outcome, .x$pred)
    cal <- calc_calibration(.x$outcome, .x$pred)
    thresh <- calc_threshold_metrics(.x$outcome, .x$pred)

    tibble(
      pooled_auc = as.numeric(pROC::auc(roc_obj)),
      auc_ci_low = as.numeric(roc_ci[1]),
      auc_ci_high = as.numeric(roc_ci[3]),
      pooled_auprc = pooled_auprc
    ) %>%
      bind_cols(cal) %>%
      bind_cols(thresh)
  }) %>%
  ungroup() %>%
  left_join(clinical_fold_summary, by = "model") %>%
  arrange(desc(pooled_auc))

full_reference <- pooled_clinical_performance %>%
  filter(model == "clinical_expanded_full") %>%
  slice(1)

baseline_reference <- pooled_clinical_performance %>%
  filter(model == "clinical_baseline") %>%
  slice(1)

full_auc <- full_reference$pooled_auc
full_auprc <- full_reference$pooled_auprc
full_brier <- full_reference$brier_score

baseline_auc <- baseline_reference$pooled_auc
baseline_auprc <- baseline_reference$pooled_auprc

clinical_decision_table <- pooled_clinical_performance %>%
  mutate(
    delta_auc_vs_full_expanded = pooled_auc - full_auc,
    delta_auprc_vs_full_expanded = pooled_auprc - full_auprc,
    delta_brier_vs_full_expanded = brier_score - full_brier,
    delta_auc_vs_clinical_baseline = pooled_auc - baseline_auc,
    delta_auprc_vs_clinical_baseline = pooled_auprc - baseline_auprc
  ) %>%
  arrange(desc(pooled_auc))

clinical_feature_votes <- clinical_selected_features %>%
  group_by(model, feature) %>%
  summarise(
    selected_n = n(),
    selection_frequency = selected_n / 5,
    mean_coefficient = mean(coefficient),
    sd_coefficient = sd(coefficient),
    median_coefficient = median(coefficient),
    mean_abs_coefficient = mean(abs(coefficient)),
    sign_consistency = max(mean(coefficient > 0), mean(coefficient < 0)),
    direction = case_when(
      mean_coefficient > 0 ~ "Increase risk",
      mean_coefficient < 0 ~ "Decrease risk",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  arrange(model, desc(selection_frequency), desc(mean_abs_coefficient))

# Univariate clinical biomarker screening
#
# This section treats each clinical variable as a candidate clinical biomarker.
# It estimates individual discriminatory signal using ROC AUC, matching the
# univariate AUC logic used in the broader feature-selection pipeline.
prepare_univariate_predictor <- function(data, variable) {
  x <- data[[variable]]

  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  x_chr <- as.character(x)

  # Binary factor/character variables can be scored directly as 0/1.
  non_missing_levels <- unique(x_chr[!is.na(x_chr) & x_chr != ""])

  if (length(non_missing_levels) == 2) {
    return(as.numeric(factor(x_chr)))
  }

  # Multi-level categorical variables are represented by predicted probabilities
  # from a univariate logistic model. This gives one risk score per patient while
  # allowing variables such as primary_aetiology or cirrhosis_stage to be ranked
  # alongside binary and numeric variables.
  model_df <- tibble(
    outcome = data$mace,
    x = factor(replace_na(x_chr, "Missing"))
  )

  fit <- tryCatch(
    glm(outcome ~ x, data = model_df, family = binomial),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(rep(NA_real_, nrow(data)))
  }

  as.numeric(predict(fit, type = "response"))
}

calc_univariate_clinical_auc <- function(variable) {
  predictor <- prepare_univariate_predictor(model_source, variable)
  outcome <- model_source$mace

  valid <- !is.na(predictor) & !is.na(outcome)

  if (sum(valid) < 10 || length(unique(outcome[valid])) < 2 || length(unique(predictor[valid])) < 2) {
    return(
      tibble(
        variable = variable,
        n_complete = sum(valid),
        missing_prop = mean(!valid),
        auc = NA_real_,
        auc_ci_low = NA_real_,
        auc_ci_high = NA_real_,
        auc_strength = NA_real_,
        auprc = NA_real_,
        threshold = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        ppv = NA_real_,
        npv = NA_real_,
        youden_j = NA_real_
      )
    )
  }

  roc_obj <- tryCatch(
    pROC::roc(outcome[valid], predictor[valid], quiet = TRUE),
    error = function(e) NULL
  )

  if (is.null(roc_obj)) {
    return(
      tibble(
        variable = variable,
        n_complete = sum(valid),
        missing_prop = mean(!valid),
        auc = NA_real_,
        auc_ci_low = NA_real_,
        auc_ci_high = NA_real_,
        auc_strength = NA_real_,
        auprc = NA_real_,
        threshold = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        ppv = NA_real_,
        npv = NA_real_,
        youden_j = NA_real_
      )
    )
  }

  auc_value <- as.numeric(pROC::auc(roc_obj))
  auc_ci <- as.numeric(pROC::ci.auc(roc_obj))
  auprc_value <- calc_auprc(outcome[valid], predictor[valid])

  best <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
    transpose = FALSE
  )

  tibble(
    variable = variable,
    n_complete = sum(valid),
    missing_prop = mean(!valid),
    auc = auc_value,
    auc_ci_low = auc_ci[1],
    auc_ci_high = auc_ci[3],
    auc_strength = max(auc_value, 1 - auc_value),
    auprc = auprc_value,
    threshold = best$threshold,
    sensitivity = best$sensitivity,
    specificity = best$specificity,
    ppv = best$ppv,
    npv = best$npv,
    youden_j = best$sensitivity + best$specificity - 1
  )
}

clinical_univariate_biomarker_auc <- map_dfr(
  clinical_variable_universe,
  calc_univariate_clinical_auc
) %>%
  left_join(clinical_variable_metadata, by = "variable") %>%
  arrange(desc(auc_strength))

# Optional univariate odds-ratio table
#
# AUC is the main single-variable screening metric. Odds ratios are included as
# supporting clinical interpretation, especially for binary predictors.
calc_univariate_or <- function(variable) {
  model_df <- model_source %>%
    select(mace, all_of(variable)) %>%
    rename(x = all_of(variable)) %>%
    mutate(
      x = if (is.numeric(x)) x else factor(replace_na(as.character(x), "Missing"))
    ) %>%
    filter(!is.na(mace), !is.na(x))

  if (nrow(model_df) < 10 || length(unique(model_df$mace)) < 2 || length(unique(model_df$x)) < 2) {
    return(
      tibble(
        variable = variable,
        term = NA_character_,
        odds_ratio = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_
      )
    )
  }

  fit <- tryCatch(
    glm(mace ~ x, data = model_df, family = binomial),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(
      tibble(
        variable = variable,
        term = NA_character_,
        odds_ratio = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_
      )
    )
  }

  coef_tbl <- summary(fit)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term != "(Intercept)") %>%
    transmute(
      variable = variable,
      term = term,
      odds_ratio = exp(Estimate),
      ci_low = exp(Estimate - 1.96 * `Std. Error`),
      ci_high = exp(Estimate + 1.96 * `Std. Error`),
      p_value = `Pr(>|z|)`
    )

  coef_tbl
}

clinical_univariate_or <- map_dfr(
  clinical_variable_universe,
  calc_univariate_or
)

# Clinical biomarker identification evidence score
#
# This integrates several lines of evidence:
#   1. Univariate discriminatory strength
#   2. Selection stability in expanded clinical elastic-net model
#   3. Performance loss when removed from the expanded clinical model
#   4. Clinical domain relevance
#
# This is not a formal statistical test. It is an evidence-ranking table for
# candidate clinical biomarker identification and later CCRVI interpretation.
make_feature_lookup <- function(features, source_variables) {
  tibble(feature = unique(features)) %>%
    mutate(
      variable = map_chr(feature, function(f) {
        matches <- source_variables[
          startsWith(f, source_variables)
        ]
        
        if (length(matches) == 0) {
          return(NA_character_)
        }
        
        matches[which.max(nchar(matches))]
      })
    )
}

feature_lookup <- make_feature_lookup(
  features = clinical_feature_votes$feature,
  source_variables = clinical_variable_universe
)

expanded_feature_stability <- clinical_feature_votes %>%
  filter(model == "clinical_expanded_full") %>%
  left_join(feature_lookup, by = "feature") %>%
  filter(!is.na(variable)) %>%
  group_by(variable) %>%
  summarise(
    expanded_selection_frequency = max(selection_frequency, na.rm = TRUE),
    expanded_mean_abs_coefficient = max(mean_abs_coefficient, na.rm = TRUE),
    expanded_direction = direction[which.max(mean_abs_coefficient)],
    .groups = "drop"
  )

variable_ablation_effect <- clinical_decision_table %>%
  filter(interrogation_family == "variable_ablation") %>%
  transmute(
    variable = ablated_variable,
    delta_auc_when_removed = delta_auc_vs_full_expanded,
    delta_auprc_when_removed = delta_auprc_vs_full_expanded,
    performance_loss_auc_when_removed = -delta_auc_vs_full_expanded,
    performance_loss_auprc_when_removed = -delta_auprc_vs_full_expanded
  )

clinical_biomarker_identification <- clinical_univariate_biomarker_auc %>%
  left_join(expanded_feature_stability, by = "variable") %>%
  left_join(variable_ablation_effect, by = "variable") %>%
  mutate(
    expanded_selection_frequency = replace_na(expanded_selection_frequency, 0),
    expanded_mean_abs_coefficient = replace_na(expanded_mean_abs_coefficient, 0),
    performance_loss_auc_when_removed = replace_na(performance_loss_auc_when_removed, 0),
    performance_loss_auprc_when_removed = replace_na(performance_loss_auprc_when_removed, 0),

    evidence_high_univariate_auc = auc_strength >= 0.65,
    evidence_selected_stably = expanded_selection_frequency >= 0.60,
    evidence_selected_any = expanded_selection_frequency > 0,
    evidence_performance_loss_when_removed = performance_loss_auc_when_removed >= 0.01,
    evidence_clinically_targeted_domain = clinical_block %in% c("liver_severity", "medication", "cardiometabolic"),

    clinical_biomarker_evidence_score =
      as.integer(evidence_high_univariate_auc) +
      as.integer(evidence_selected_stably) +
      as.integer(evidence_performance_loss_when_removed) +
      as.integer(evidence_clinically_targeted_domain),

    biomarker_priority = case_when(
      clinical_biomarker_evidence_score >= 4 ~ "high_priority_candidate_biomarker",
      clinical_biomarker_evidence_score == 3 ~ "moderate_priority_candidate_biomarker",
      clinical_biomarker_evidence_score == 2 ~ "supporting_candidate_biomarker",
      TRUE ~ "low_priority_or_contextual_variable"
    )
  ) %>%
  arrange(
    desc(clinical_biomarker_evidence_score),
    desc(auc_strength),
    desc(expanded_selection_frequency),
    desc(performance_loss_auc_when_removed)
  )

# Save tables
save_table(clinical_fold_performance, "03D_clinical_fold_performance.csv")
save_table(clinical_fold_summary, "03D_clinical_fold_summary.csv")
save_table(pooled_clinical_performance, "03D_clinical_pooled_performance.csv")
save_table(clinical_decision_table, "03D_clinical_decision_table.csv")
save_table(clinical_nested_results, "03D_clinical_nested_results.csv")
save_table(clinical_feature_votes, "03D_clinical_feature_votes.csv")
save_table(clinical_univariate_biomarker_auc, "03D_clinical_univariate_biomarker_auc.csv")
save_table(clinical_univariate_or, "03D_clinical_univariate_odds_ratios.csv")
save_table(clinical_biomarker_identification, "03D_clinical_biomarker_identification.csv")

block_ablation_summary <- clinical_decision_table %>%
  filter(interrogation_family == "block_ablation") %>%
  transmute(
    removed_block = ablated_block,
    model,
    pooled_auc,
    delta_auc_vs_full_expanded,
    pooled_auprc,
    delta_auprc_vs_full_expanded,
    brier_score,
    delta_brier_vs_full_expanded,
    calibration_slope,
    mean_features_selected
  ) %>%
  arrange(delta_auc_vs_full_expanded)

variable_ablation_summary <- clinical_decision_table %>%
  filter(interrogation_family == "variable_ablation") %>%
  transmute(
    variable = ablated_variable,
    block = ablated_block,
    model,
    pooled_auc,
    delta_auc_vs_full_expanded,
    pooled_auprc,
    delta_auprc_vs_full_expanded,
    brier_score,
    delta_brier_vs_full_expanded,
    calibration_slope,
    mean_features_selected
  ) %>%
  arrange(delta_auc_vs_full_expanded)

progressive_build_summary <- clinical_decision_table %>%
  filter(interrogation_family == "progressive_build" | model == "clinical_baseline") %>%
  mutate(
    progressive_order = case_when(
      model == "clinical_baseline" ~ 1,
      model == "progressive_baseline_plus_liver" ~ 2,
      model == "progressive_baseline_plus_medication" ~ 3,
      model == "progressive_baseline_plus_liver_medication" ~ 4,
      model == "progressive_baseline_plus_liver_medication_cardiometabolic" ~ 5,
      model %in% c("clinical_expanded_full", "progressive_full_expanded") ~ 6,
      TRUE ~ 99
    )
  ) %>%
  arrange(progressive_order) %>%
  mutate(
    incremental_auc_from_previous = pooled_auc - lag(pooled_auc),
    incremental_auprc_from_previous = pooled_auprc - lag(pooled_auprc)
  )

save_table(block_ablation_summary, "03D_block_ablation_summary.csv")
save_table(variable_ablation_summary, "03D_variable_ablation_summary.csv")
save_table(progressive_build_summary, "03D_progressive_build_summary.csv")

# Clinical correlation screening
#
# This creates a simple encoded clinical matrix and calculates Spearman
# correlations. It is intended to identify obvious redundancy, not to replace
# formal causal interpretation.
encoded_clinical_for_correlation <- model_source %>%
  select(all_of(clinical_variable_universe)) %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.factor), ~ as.numeric(.x))) %>%
  mutate(across(everything(), as.numeric))

clinical_correlation_matrix <- cor(
  encoded_clinical_for_correlation,
  use = "pairwise.complete.obs",
  method = "spearman"
)

clinical_correlation_long <- as.data.frame(as.table(clinical_correlation_matrix)) %>%
  as_tibble() %>%
  rename(
    variable_1 = Var1,
    variable_2 = Var2,
    spearman_r = Freq
  ) %>%
  filter(variable_1 != variable_2) %>%
  mutate(abs_spearman_r = abs(spearman_r)) %>%
  arrange(desc(abs_spearman_r))

save_table(clinical_correlation_long, "03D_clinical_correlation_long.csv")

# Figures

# 1. Benchmark and ablation AUCs
clinical_auc_plot <- clinical_decision_table %>%
  filter(interrogation_family %in% c("benchmark", "block_ablation")) %>%
  mutate(model = factor(model, levels = rev(model))) %>%
  ggplot(aes(x = model, y = pooled_auc, fill = interrogation_family)) +
  geom_col() +
  geom_errorbar(aes(ymin = auc_ci_low, ymax = auc_ci_high), width = 0.25) +
  coord_flip() +
  labs(
    title = "Clinical model benchmark and block ablation",
    subtitle = "Pooled out-of-fold AUC with confidence intervals",
    x = NULL,
    y = "Pooled AUC",
    fill = "Model family"
  ) +
  theme_minimal(base_size = 12)

save_plot(clinical_auc_plot, "03D_clinical_block_ablation_auc.png", width = 9, height = 6)

# 2. Block ablation delta AUC
block_delta_plot <- block_ablation_summary %>%
  mutate(removed_block = factor(removed_block, levels = removed_block[order(delta_auc_vs_full_expanded)])) %>%
  ggplot(aes(x = removed_block, y = delta_auc_vs_full_expanded)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Effect of removing each clinical block",
    subtitle = "Negative values indicate performance loss versus full expanded clinical",
    x = "Removed block",
    y = "Delta pooled AUC"
  ) +
  theme_minimal(base_size = 12)

save_plot(block_delta_plot, "03D_block_ablation_delta_auc.png", width = 8, height = 5)

# 3. Individual variable ablation delta AUC
variable_delta_plot <- variable_ablation_summary %>%
  mutate(variable = factor(variable, levels = variable[order(delta_auc_vs_full_expanded)])) %>%
  ggplot(aes(x = variable, y = delta_auc_vs_full_expanded, fill = block)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Effect of removing medication and liver-severity variables",
    subtitle = "Negative values indicate performance loss versus full expanded clinical",
    x = "Removed variable",
    y = "Delta pooled AUC",
    fill = "Clinical block"
  ) +
  theme_minimal(base_size = 12)

save_plot(variable_delta_plot, "03D_variable_ablation_delta_auc.png", width = 9, height = 7)

# 4. Progressive build-up
progressive_plot <- progressive_build_summary %>%
  mutate(
    model_label = case_when(
      model == "clinical_baseline" ~ "Baseline",
      model == "progressive_baseline_plus_liver" ~ "+ Liver severity",
      model == "progressive_baseline_plus_medication" ~ "+ Medication",
      model == "progressive_baseline_plus_liver_medication" ~ "+ Liver + medication",
      model == "progressive_baseline_plus_liver_medication_cardiometabolic" ~ "+ Liver + medication + cardiometabolic",
      model %in% c("clinical_expanded_full", "progressive_full_expanded") ~ "Full expanded",
      TRUE ~ model
    ),
    model_label = factor(model_label, levels = model_label)
  ) %>%
  ggplot(aes(x = model_label, y = pooled_auc, group = 1)) +
  geom_line() +
  geom_point(size = 2.8) +
  coord_flip() +
  labs(
    title = "Progressive clinical feature build-up",
    subtitle = "How discrimination changes as clinical information is added",
    x = NULL,
    y = "Pooled AUC"
  ) +
  theme_minimal(base_size = 12)

save_plot(progressive_plot, "03D_progressive_clinical_build_auc.png", width = 8, height = 5)

# 5. Univariate clinical biomarker AUC ranking
univariate_auc_plot <- clinical_univariate_biomarker_auc %>%
  filter(!is.na(auc_strength)) %>%
  arrange(desc(auc_strength)) %>%
  distinct(variable, .keep_all = TRUE) %>%
  slice_head(n = 25) %>%
  mutate(
    variable = factor(
      variable,
      levels = rev(unique(variable))
    )
  ) %>%
  ggplot(aes(x = variable,
             y = auc_strength,
             fill = clinical_block)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top candidate clinical biomarkers by univariate AUC strength",
    subtitle = "AUC strength = max(AUC, 1 − AUC)",
    x = NULL,
    y = "Univariate AUC strength",
    fill = "Clinical block"
  ) +
  theme_minimal(base_size = 12)

save_plot(univariate_auc_plot, "03D_univariate_clinical_biomarker_auc.png", width = 9, height = 7)

# 6. Clinical biomarker evidence score
biomarker_score_plot <- clinical_biomarker_identification %>%
  slice_head(n = 25) %>%
  mutate(variable = factor(variable, levels = rev(variable))) %>%
  ggplot(aes(x = variable, y = clinical_biomarker_evidence_score, fill = biomarker_priority)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Candidate clinical biomarker identification evidence score",
    subtitle = "Combines univariate AUC, ENet stability, ablation loss, and clinical domain relevance",
    x = NULL,
    y = "Evidence score",
    fill = "Biomarker priority"
  ) +
  theme_minimal(base_size = 12)

save_plot(biomarker_score_plot, "03D_clinical_biomarker_evidence_score.png", width = 9, height = 7)

# 7. Feature stability in expanded clinical model
expanded_stability_plot <- clinical_feature_votes %>%
  filter(model == "clinical_expanded_full") %>%
  slice_max(order_by = selection_frequency, n = 25, with_ties = FALSE) %>%
  mutate(feature = factor(feature, levels = rev(feature))) %>%
  ggplot(aes(x = feature, y = selection_frequency, fill = direction)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Feature-selection stability in expanded clinical model",
    x = NULL,
    y = "Selection frequency across validation folds",
    fill = "Coefficient direction"
  ) +
  theme_minimal(base_size = 12)

save_plot(expanded_stability_plot, "03D_expanded_clinical_feature_stability.png", width = 9, height = 7)

# 8. Clinical correlation heatmap
top_corr_vars <- clinical_biomarker_identification %>%
  slice_head(n = min(25, n())) %>%
  pull(variable)

corr_plot_df <- as.data.frame(as.table(clinical_correlation_matrix[top_corr_vars, top_corr_vars])) %>%
  as_tibble() %>%
  rename(
    variable_1 = Var1,
    variable_2 = Var2,
    spearman_r = Freq
  )

clinical_corr_plot <- ggplot(
  corr_plot_df,
  aes(x = variable_1, y = variable_2, fill = spearman_r)
) +
  geom_tile() +
  coord_equal() +
  labs(
    title = "Correlation among top candidate clinical biomarkers",
    x = NULL,
    y = NULL,
    fill = "Spearman r"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_plot(clinical_corr_plot, "03D_candidate_clinical_biomarker_correlation.png", width = 8, height = 7)

message("Clinical model interrogation complete.")
message("Full expanded clinical pooled AUC: ", round(full_auc, 3))
message("Clinical baseline pooled AUC: ", round(baseline_auc, 3))
message("Clinical biomarker identification table saved as 03D_clinical_biomarker_identification.csv")
