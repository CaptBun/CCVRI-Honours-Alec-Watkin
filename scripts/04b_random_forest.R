# Repository stage: rf
# Source supplied as: 04_random_forest_complete_modality_no_aspirin(1).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# Random-forest complete-modality model using leakage-resistant nested cross-validation
#
# Project:
# Cardiovascular Risk Stratification in Liver Cirrhosis: Multiomics integration
#
# Purpose:
# This script evaluates one complete multimodality random-forest model containing:
#   - expanded clinical predictors with aspirin excluded;
#   - lipoproteins;
#   - lipids;
#   - small molecules;
#   - cytokines;
#   - bile acids; and
#   - pathology.
#
# The validation structure mirrors the elastic-net comparison pipeline as closely
# as possible while using random-forest-specific tuning parameters. A global
# 5-fold outer validation loop estimates held-out performance. Within each outer
# development set, a 5-fold inner cross-validation loop selects the number of
# omics features retained by univariate AUC filtering, mtry, and minimum node size.
# All outcome-informed filtering and preprocessing are learned from training data
# only. The held-out outer fold is used once, after tuning, for final evaluation.
# Load packages
library(tidyverse)
library(caret)
library(pROC)
library(PRROC)
library(ranger)

set.seed(123)

# Reproducibility and computational settings
# ranger parallelises tree construction internally. Because this script runs only
# one model specification, model-level parallelism is unnecessary.
max_rf_threads <- ccrvi_config$max_rf_threads
rf_threads <- min(
  max_rf_threads,
  max(1, parallel::detectCores(logical = TRUE) - 1)
)

# A smaller forest is sufficient while comparing tuning combinations. The final
# held-out model uses more trees to reduce Monte Carlo variation in predictions.
tuning_num_trees <- 500
final_num_trees <- 2000

# Create output folders
# Fold-level checkpoints allow the analysis to resume after interruption without
# rerunning outer folds that have already completed successfully.
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/models", showWarnings = FALSE, recursive = TRUE)
dir.create(
  "outputs/models/checkpoints_rf_complete_no_aspirin",
  showWarnings = FALSE,
  recursive = TRUE
)

checkpoint_dir <- ccrvi_checkpoint_directory("outputs/models/checkpoints_rf_complete_no_aspirin")

# Load data
# The compiled file is expected to contain one row per participant, the binary
# MACE outcome, clinical variables, and all omics variables with modality-specific
# prefixes.
model_source <- read.csv(ccrvi_config$files$compiled_cohort)

# Standardise outcome coding
# MACE is represented as 0/1, with 1 indicating a MACE-positive participant.
model_source <- model_source %>%
  mutate(mace = as.numeric(mace)) %>%
  filter(!is.na(mace))

if (!all(model_source$mace %in% c(0, 1))) {
  stop("The mace outcome must contain only 0 and 1 after numeric conversion.")
}

# Define expanded clinical predictors
# This matches the expanded clinical variable set used in the aspirin-ablation
# elastic-net analysis before aspirin is removed below.
clinical_expanded_vars <- c(
  "age",
  "gender",
  "ethnicity",
  "bmi",
  "smoking_bin",
  "alcohol_frequency",
  "alcohol_drinks",
  "alcohol_dependency",
  "masld",
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
  "med_aspirin",
  "med_ppi",
  "med_glp1_gip",
  "med_none"
)

clinical_expanded_vars <- intersect(
  clinical_expanded_vars,
  names(model_source)
)

# Remove aspirin only
# This retains the other medication indicators while excluding med_aspirin because
# aspirin may partly reflect treatment initiated after a cardiovascular event.
aspirin_var <- "med_aspirin"

if (!aspirin_var %in% names(model_source)) {
  stop("The expected aspirin variable 'med_aspirin' is not present in the cohort.")
}

clinical_expanded_no_aspirin_vars <- setdiff(
  clinical_expanded_vars,
  aspirin_var
)

# Define all omics modalities
# Prefixes are inherited from the compiled multimodal cohort.
omics_definitions <- tribble(
  ~modality,         ~prefix,
  "lipoproteins",    "lipo_",
  "lipids",          "lipid_",
  "small_molecules", "sm_",
  "cytokines",       "cyto_",
  "bile_acids",      "ba_",
  "pathology",       "path_"
)

omics_prefixes <- omics_definitions$prefix

# Extract biological variables for each modality
# Sample-handling and assay metadata are excluded because they should not be used
# as biological predictors of MACE.
omics_list <- omics_definitions %>%
  mutate(
    variables = map(prefix, function(prefix_value) {
      names(model_source) %>%
        str_subset(paste0("^", prefix_value)) %>%
        str_subset(
          "Visit|Timepoint|Label|sampleID|plateID|tubeLabel|sourceID|Batch|Location|Sample|Study\\.ID|Study.ID",
          negate = TRUE
        )
    })
  ) %>%
  select(modality, variables) %>%
  deframe()

# Confirm that every expected modality contributes at least one predictor
missing_modalities <- names(omics_list)[map_int(omics_list, length) == 0]

if (length(missing_modalities) > 0) {
  stop(
    "No predictors were detected for the following modalities: ",
    paste(missing_modalities, collapse = ", ")
  )
}

# Define the single complete-modality model
# No alternative modality combinations are generated in this script.
model_name <- "clinical_expanded_no_aspirin_all_modalities_rf"
model_omics_modalities <- names(omics_list)
model_omics_vars <- unlist(omics_list, use.names = FALSE)
model_predictors <- unique(c(
  clinical_expanded_no_aspirin_vars,
  model_omics_vars
))

# Record the number of variables available before fold-specific preprocessing
variable_count_table <- tibble(
  component = c(
    "clinical_expanded_no_aspirin",
    names(omics_list),
    "complete_model_total"
  ),
  n_variables = c(
    length(clinical_expanded_no_aspirin_vars),
    map_int(omics_list, length),
    length(model_predictors)
  )
)

write.csv(
  variable_count_table,
  "outputs/tables/rf_complete_no_aspirin_variable_counts.csv",
  row.names = FALSE
)

# Random-forest tuning grid
# top_k controls the number of omics features retained by univariate AUC filtering.
# Clinical variables are always retained and are not outcome-filtered.
#
# mtry_rule is converted to an integer after the fold-specific predictor count is
# known. The default random-forest behaviour is approximately sqrt(p); quarter
# and half allow the forest to consider broader correlated feature subsets.
#
# min_node_size controls terminal-node regularisation. Larger values produce
# smoother trees and may reduce overfitting in this small cohort.
top_k_grid <- ccrvi_config$top_k_grid
mtry_rule_grid <- c("sqrt", "quarter", "half")
min_node_size_grid <- c(1, 5, 10)

tuning_grid <- crossing(
  top_k = top_k_grid,
  mtry_rule = mtry_rule_grid,
  min_node_size = min_node_size_grid
)

# Helper: class weights
# The positive MACE class receives the same inverse-frequency weighting logic used
# in the elastic-net pipeline. ranger expects weights aligned to factor levels.
make_class_weights <- function(y) {
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)

  if (n_pos == 0 || n_neg == 0) {
    return(c("0" = 1, "1" = 1))
  }

  c(
    "0" = 1,
    "1" = n_neg / n_pos
  )
}

# Helper: prepare the raw model dataframe
# Missingness is deliberately retained here because all imputation parameters must
# be learned from fold-specific training data.
prepare_model_data <- function(data, predictors) {
  data %>%
    select(mace, all_of(predictors)) %>%
    mutate(outcome = as.numeric(mace)) %>%
    select(-mace)
}

# Helper: leakage-resistant preprocessing
# Numeric values are median-imputed using training data only. Categorical values
# are converted to dummy variables using training-derived levels. Zero-variance
# variables are removed based only on the training data.
#
# Unlike elastic net, random forests do not require centring or scaling because
# tree splits depend on variable ordering rather than Euclidean scale.
prep_train_test_rf <- function(train_df, test_df) {
  y_train <- train_df$outcome
  y_test <- test_df$outcome

  x_train_raw <- train_df %>% select(-outcome)
  x_test_raw <- test_df %>% select(-outcome)

  numeric_vars <- names(x_train_raw)[map_lgl(x_train_raw, is.numeric)]
  categorical_vars <- setdiff(names(x_train_raw), numeric_vars)

  # Numeric preprocessing
  if (length(numeric_vars) > 0) {
    x_train_num <- x_train_raw %>% select(all_of(numeric_vars))
    x_test_num <- x_test_raw %>% select(all_of(numeric_vars))

    train_medians <- vapply(
      x_train_num,
      function(x) median(x, na.rm = TRUE),
      numeric(1)
    )
    train_medians[!is.finite(train_medians)] <- 0

    for (nm in numeric_vars) {
      x_train_num[[nm]][is.na(x_train_num[[nm]])] <- train_medians[[nm]]
      x_test_num[[nm]][is.na(x_test_num[[nm]])] <- train_medians[[nm]]
    }
  } else {
    x_train_num <- tibble(.rows = nrow(x_train_raw))
    x_test_num <- tibble(.rows = nrow(x_test_raw))
  }

  # Categorical preprocessing
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

    x_train_cat <- model.matrix(
      ~ . - 1,
      data = x_train_cat
    ) %>%
      as.data.frame(check.names = FALSE)

    x_test_cat <- model.matrix(
      ~ . - 1,
      data = x_test_cat
    ) %>%
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

  # Remove zero-variance variables based on training data only
  train_sds <- vapply(x_train, sd, numeric(1), na.rm = TRUE)
  keep_vars <- names(train_sds)[is.finite(train_sds) & train_sds > 0]

  if (length(keep_vars) == 0) {
    stop("No non-zero-variance predictors remained after preprocessing.")
  }

  x_train <- x_train[, keep_vars, drop = FALSE]
  x_test <- x_test[, keep_vars, drop = FALSE]

  list(
    x_train = as.data.frame(x_train, check.names = FALSE),
    x_test = as.data.frame(x_test, check.names = FALSE),
    y_train = y_train,
    y_test = y_test,
    feature_names = keep_vars
  )
}

# Helper: identify omics and clinical columns after dummy encoding
get_omics_features <- function(x_train, prefixes = omics_prefixes) {
  colnames(x_train)[
    str_detect(
      colnames(x_train),
      paste0("^(", paste(prefixes, collapse = "|"), ")")
    )
  ]
}

get_clinical_features <- function(x_train, prefixes = omics_prefixes) {
  colnames(x_train)[
    !str_detect(
      colnames(x_train),
      paste0("^(", paste(prefixes, collapse = "|"), ")")
    )
  ]
}

# Helper: rank omics predictors by training-fold univariate AUC
# AUC strength is max(AUC, 1 - AUC), so both positive and inverse associations can
# be retained. This function is called only on training data within an inner or
# outer development fold.
rank_auc_features <- function(x_train, y_train, prefixes = omics_prefixes) {
  omics_cols <- get_omics_features(x_train, prefixes)

  if (length(omics_cols) == 0) {
    return(tibble(
      feature = character(),
      auc_value = numeric(),
      auc_strength = numeric()
    ))
  }

  map_dfr(omics_cols, function(feature_name) {
    feature_values <- x_train[[feature_name]]

    if (length(unique(feature_values[!is.na(feature_values)])) < 2) {
      return(tibble(
        feature = feature_name,
        auc_value = NA_real_,
        auc_strength = NA_real_
      ))
    }

    auc_val <- tryCatch(
      as.numeric(
        pROC::auc(
          pROC::roc(
            response = y_train,
            predictor = feature_values,
            quiet = TRUE
          )
        )
      ),
      error = function(e) NA_real_
    )

    tibble(
      feature = feature_name,
      auc_value = auc_val,
      auc_strength = ifelse(
        is.na(auc_val),
        NA_real_,
        max(auc_val, 1 - auc_val)
      )
    )
  }) %>%
    filter(!is.na(auc_strength)) %>%
    arrange(desc(auc_strength))
}

# Helper: convert an mtry rule into an integer for the current feature count
resolve_mtry <- function(mtry_rule, p) {
  if (p < 1) {
    stop("mtry cannot be resolved because the predictor count is < 1.")
  }

  mtry_value <- case_when(
    mtry_rule == "sqrt" ~ floor(sqrt(p)),
    mtry_rule == "quarter" ~ floor(p / 4),
    mtry_rule == "half" ~ floor(p / 2),
    TRUE ~ NA_real_
  )

  if (is.na(mtry_value)) {
    stop("Unknown mtry rule: ", mtry_rule)
  }

  as.integer(max(1, min(p, mtry_value)))
}

# Helper: fit a weighted probability random forest and generate probabilities
fit_predict_rf <- function(
  x_train,
  y_train,
  x_test,
  mtry_value,
  min_node_size,
  num_trees,
  importance_mode = "none",
  seed_value = 123
) {
  y_factor <- factor(y_train, levels = c(0, 1))
  class_weights <- make_class_weights(y_train)

  rf_fit <- ranger::ranger(
    x = x_train,
    y = y_factor,
    num.trees = num_trees,
    mtry = mtry_value,
    min.node.size = min_node_size,
    probability = TRUE,
    replace = TRUE,
    sample.fraction = 1,
    class.weights = class_weights,
    splitrule = "gini",
    importance = importance_mode,
    write.forest = TRUE,
    num.threads = rf_threads,
    seed = seed_value,
    verbose = FALSE
  )

  pred_matrix <- predict(
    rf_fit,
    data = x_test,
    num.threads = rf_threads
  )$predictions

  if (!"1" %in% colnames(pred_matrix)) {
    stop("The fitted probability forest did not return a probability column for class 1.")
  }

  list(
    pred = as.numeric(pred_matrix[, "1"]),
    fit = rf_fit
  )
}

# Helper: calibration metrics
# The calibration intercept and slope are descriptive in this internal-validation
# setting. Class weighting can improve discrimination while altering probability
# calibration, so Brier score and calibration should be interpreted alongside AUC.
calc_calibration <- function(y_true, pred) {
  pred_clipped <- pmin(pmax(pred, 1e-6), 1 - 1e-6)
  logit_pred <- qlogis(pred_clipped)

  cal_fit <- tryCatch(
    glm(y_true ~ logit_pred, family = binomial),
    error = function(e) NULL
  )

  if (is.null(cal_fit)) {
    return(tibble(
      calibration_intercept = NA_real_,
      calibration_slope = NA_real_,
      brier_score = mean((y_true - pred)^2)
    ))
  }

  tibble(
    calibration_intercept = as.numeric(coef(cal_fit)[1]),
    calibration_slope = as.numeric(coef(cal_fit)[2]),
    brier_score = mean((y_true - pred)^2)
  )
}

# Helper: threshold metrics
# Youden's J is retained as a descriptive threshold to keep reporting compatible
# with the existing project. It is not treated as a clinically validated cutoff.
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

# Helper: AUPRC
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

# Build the single model dataframe
model_data <- prepare_model_data(
  data = model_source,
  predictors = model_predictors
)

# Create global outer validation folds
# These folds are generated with seed 123, matching the uploaded elastic-net script,
# so the same cohort ordering should produce the same five held-out partitions.
validation_folds_global <- caret::createFolds(
  factor(model_data$outcome),
  k = 5,
  returnTrain = FALSE
)

# Helper: validate a completed outer-fold checkpoint
validate_fold_checkpoint <- function(result, validation_i) {
  if (!is.list(result)) {
    return(FALSE)
  }

  required_names <- c(
    "fold_performance",
    "predictions",
    "feature_importance",
    "prefilter_features",
    "tuning_results",
    "final_fit"
  )

  if (!all(required_names %in% names(result))) {
    return(FALSE)
  }

  identical(
    as.integer(result$fold_performance$validation_fold[[1]]),
    as.integer(validation_i)
  )
}

# Run one outer validation fold
run_outer_fold <- function(validation_i) {
  message("Outer validation fold ", validation_i, "/", length(validation_folds_global))

  validation_idx <- validation_folds_global[[validation_i]]
  development_df <- model_data[-validation_idx, , drop = FALSE]
  validation_df <- model_data[validation_idx, , drop = FALSE]

  # Inner 5-fold tuning is performed only inside the outer development set.
  set.seed(123 + validation_i)
  inner_folds <- caret::createFolds(
    factor(development_df$outcome),
    k = 5,
    returnTrain = FALSE
  )

  # Evaluate every tuning combination across the same five inner folds.
  # Performance is aggregated across folds before selecting the final combination,
  # which avoids choosing hyperparameters from a single favourable split.
  inner_fold_results <- map_dfr(seq_along(inner_folds), function(inner_i) {
    message("  Inner fold ", inner_i, "/", length(inner_folds))

    inner_test_idx <- inner_folds[[inner_i]]
    inner_train_df <- development_df[-inner_test_idx, , drop = FALSE]
    inner_test_df <- development_df[inner_test_idx, , drop = FALSE]

    processed_inner <- prep_train_test_rf(
      train_df = inner_train_df,
      test_df = inner_test_df
    )

    clinical_features <- get_clinical_features(processed_inner$x_train)
    auc_ranking <- rank_auc_features(
      x_train = processed_inner$x_train,
      y_train = processed_inner$y_train
    )

    pmap_dfr(
      tuning_grid,
      function(top_k, mtry_rule, min_node_size) {
        omics_features <- auc_ranking %>%
          slice_head(n = min(top_k, nrow(auc_ranking))) %>%
          pull(feature)

        selected_features <- unique(c(clinical_features, omics_features))

        if (length(selected_features) < 2) {
          return(tibble(
            validation_fold = validation_i,
            inner_fold = inner_i,
            top_k = top_k,
            mtry_rule = mtry_rule,
            mtry = NA_integer_,
            min_node_size = min_node_size,
            n_predictors = length(selected_features),
            auc = NA_real_
          ))
        }

        x_train_sel <- processed_inner$x_train[, selected_features, drop = FALSE]
        x_test_sel <- processed_inner$x_test[, selected_features, drop = FALSE]
        mtry_value <- resolve_mtry(mtry_rule, ncol(x_train_sel))

        grid_row <- which(
          tuning_grid$top_k == top_k &
            tuning_grid$mtry_rule == mtry_rule &
            tuning_grid$min_node_size == min_node_size
        )[1]

        seed_value <-
          100000L +
          validation_i * 10000L +
          inner_i * 1000L +
          as.integer(grid_row)

        fit_out <- tryCatch(
          fit_predict_rf(
            x_train = x_train_sel,
            y_train = processed_inner$y_train,
            x_test = x_test_sel,
            mtry_value = mtry_value,
            min_node_size = min_node_size,
            num_trees = tuning_num_trees,
            importance_mode = "none",
            seed_value = seed_value
          ),
          error = function(e) {
            message(
              "    Fit failed: top_k=", top_k,
              ", mtry_rule=", mtry_rule,
              ", min_node_size=", min_node_size,
              ": ", conditionMessage(e)
            )
            NULL
          }
        )

        if (is.null(fit_out)) {
          auc_value <- NA_real_
        } else {
          auc_value <- as.numeric(
            pROC::auc(
              pROC::roc(
                processed_inner$y_test,
                fit_out$pred,
                quiet = TRUE
              )
            )
          )
        }

        tibble(
          validation_fold = validation_i,
          inner_fold = inner_i,
          top_k = top_k,
          mtry_rule = mtry_rule,
          mtry = mtry_value,
          min_node_size = min_node_size,
          n_predictors = ncol(x_train_sel),
          auc = auc_value
        )
      }
    )
  })

  # Aggregate inner-fold performance for each hyperparameter combination.
  tuning_results <- inner_fold_results %>%
    group_by(top_k, mtry_rule, min_node_size) %>%
    summarise(
      mean_inner_auc = mean(auc, na.rm = TRUE),
      sd_inner_auc = sd(auc, na.rm = TRUE),
      n_valid_inner_folds = sum(!is.na(auc)),
      mean_mtry = mean(mtry, na.rm = TRUE),
      mean_predictors = mean(n_predictors, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(validation_fold = validation_i) %>%
    relocate(validation_fold)

  # Require all five inner folds to have produced a valid AUC for a tuning
  # combination. Ties favour lower variability and then fewer prefiltered omics
  # features to retain the most parsimonious equally performing solution.
  best_params <- tuning_results %>%
    filter(n_valid_inner_folds == length(inner_folds)) %>%
    arrange(
      desc(mean_inner_auc),
      sd_inner_auc,
      top_k,
      desc(min_node_size)
    ) %>%
    slice(1)

  if (nrow(best_params) == 0) {
    stop(
      "No hyperparameter combination completed all inner folds in outer fold ",
      validation_i,
      "."
    )
  }

  final_top_k <- best_params$top_k[[1]]
  final_mtry_rule <- best_params$mtry_rule[[1]]
  final_min_node_size <- best_params$min_node_size[[1]]

  # Preprocess the full development set and untouched outer validation fold.
  processed_final <- prep_train_test_rf(
    train_df = development_df,
    test_df = validation_df
  )

  clinical_features <- get_clinical_features(processed_final$x_train)
  auc_ranking_final <- rank_auc_features(
    x_train = processed_final$x_train,
    y_train = processed_final$y_train
  )

  omics_features <- auc_ranking_final %>%
    slice_head(n = min(final_top_k, nrow(auc_ranking_final))) %>%
    pull(feature)

  selected_final <- unique(c(clinical_features, omics_features))

  if (length(selected_final) < 2) {
    stop(
      "Fewer than two predictors remained for the final forest in outer fold ",
      validation_i,
      "."
    )
  }

  x_train_final <- processed_final$x_train[, selected_final, drop = FALSE]
  x_val_final <- processed_final$x_test[, selected_final, drop = FALSE]
  final_mtry <- resolve_mtry(final_mtry_rule, ncol(x_train_final))

  final_fit <- fit_predict_rf(
    x_train = x_train_final,
    y_train = processed_final$y_train,
    x_test = x_val_final,
    mtry_value = final_mtry,
    min_node_size = final_min_node_size,
    num_trees = final_num_trees,
    importance_mode = "permutation",
    seed_value = 500000L + validation_i
  )

  y_val <- processed_final$y_test
  pred_val <- final_fit$pred

  validation_auc <- as.numeric(
    pROC::auc(
      pROC::roc(y_val, pred_val, quiet = TRUE)
    )
  )

  validation_auprc <- calc_auprc(y_val, pred_val)

  # Permutation importance is reported instead of coefficients because random
  # forests do not estimate a single linear effect size for each predictor.
  feature_importance <- enframe(
    final_fit$fit$variable.importance,
    name = "feature",
    value = "permutation_importance"
  ) %>%
    mutate(
      model = model_name,
      validation_fold = validation_i,
      feature_type = if_else(
        str_detect(feature, paste0("^(", paste(omics_prefixes, collapse = "|"), ")")),
        "omics",
        "clinical"
      )
    ) %>%
    arrange(desc(permutation_importance))

  # Record every predictor passed into the final forest. Clinical predictors are
  # always retained; omics predictors are the fold-specific top-k AUC features.
  prefilter_features <- tibble(
    feature = selected_final,
    model = model_name,
    validation_fold = validation_i,
    feature_type = if_else(
      str_detect(selected_final, paste0("^(", paste(omics_prefixes, collapse = "|"), ")")),
      "omics_auc_filtered",
      "clinical_always_retained"
    )
  )

  fold_performance <- tibble(
    model = model_name,
    validation_fold = validation_i,
    omics_modalities = paste(model_omics_modalities, collapse = "+"),
    n_omics_modalities = length(model_omics_modalities),
    final_top_k = final_top_k,
    final_mtry_rule = final_mtry_rule,
    final_mtry = final_mtry,
    final_min_node_size = final_min_node_size,
    tuning_num_trees = tuning_num_trees,
    final_num_trees = final_num_trees,
    n_features_available = ncol(processed_final$x_train),
    n_clinical_features = length(clinical_features),
    n_omics_prefilter_features = length(omics_features),
    n_total_prefilter_features = length(selected_final),
    validation_auc = validation_auc,
    validation_auprc = validation_auprc,
    oob_brier_score = final_fit$fit$prediction.error
  )

  predictions <- tibble(
    model = model_name,
    validation_fold = validation_i,
    row_index = validation_idx,
    outcome = y_val,
    pred = pred_val
  )

  list(
    fold_performance = fold_performance,
    predictions = predictions,
    feature_importance = feature_importance,
    prefilter_features = prefilter_features,
    tuning_results = tuning_results,
    inner_fold_results = inner_fold_results,
    final_fit = final_fit$fit
  )
}

# Run the five outer folds with checkpointing
fold_result_list <- map(
  seq_along(validation_folds_global),
  function(validation_i) {
    checkpoint_file <- file.path(
      checkpoint_dir,
      paste0("outer_fold_", validation_i, ".rds")
    )

    if (file.exists(checkpoint_file)) {
      cached_result <- readRDS(checkpoint_file)

      if (validate_fold_checkpoint(cached_result, validation_i)) {
        message("Loading valid checkpoint for outer fold ", validation_i)
        return(cached_result)
      }

      message("Ignoring invalid checkpoint for outer fold ", validation_i)
    }

    start_time <- Sys.time()
    result <- run_outer_fold(validation_i)
    elapsed_minutes <- as.numeric(
      difftime(Sys.time(), start_time, units = "mins")
    )

    result$fold_performance <- result$fold_performance %>%
      mutate(elapsed_minutes = elapsed_minutes)

    saveRDS(result, checkpoint_file)
    message(
      "Saved outer-fold checkpoint ", validation_i,
      " (", round(elapsed_minutes, 2), " min)"
    )

    result
  }
)

# Combine fold-level outputs
rf_fold_performance <- map_dfr(
  fold_result_list,
  "fold_performance"
)

rf_predictions <- map_dfr(
  fold_result_list,
  "predictions"
) %>%
  arrange(row_index)

rf_feature_importance <- map_dfr(
  fold_result_list,
  "feature_importance"
)

rf_prefilter_features <- map_dfr(
  fold_result_list,
  "prefilter_features"
)

rf_tuning_results <- map_dfr(
  fold_result_list,
  "tuning_results"
)

rf_inner_fold_results <- map_dfr(
  fold_result_list,
  "inner_fold_results"
)

# Summarise outer-fold discrimination
rf_performance_summary <- rf_fold_performance %>%
  summarise(
    model = model_name,
    mean_auc = mean(validation_auc, na.rm = TRUE),
    sd_auc = sd(validation_auc, na.rm = TRUE),
    min_auc = min(validation_auc, na.rm = TRUE),
    max_auc = max(validation_auc, na.rm = TRUE),
    mean_auprc = mean(validation_auprc, na.rm = TRUE),
    sd_auprc = sd(validation_auprc, na.rm = TRUE)
  )

# Calculate pooled out-of-fold performance
# Each participant appears once because predictions come only from their held-out
# outer validation fold.
pooled_roc <- pROC::roc(
  rf_predictions$outcome,
  rf_predictions$pred,
  quiet = TRUE
)

pooled_auc_ci <- pROC::ci.auc(pooled_roc)
pooled_auprc <- calc_auprc(
  rf_predictions$outcome,
  rf_predictions$pred
)

pooled_calibration <- calc_calibration(
  rf_predictions$outcome,
  rf_predictions$pred
)

pooled_threshold <- calc_threshold_metrics(
  rf_predictions$outcome,
  rf_predictions$pred
)

rf_pooled_performance <- tibble(
  model = model_name,
  pooled_auc = as.numeric(pROC::auc(pooled_roc)),
  auc_ci_low = as.numeric(pooled_auc_ci[1]),
  auc_ci_high = as.numeric(pooled_auc_ci[3]),
  pooled_auprc = pooled_auprc
) %>%
  bind_cols(pooled_calibration) %>%
  bind_cols(pooled_threshold)

# Summarise how consistently predictors entered the final forests
# For random forest, inclusion frequency reflects the AUC-prefilter stage rather
# than coefficient-based variable selection. Clinical predictors should usually
# appear in all five folds unless removed as zero-variance after dummy encoding.
rf_prefilter_votes <- rf_prefilter_features %>%
  group_by(feature, feature_type) %>%
  summarise(
    included_n = n_distinct(validation_fold),
    inclusion_frequency = included_n / length(validation_folds_global),
    .groups = "drop"
  ) %>%
  arrange(desc(inclusion_frequency), feature_type, feature)

# Summarise permutation importance across folds
# Importance values can vary in scale across fitted forests, so both the raw mean
# and within-fold rank are reported. Rank is often more interpretable for comparing
# stability across folds than directly averaging unnormalised importance values.
rf_feature_importance_summary <- rf_feature_importance %>%
  group_by(validation_fold) %>%
  mutate(
    importance_rank = min_rank(desc(permutation_importance))
  ) %>%
  ungroup() %>%
  group_by(feature, feature_type) %>%
  summarise(
    folds_present = n_distinct(validation_fold),
    mean_permutation_importance = mean(permutation_importance, na.rm = TRUE),
    median_permutation_importance = median(permutation_importance, na.rm = TRUE),
    mean_importance_rank = mean(importance_rank, na.rm = TRUE),
    best_importance_rank = min(importance_rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_importance_rank, desc(mean_permutation_importance))

# Save analysis objects
rf_complete_modality_outputs <- list(
  model_name = model_name,
  model_omics_modalities = model_omics_modalities,
  clinical_expanded_no_aspirin_vars = clinical_expanded_no_aspirin_vars,
  omics_definitions = omics_definitions,
  omics_list = omics_list,
  variable_count_table = variable_count_table,
  tuning_grid = tuning_grid,
  tuning_num_trees = tuning_num_trees,
  final_num_trees = final_num_trees,
  validation_folds_global = validation_folds_global,
  fold_result_list = fold_result_list,
  rf_fold_performance = rf_fold_performance,
  rf_performance_summary = rf_performance_summary,
  rf_predictions = rf_predictions,
  rf_pooled_performance = rf_pooled_performance,
  rf_feature_importance = rf_feature_importance,
  rf_feature_importance_summary = rf_feature_importance_summary,
  rf_prefilter_features = rf_prefilter_features,
  rf_prefilter_votes = rf_prefilter_votes,
  rf_tuning_results = rf_tuning_results,
  rf_inner_fold_results = rf_inner_fold_results
)

saveRDS(
  rf_complete_modality_outputs,
  "outputs/models/random_forest_complete_modality_no_aspirin_nested_cv_outputs.rds"
)

# Save flat tables for thesis/manuscript review
write.csv(
  rf_fold_performance,
  "outputs/tables/rf_complete_no_aspirin_fold_performance.csv",
  row.names = FALSE
)

write.csv(
  rf_performance_summary,
  "outputs/tables/rf_complete_no_aspirin_performance_summary.csv",
  row.names = FALSE
)

write.csv(
  rf_pooled_performance,
  "outputs/tables/rf_complete_no_aspirin_pooled_performance.csv",
  row.names = FALSE
)

write.csv(
  rf_predictions,
  "outputs/tables/rf_complete_no_aspirin_predictions.csv",
  row.names = FALSE
)

write.csv(
  rf_feature_importance,
  "outputs/tables/rf_complete_no_aspirin_feature_importance_by_fold.csv",
  row.names = FALSE
)

write.csv(
  rf_feature_importance_summary,
  "outputs/tables/rf_complete_no_aspirin_feature_importance_summary.csv",
  row.names = FALSE
)

write.csv(
  rf_prefilter_votes,
  "outputs/tables/rf_complete_no_aspirin_prefilter_votes.csv",
  row.names = FALSE
)

write.csv(
  rf_tuning_results,
  "outputs/tables/rf_complete_no_aspirin_tuning_results.csv",
  row.names = FALSE
)

write.csv(
  rf_inner_fold_results,
  "outputs/tables/rf_complete_no_aspirin_inner_fold_results.csv",
  row.names = FALSE
)

# Print primary results
print(rf_performance_summary)
print(rf_pooled_performance)
print(rf_fold_performance)

message("Random-forest complete-modality analysis complete.")
message(
  "Primary pooled results: outputs/tables/",
  "rf_complete_no_aspirin_pooled_performance.csv"
)
message(
  "Full output object: outputs/models/",
  "random_forest_complete_modality_no_aspirin_nested_cv_outputs.rds"
)
