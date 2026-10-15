# Repository stage: sensitivity
# Source supplied as: 03F_medication_sensitivity_from_original_script03_figures_fixed(5).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# Modality comparison using two-layer nested cross-validation
#
# Project:
#   Cardiovascular Risk Stratification in Liver Cirrhosis: Multiomics integration
#
# Purpose:
#   Focused medication sensitivity analysis using the original Script 03 nested-CV
#   pipeline as closely as possible. This compares the current leading candidate
#   model against medication-restricted versions.
#
# Major design principles:
#   1. Validation performance is estimated using an outer 5-fold validation loop.
#   2. Hyperparameters are selected inside the development data only.
#   3. Univariate AUC feature filtering is performed inside each training fold only.
#   4. Final validation-fold predictions are generated only for patients held out
#      from the full model-development process.
#   5. Clinical baseline and expanded clinical versions are run in parallel.
#
# Important interpretation note:
#   This script compares modelling strategies. It should not be interpreted as a
#   final locked clinical model unless the modelling plan, cohort definition,
#   feature filtering strategy, and validation approach have been pre-specified.

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(pROC)
library(PRROC)

# Limit lower-level numerical libraries to one thread per R worker.
# This prevents accidental nested threading when several model-level workers are
# running at once. It usually improves total CPU utilisation and stability.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)


use_model_parallel <- ccrvi_config$parallel_models
# Use model-level parallelism: one worker runs one complete model specification.
# Increase max_model_workers if your machine has enough RAM and you want to use
# more CPU. A value between 6 and 12 is usually reasonable for this workflow.
max_model_workers <- ccrvi_config$max_model_workers
n_model_workers <- min(max_model_workers, max(1, parallel::detectCores(logical = TRUE) - 1))

# Increase the allowed size of objects exported to parallel workers.
# The model source data and model-set definitions are modest, but this prevents
# avoidable future.globals.maxSize errors on some systems.
options(future.globals.maxSize = 8 * 1024^3)

if (use_model_parallel &&
    requireNamespace("future", quietly = TRUE) &&
    requireNamespace("furrr", quietly = TRUE)) {
  future::plan(future::multisession, workers = n_model_workers)
  message("Model-level parallel backend registered with ", n_model_workers, " workers.")
} else {
  use_model_parallel <- FALSE
  message("Running model loop sequentially. Install future and furrr for model-level parallelism.")
}

# Keep glmnet cross-validation sequential inside each model to avoid nested
# parallelism. Nested parallelism often reduces performance and can cause
# difficult-to-debug worker errors.
use_glmnet_parallel <- FALSE

set.seed(123)

# Create output folders
# Outputs are separated from source data so that tables and model objects can be
# reviewed, versioned, and imported into thesis/manuscript figures.
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/models", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/models/checkpoints", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/models/checkpoints_v6_no_mets", showWarnings = FALSE, recursive = TRUE)

# Medication sensitivity outputs are isolated so they do not overwrite the full
# 191-model comparison outputs.
dir.create("outputs/tables/medication_sensitivity_script03", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures/medication_sensitivity_script03", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/models/medication_sensitivity_script03", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/models/checkpoints_medication_sensitivity_script03", showWarnings = FALSE, recursive = TRUE)

# Load data
# This file is expected to be produced by the compilation script. It should
# contain one row per participant, MACE/CVD outcomes, clinical variables, and
# all available omics variables with modality-specific prefixes.
model_source <- read.csv(ccrvi_config$files$compiled_cohort)

# Standardise outcome coding
# glmnet requires a numeric binary outcome for binomial models. Here, MACE is
# coded as 0/1, where 1 indicates that the participant experienced MACE.
model_source <- model_source %>%
  mutate(
    mace = as.numeric(mace)
  ) %>%
  filter(!is.na(mace))

# Define clinical variable sets
#
# The baseline clinical model is intentionally compact. It contains the core
# cardiometabolic and demographic variables that approximate conventional risk
# assessment logic.
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

# The expanded clinical model includes broader clinical context relevant to
# cirrhosis and cardiovascular risk. This includes liver disease severity,
# decompensation history, medication indicators, alcohol variables, and primary
# liver disease aetiology.
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

# Retain only variables that are actually present in the compiled dataset.
# This makes the script robust to minor changes in the clinical extract.
clinical_baseline_vars <- intersect(clinical_baseline_vars, names(model_source))
clinical_expanded_vars <- intersect(clinical_expanded_vars, names(model_source))

# Medication sensitivity variants. These are the only intentional changes to the
# clinical variable set. All preprocessing, feature filtering, tuning, folds,
# weighting, and model fitting below remain as in Script 03.
medication_vars <- c(
  "med_metformin",
  "med_statin",
  "med_aspirin",
  "med_ppi",
  "med_glp1_gip",
  "med_none"
) %>%
  intersect(names(model_source))

target_medication_vars <- c(
  "med_statin",
  "med_aspirin"
) %>%
  intersect(names(model_source))

clinical_expanded_no_all_medications_vars <- setdiff(
  clinical_expanded_vars,
  medication_vars
)

clinical_expanded_no_statin_aspirin_vars <- setdiff(
  clinical_expanded_vars,
  target_medication_vars
)

# Define omics modalities
#
# Each omics modality is identified by the prefix added during compilation.
# Bile acids are included here using the ba_ prefix.
#
# To add another modality later, add one row to this table only.
omics_definitions <- tribble(
  ~modality,          ~prefix,
  "lipoproteins",     "lipo_",
  "lipids",           "lipid_",
  "small_molecules",  "sm_",
  "cytokines",        "cyto_",
  "bile_acids",       "ba_",
  "pathology",        "path_"
)

omics_prefixes <- omics_definitions$prefix

# Extract feature names for each omics modality.
# Metadata columns are excluded because they describe sample handling rather
# than biological measurements and should not enter the predictive model.
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

# Record variable counts for audit purposes.
variable_count_table <- tibble(
  feature_set = c(
    "clinical_baseline",
    "clinical_expanded",
    "clinical_expanded_no_all_medications",
    "clinical_expanded_no_statin_aspirin",
    names(omics_list)
  ),
  n_variables = c(
    length(clinical_baseline_vars),
    length(clinical_expanded_vars),
    length(clinical_expanded_no_all_medications_vars),
    length(clinical_expanded_no_statin_aspirin_vars),
    map_int(omics_list, length)
  )
)

write.csv(
  variable_count_table,
  "outputs/tables/medication_sensitivity_script03/03F_model_comparison_variable_counts.csv",
  row.names = FALSE
)

# Generate focused medication sensitivity model sets
#
# This is intentionally not an exhaustive modality screen. It reuses the original
# Script 03 modelling functions to compare only the current leading candidate
# model and medication-restricted variants.
#
# Candidate modality set:
#   Expanded clinical + Lipids + Small molecules + Pathology
#
# Variants:
#   1. full expanded clinical variables, including medications
#   2. expanded clinical variables with all medication variables removed
#   3. expanded clinical variables with only statin and aspirin removed
target_omics_modalities <- c("lipids", "small_molecules", "pathology")

missing_target_modalities <- setdiff(target_omics_modalities, names(omics_list))

if (length(missing_target_modalities) > 0) {
  stop(
    "Target omics modalities are missing from omics_list: ",
    paste(missing_target_modalities, collapse = ", ")
  )
}

make_model_spec <- function(model_name, clinical_set, omics_modalities) {
  list(
    model_name = model_name,
    clinical = clinical_set,
    omics_modalities = omics_modalities,
    omics = unlist(omics_list[omics_modalities], use.names = FALSE)
  )
}

model_sets <- list(
  candidate_full_with_medications = make_model_spec(
    model_name = "candidate_full_with_medications",
    clinical_set = clinical_expanded_vars,
    omics_modalities = target_omics_modalities
  ),
  candidate_no_all_medications = make_model_spec(
    model_name = "candidate_no_all_medications",
    clinical_set = clinical_expanded_no_all_medications_vars,
    omics_modalities = target_omics_modalities
  ),
  candidate_no_statin_aspirin = make_model_spec(
    model_name = "candidate_no_statin_aspirin",
    clinical_set = clinical_expanded_no_statin_aspirin_vars,
    omics_modalities = target_omics_modalities
  )
)

model_run_mode <- "focused_medication_sensitivity"

# Hyperparameter grids
#
# alpha controls the elastic-net penalty:
#   alpha = 1.00 approximates LASSO
#   alpha = 0.50 balances LASSO and ridge behaviour
#   alpha = 0.25 is closer to ridge but still permits feature selection
#
# alpha = 0 is excluded because pure ridge does not perform feature selection
# and previously produced difficult-to-interpret coefficient patterns.
alpha_grid <- ccrvi_config$alpha_grid

# top_k controls the number of omics features retained after training-fold-only
# univariate AUC screening. Clinical features are not subject to this AUC screen.
top_k_grid <- ccrvi_config$top_k_grid

# Helper function: class weights
#
# MACE is relatively uncommon in the complete cohort. Positive-class weighting
# partially compensates for imbalance during model fitting by increasing the
# contribution of MACE cases to the binomial loss function.
make_weights <- function(y) {
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)

  if (n_pos == 0 || n_neg == 0) {
    return(rep(1, length(y)))
  }

  ifelse(y == 1, n_neg / n_pos, 1)
}

# Helper function: prepare raw model dataframe
#
# This selects the outcome and requested predictors. Missingness is not handled
# at this stage because imputation must be learned within training data only.
prepare_model_data <- function(data, clinical_vars, omics_vars) {
  selected_vars <- unique(c(clinical_vars, omics_vars))

  data %>%
    select(mace, all_of(selected_vars)) %>%
    mutate(outcome = as.numeric(mace)) %>%
    select(-mace)
}

# Helper function: safe checkpoint file names
#
# Model names are converted to file-safe strings before being used in checkpoint
# paths. Each model has a separate checkpoint file, allowing the script to resume
# after interruption without rerunning completed models.
safe_checkpoint_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_]+", "_") %>%
    str_replace_all("_+", "_")
}


# Helper function: identify omics columns after preprocessing
#
# Categorical clinical variables become dummy variables during preprocessing.
# Omics variables remain identifiable because their original modality prefixes
# are preserved in the column names.
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

# Helper function: select top omics features by univariate AUC
#
# This function is called only on training data inside each nested fold.
# It therefore avoids information leakage from validation or test folds.
#
# AUC strength is defined as max(AUC, 1 - AUC), so features with inverse
# association are not discarded simply because their raw AUC is below 0.5.
select_auc_features <- function(x_train, y_train, top_k, prefixes = omics_prefixes) {
  omics_cols <- get_omics_features(x_train, prefixes)

  if (length(omics_cols) == 0 || is.na(top_k)) {
    return(character())
  }

  auc_tbl <- map_dfr(omics_cols, function(feature_name) {
    feature_values <- x_train[, feature_name]

    if (length(unique(feature_values[!is.na(feature_values)])) < 2) {
      return(
        tibble(
          feature = feature_name,
          auc_value = NA_real_,
          auc_strength = NA_real_
        )
      )
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
  })

  auc_tbl %>%
    filter(!is.na(auc_strength)) %>%
    arrange(desc(auc_strength)) %>%
    slice_head(n = min(top_k, nrow(.))) %>%
    pull(feature)
}

# Helper function: preprocessing
#
# All preprocessing parameters are learned from the training data and applied to
# the corresponding test/validation data.
#
# Steps:
#   1. Separate numeric and categorical predictors.
#   2. Impute numeric variables using training medians.
#   3. Encode categorical variables using dummy variables trained on training data.
#   4. Remove zero-variance predictors based on training data.
#   5. Centre and scale predictors using training means and standard deviations.
#
# This structure is deliberately conservative because preprocessing leakage can
# inflate validation performance in small biomedical cohorts.
prep_train_test <- function(train_df, test_df) {
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

    # Convert categorical predictors to dummy variables.
    #
    # caret::dummyVars can fail inside future workers for some factor/formula
    # environments. model.matrix is used here instead because it is base R,
    # deterministic, and easier to align between training and validation data.
    #
    # Levels are learned from the training data only. Any unseen validation
    # category is mapped to "Other" before dummy encoding.
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

    # Ensure test data has exactly the same dummy columns as training data.
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

  # Remove zero-variance variables using training data only.
  train_sds <- apply(x_train, 2, sd, na.rm = TRUE)
  keep_vars <- names(train_sds)[is.finite(train_sds) & train_sds > 0]

  x_train <- x_train[, keep_vars, drop = FALSE]
  x_test <- x_test[, keep_vars, drop = FALSE]

  # Scale using training parameters only.
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

# Helper function: fit weighted elastic-net logistic regression
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
    parallel = use_glmnet_parallel
  )

  final_fit <- glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = alpha_val,
    weights = w_train,
    lambda = cv_fit$lambda.min
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

# Helper function: calibration metrics
#
# Calibration describes whether predicted risks correspond to observed risk.
# Brier score is the mean squared error of probabilistic predictions.
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

# Helper function: threshold metrics
#
# Youden's J is used here as a descriptive threshold. It should not be treated
# as a clinically final decision threshold without further validation and a
# clinical utility analysis.
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

# Helper function: AUPRC
#
# AUPRC is useful for imbalanced outcomes because it focuses on positive-class
# retrieval. Here, MACE cases are treated as the positive class.
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

# Create global validation folds
#
# These folds are reused across all model comparisons so that each model is
# evaluated on identical held-out patients. This makes fold-level comparisons
# more interpretable.
global_model_data <- model_source %>%
  select(mace) %>%
  mutate(outcome = as.numeric(mace))

validation_folds_global <- caret::createFolds(
  factor(global_model_data$outcome),
  k = 5,
  returnTrain = FALSE
)

# Main model-running function
#
# For each validation fold:
#   1. Hold out one-fifth of the cohort for validation.
#   2. Use the remaining development data for nested tuning.
#   3. Select alpha and top_k based on development-only nested performance.
#   4. Refit on the full development set using selected hyperparameters.
#   5. Predict the held-out validation fold.
run_aucfilter_model <- function(model_name, model_spec) {
  message("Running model: ", model_name)

  model_data <- prepare_model_data(
    data = model_source,
    clinical_vars = model_spec$clinical,
    omics_vars = model_spec$omics
  )

  has_omics <- length(model_spec$omics) > 0

  map_dfr(seq_along(validation_folds_global), function(validation_i) {
    message("  [", model_name, "] Validation fold ", validation_i)

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

      tuning_grid <- if (has_omics) {
        crossing(
          top_k = top_k_grid,
          alpha = alpha_grid
        )
      } else {
        crossing(
          top_k = NA_real_,
          alpha = alpha_grid
        )
      }

      tuning_results <- tuning_grid %>%
        mutate(
          auc = pmap_dbl(
            list(top_k, alpha),
            function(top_k, alpha) {
              clinical_features <- get_clinical_features(processed_outer$x_train)

              omics_features <- if (has_omics) {
                select_auc_features(
                  x_train = processed_outer$x_train,
                  y_train = processed_outer$y_train,
                  top_k = top_k,
                  prefixes = omics_prefixes
                )
              } else {
                character()
              }

              selected_features <- unique(c(clinical_features, omics_features))

              if (length(selected_features) < 2) {
                return(NA_real_)
              }

              x_train_sel <- processed_outer$x_train[, selected_features, drop = FALSE]
              x_test_sel <- processed_outer$x_test[, selected_features, drop = FALSE]

              fit_out <- tryCatch(
                fit_predict_enet(
                  x_train = x_train_sel,
                  y_train = processed_outer$y_train,
                  x_test = x_test_sel,
                  alpha_val = alpha
                ),
                error = function(e) {
                  message(
                    "    Fit failed in ", model_name,
                    ", validation fold ", validation_i,
                    ", nested fold ", outer_i,
                    ", top_k=", top_k,
                    ", alpha=", alpha,
                    ": ", conditionMessage(e)
                  )
                  NULL
                }
              )

              if (is.null(fit_out)) {
                return(NA_real_)
              }

              as.numeric(
                pROC::auc(
                  pROC::roc(
                    processed_outer$y_test,
                    fit_out$pred,
                    quiet = TRUE
                  )
                )
              )
            }
          )
        )

      best <- tuning_results %>%
        filter(!is.na(auc)) %>%
        arrange(desc(auc)) %>%
        slice(1)

      if (nrow(best) == 0) {
        return(
          tibble(
            validation_fold = validation_i,
            nested_outer_fold = outer_i,
            model = model_name,
            best_top_k = NA_real_,
            best_alpha = NA_real_,
            nested_auc = NA_real_
          )
        )
      }

      tibble(
        validation_fold = validation_i,
        nested_outer_fold = outer_i,
        model = model_name,
        best_top_k = best$top_k,
        best_alpha = best$alpha,
        nested_auc = best$auc
      )
    })

    # Select the final paired hyperparameter combination.
    # The combination with the highest mean nested AUC is selected. If tied, the
    # lower nested-fold variability is preferred, then the more frequently
    # selected combination.
    final_params <- nested_results %>%
      filter(!is.na(nested_auc), !is.na(best_alpha)) %>%
      group_by(best_top_k, best_alpha) %>%
      summarise(
        mean_nested_auc = mean(nested_auc, na.rm = TRUE),
        sd_nested_auc = sd(nested_auc, na.rm = TRUE),
        n_selected = n(),
        .groups = "drop"
      ) %>%
      arrange(
        desc(mean_nested_auc),
        sd_nested_auc,
        desc(n_selected)
      ) %>%
      slice(1)

    if (nrow(final_params) == 0) {
      stop(
        paste0(
          "No valid hyperparameter combination for model '", model_name,
          "' in validation fold ", validation_i,
          ". Nested results rows: ", nrow(nested_results),
          ". Valid nested AUC rows: ", sum(!is.na(nested_results$nested_auc))
        )
      )
    }

    final_top_k <- final_params$best_top_k
    final_alpha <- final_params$best_alpha

    processed_final <- prep_train_test(
      train_df = development_df,
      test_df = validation_df
    )

    clinical_features <- get_clinical_features(processed_final$x_train)

    omics_features <- if (has_omics) {
      select_auc_features(
        x_train = processed_final$x_train,
        y_train = processed_final$y_train,
        top_k = final_top_k,
        prefixes = omics_prefixes
      )
    } else {
      character()
    }

    selected_final <- unique(c(clinical_features, omics_features))

    if (length(selected_final) < 2) {
      stop(
        paste0(
          "Fewer than two predictors selected for model '", model_name,
          "' in validation fold ", validation_i
        )
      )
    }

    x_train_final <- processed_final$x_train[, selected_final, drop = FALSE]
    x_val_final <- processed_final$x_test[, selected_final, drop = FALSE]

    final_fit <- tryCatch(
      fit_predict_enet(
        x_train = x_train_final,
        y_train = processed_final$y_train,
        x_test = x_val_final,
        alpha_val = final_alpha
      ),
      error = function(e) {
        message(
          "  Final model failed in ", model_name,
          ", validation fold ", validation_i,
          ": ", conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(final_fit)) {
      stop(
        paste0(
          "Final model fit failed for model '", model_name,
          "' in validation fold ", validation_i
        )
      )
    }

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
      model_family = case_when(
        length(model_spec$clinical) == 0 ~ "omics_only",
        identical(model_spec$clinical, clinical_baseline_vars) ~ "clinical_baseline_plus_or_only",
        identical(model_spec$clinical, clinical_expanded_vars) ~ "clinical_expanded_plus_or_only",
        TRUE ~ "other"
      ),
      omics_modalities = paste(model_spec$omics_modalities, collapse = "+"),
      n_omics_modalities = length(model_spec$omics_modalities),
      final_top_k = final_top_k,
      final_alpha = final_alpha,
      final_lambda_min = final_fit$lambda_min,
      final_lambda_1se = final_fit$lambda_1se,
      n_features_available = ncol(processed_final$x_train),
      n_clinical_features = length(clinical_features),
      n_omics_prefilter_features = length(omics_features),
      n_total_prefilter_features = length(selected_final),
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

# Run all models with checkpointing
#
# Each model specification is saved as soon as it finishes. This is critical for
# long exhaustive runs: if a later aggregation step fails, completed model files
# remain recoverable in outputs/models/checkpoints.
checkpoint_dir <- ccrvi_checkpoint_directory("outputs/models/checkpoints_medication_sensitivity_script03")
dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)

validate_checkpoint_result <- function(result, model_name) {
  required_columns <- c(
    "model",
    "validation_fold",
    "predictions",
    "selected_features",
    "nested_results"
  )

  if (!inherits(result, "data.frame")) {
    return(FALSE)
  }

  if ("model_run_error" %in% names(result)) {
    return(FALSE)
  }

  if (!all(required_columns %in% names(result))) {
    return(FALSE)
  }

  if (nrow(result) != length(validation_folds_global)) {
    return(FALSE)
  }

  if (!setequal(result$validation_fold, seq_along(validation_folds_global))) {
    return(FALSE)
  }

  if (!all(result$model == model_name)) {
    return(FALSE)
  }

  TRUE
}

run_one_model_checkpointed <- function(model_name, model_index = NA_integer_, total_models = NA_integer_) {
  checkpoint_file <- file.path(
    checkpoint_dir,
    paste0(safe_checkpoint_name(model_name), ".rds")
  )

  progress_label <- if (!is.na(model_index) && !is.na(total_models)) {
    paste0("[", model_index, "/", total_models, "] ")
  } else {
    ""
  }

  if (file.exists(checkpoint_file)) {
    cached_result <- readRDS(checkpoint_file)

    if (validate_checkpoint_result(cached_result, model_name)) {
      message(progress_label, "Loading valid checkpoint: ", model_name)
      return(cached_result)
    }

    message(progress_label, "Ignoring invalid checkpoint and rerunning: ", model_name)
  }

  message(progress_label, "Running checkpointed model: ", model_name)
  start_time <- Sys.time()

  result <- tryCatch(
    run_aucfilter_model(
      model_name = model_name,
      model_spec = model_sets[[model_name]]
    ),
    error = function(e) {
      tibble(
        model = model_name,
        model_run_error = conditionMessage(e)
      )
    }
  )

  elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  if (!validate_checkpoint_result(result, model_name)) {
    if (!"model_run_error" %in% names(result)) {
      result <- tibble(
        model = model_name,
        model_run_error = paste0(
          "Invalid model result structure. Rows: ", nrow(result),
          "; columns: ", paste(names(result), collapse = ", ")
        )
      )
    }
  }

  result <- result %>%
    mutate(model_elapsed_minutes = elapsed_minutes)

  saveRDS(result, checkpoint_file)
  message(progress_label, "Saved checkpoint: ", model_name, " (", round(elapsed_minutes, 2), " min)")
  result
}

model_names <- names(model_sets)
total_models <- length(model_names)
message("Prepared ", total_models, " model specifications for run mode: ", model_run_mode)

if (use_model_parallel) {
  model_result_list <- furrr::future_imap(
    model_names,
    ~ run_one_model_checkpointed(
      model_name = .x,
      model_index = .y,
      total_models = total_models
    ),
    .options = furrr::furrr_options(seed = TRUE)
  )
} else {
  model_result_list <- imap(
    model_names,
    ~ run_one_model_checkpointed(
      model_name = .x,
      model_index = .y,
      total_models = total_models
    )
  )
}

names(model_result_list) <- model_names

# Save the raw list as an additional recovery object before any flattening.
saveRDS(
  model_result_list,
  "outputs/models/medication_sensitivity_script03/03F_model_result_list_checkpointed.rds"
)

aucfilter_modality_results <- bind_rows(model_result_list)

# Separate successful and failed model outputs before post-processing.
if ("model_run_error" %in% names(aucfilter_modality_results)) {
  failed_model_runs <- aucfilter_modality_results %>%
    filter(!is.na(model_run_error))

  successful_model_results <- aucfilter_modality_results %>%
    filter(is.na(model_run_error)) %>%
    select(-model_run_error)
} else {
  failed_model_runs <- tibble()
  successful_model_results <- aucfilter_modality_results
}

if (nrow(failed_model_runs) > 0) {
  failed_model_runs %>%
    select(where(~ !is.list(.x))) %>%
    write.csv(
      "outputs/tables/medication_sensitivity_script03/03F_failed_model_runs.csv",
      row.names = FALSE
    )

  saveRDS(
    failed_model_runs,
    "outputs/models/medication_sensitivity_script03/03F_failed_model_runs.rds"
  )
}

if (!"predictions" %in% names(successful_model_results)) {
  stop(
    "No successful model results contained a predictions column. ",
    "Check outputs/tables/aucfilter_failed_model_runs.csv and checkpoint files."
  )
}

# Extract nested list-columns into flat tables
modality_predictions <- successful_model_results %>%
  select(predictions) %>%
  tidyr::unnest(cols = c(predictions))

modality_selected_features <- successful_model_results %>%
  select(selected_features) %>%
  tidyr::unnest(cols = c(selected_features))

modality_nested_results <- successful_model_results %>%
  select(nested_results) %>%
  tidyr::unnest(cols = c(nested_results))

modality_fold_performance <- successful_model_results %>%
  select(
    model,
    validation_fold,
    model_family,
    omics_modalities,
    n_omics_modalities,
    final_top_k,
    final_alpha,
    final_lambda_min,
    final_lambda_1se,
    n_features_available,
    n_clinical_features,
    n_omics_prefilter_features,
    n_total_prefilter_features,
    n_model_features_selected,
    validation_auc,
    validation_auprc,
    model_elapsed_minutes
  )

# Summarise validation-fold performance
modality_performance_summary <- modality_fold_performance %>%
  group_by(model, model_family, omics_modalities, n_omics_modalities) %>%
  summarise(
    mean_auc = mean(validation_auc, na.rm = TRUE),
    sd_auc = sd(validation_auc, na.rm = TRUE),
    min_auc = min(validation_auc, na.rm = TRUE),
    max_auc = max(validation_auc, na.rm = TRUE),
    mean_auprc = mean(validation_auprc, na.rm = TRUE),
    sd_auprc = sd(validation_auprc, na.rm = TRUE),
    mean_omics_prefilter_features = mean(n_omics_prefilter_features, na.rm = TRUE),
    mean_total_prefilter_features = mean(n_total_prefilter_features, na.rm = TRUE),
    mean_model_features_selected = mean(n_model_features_selected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_auc))

# Pooled out-of-fold performance
#
# Predictions from the five validation folds are pooled so that each participant
# contributes exactly one out-of-fold prediction. This is useful for estimating
# overall discrimination, AUPRC, calibration, and threshold-based metrics.
pooled_modality_performance <- modality_predictions %>%
  group_by(model) %>%
  group_modify(~ {
    roc_obj <- pROC::roc(.x$outcome, .x$pred, quiet = TRUE)
    roc_ci <- pROC::ci.auc(roc_obj)

    pr_auc <- calc_auprc(
      y_true = .x$outcome,
      pred = .x$pred
    )

    cal <- calc_calibration(
      y_true = .x$outcome,
      pred = .x$pred
    )

    thresh <- calc_threshold_metrics(
      y_true = .x$outcome,
      pred = .x$pred
    )

    tibble(
      pooled_auc = as.numeric(pROC::auc(roc_obj)),
      auc_ci_low = as.numeric(roc_ci[1]),
      auc_ci_high = as.numeric(roc_ci[3]),
      pooled_auprc = pr_auc
    ) %>%
      bind_cols(cal) %>%
      bind_cols(thresh)
  }) %>%
  ungroup() %>%
  left_join(
    modality_fold_performance %>%
      distinct(model, model_family, omics_modalities, n_omics_modalities),
    by = "model"
  ) %>%
  arrange(desc(pooled_auc))

full_candidate_reference <- pooled_modality_performance %>%
  filter(model == "candidate_full_with_medications") %>%
  slice(1)

pooled_modality_performance <- pooled_modality_performance %>%
  mutate(
    delta_auc_vs_full_candidate = pooled_auc - full_candidate_reference$pooled_auc,
    delta_auprc_vs_full_candidate = pooled_auprc - full_candidate_reference$pooled_auprc,
    delta_brier_vs_full_candidate = brier_score - full_candidate_reference$brier_score
  ) %>%
  arrange(desc(pooled_auc))

# Paired DeLong comparisons against the full medication-including candidate.
#
# These tests use matched pooled out-of-fold predictions from the same validation
# folds. They assess whether medication-restricted variants differ in ROC AUC
# from the full candidate model.
safe_delong <- function(prediction_df, model_name, benchmark_name) {
  pair_df <- prediction_df %>%
    filter(model %in% c(model_name, benchmark_name)) %>%
    select(model, outcome, pred) %>%
    group_by(model) %>%
    mutate(row_id = row_number()) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = row_id,
      names_from = model,
      values_from = pred
    )

  outcome_df <- prediction_df %>%
    filter(model == model_name) %>%
    select(outcome) %>%
    mutate(row_id = row_number())

  pair_df <- pair_df %>%
    left_join(outcome_df, by = "row_id") %>%
    filter(
      !is.na(outcome),
      !is.na(.data[[model_name]]),
      !is.na(.data[[benchmark_name]])
    )

  if (nrow(pair_df) < 10 || n_distinct(pair_df$outcome) != 2) {
    return(tibble(
      model = model_name,
      benchmark_model = benchmark_name,
      auc_model = NA_real_,
      auc_benchmark = NA_real_,
      delta_auc = NA_real_,
      delong_p_value = NA_real_
    ))
  }

  roc_model <- pROC::roc(
    pair_df$outcome,
    pair_df[[model_name]],
    quiet = TRUE
  )

  roc_benchmark <- pROC::roc(
    pair_df$outcome,
    pair_df[[benchmark_name]],
    quiet = TRUE
  )

  test_out <- tryCatch(
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
    delong_p_value = ifelse(is.null(test_out), NA_real_, test_out$p.value)
  )
}

medication_sensitivity_delong <- setdiff(
  unique(modality_predictions$model),
  "candidate_full_with_medications"
) %>%
  map_dfr(~ safe_delong(
    prediction_df = modality_predictions,
    model_name = .x,
    benchmark_name = "candidate_full_with_medications"
  )) %>%
  arrange(desc(delta_auc))

# Feature-selection stability
#
# A feature selected in many validation folds is more stable than one selected
# once. This does not prove biological importance, but it helps distinguish
# reproducible model signals from fold-specific noise.
modality_feature_votes <- modality_selected_features %>%
  group_by(model, feature) %>%
  summarise(
    selected_n = n(),
    selection_frequency = selected_n / 5,
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
  left_join(
    modality_fold_performance %>%
      distinct(model, model_family, omics_modalities, n_omics_modalities),
    by = "model"
  ) %>%
  arrange(model, desc(selection_frequency), desc(mean_abs_coefficient))

# Compare baseline clinical and expanded clinical versions directly
#
# This table helps determine whether adding broader clinical/liver-disease
# context improves performance relative to the compact clinical baseline.
clinical_family_comparison <- pooled_modality_performance %>%
  filter(model_family %in% c(
    "clinical_baseline_plus_or_only",
    "clinical_expanded_plus_or_only"
  )) %>%
  select(
    model_family,
    omics_modalities,
    n_omics_modalities,
    model,
    pooled_auc,
    auc_ci_low,
    auc_ci_high,
    pooled_auprc,
    calibration_intercept,
    calibration_slope,
    brier_score
  ) %>%
  arrange(omics_modalities, model_family)

# Save all outputs
aucfilter_modality_outputs <- list(
  model_run_mode = model_run_mode,
  model_sets = model_sets,
  omics_definitions = omics_definitions,
  clinical_baseline_vars = clinical_baseline_vars,
  clinical_expanded_vars = clinical_expanded_vars,
  variable_count_table = variable_count_table,
  aucfilter_modality_results = aucfilter_modality_results,
  successful_model_results = successful_model_results,
  failed_model_runs = failed_model_runs,
  model_result_list = model_result_list,
  modality_predictions = modality_predictions,
  modality_selected_features = modality_selected_features,
  modality_nested_results = modality_nested_results,
  modality_fold_performance = modality_fold_performance,
  modality_performance_summary = modality_performance_summary,
  pooled_modality_performance = pooled_modality_performance,
  modality_feature_votes = modality_feature_votes,
  medication_sensitivity_delong = medication_sensitivity_delong,
  clinical_family_comparison = clinical_family_comparison
)

saveRDS(
  aucfilter_modality_outputs,
  "outputs/models/medication_sensitivity_script03/03F_medication_sensitivity_script03_outputs.rds"
)

write.csv(
  modality_performance_summary,
  "outputs/tables/medication_sensitivity_script03/03F_modality_performance_summary.csv",
  row.names = FALSE
)

write.csv(
  pooled_modality_performance,
  "outputs/tables/medication_sensitivity_script03/03F_pooled_performance.csv",
  row.names = FALSE
)

write.csv(
  modality_fold_performance,
  "outputs/tables/medication_sensitivity_script03/03F_fold_performance.csv",
  row.names = FALSE
)

write.csv(
  modality_feature_votes,
  "outputs/tables/medication_sensitivity_script03/03F_feature_votes.csv",
  row.names = FALSE
)

write.csv(
  modality_nested_results,
  "outputs/tables/medication_sensitivity_script03/03F_nested_results.csv",
  row.names = FALSE
)

write.csv(
  clinical_family_comparison,
  "outputs/tables/medication_sensitivity_script03/03F_clinical_comparison.csv",
  row.names = FALSE
)

write.csv(
  medication_sensitivity_delong,
  "outputs/tables/medication_sensitivity_script03/03F_delong_vs_full_candidate.csv",
  row.names = FALSE
)

# Generate medication sensitivity figures
#
# This section is deliberately placed after all model tables have been created.
# It uses the existing objects from the original Script 03-derived pipeline and
# does not alter model fitting, feature filtering, tuning, or validation folds.

save_med_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path("outputs/figures/medication_sensitivity_script03", filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

model_label_table <- tibble(
  model = c(
    "candidate_full_with_medications",
    "candidate_no_all_medications",
    "candidate_no_statin_aspirin"
  ),
  model_label = c(
    "Full candidate\nwith medications",
    "Candidate\nno medications",
    "Candidate\nno statin/aspirin"
  ),
  short_label = c(
    "Full",
    "No medications",
    "No statin/aspirin"
  )
)

feature_modality_label <- function(feature_name) {
  case_when(
    str_starts(feature_name, "lipid_") ~ "Lipids",
    str_starts(feature_name, "sm_") ~ "Small molecules",
    str_starts(feature_name, "path_") ~ "Pathology",
    str_starts(feature_name, "med_") ~ "Medication",
    TRUE ~ "Expanded clinical"
  )
}

# AUC comparison
auc_comparison_plot <- pooled_modality_performance %>%
  left_join(model_label_table, by = "model") %>%
  mutate(
    model_label = factor(model_label, levels = model_label[order(pooled_auc)])
  ) %>%
  ggplot(aes(x = model_label, y = pooled_auc)) +
  geom_col() +
  geom_errorbar(aes(ymin = auc_ci_low, ymax = auc_ci_high), width = 0.2) +
  coord_flip() +
  labs(
    title = "Candidate model performance after medication removal",
    subtitle = "Pooled out-of-fold AUC with 95% confidence intervals",
    x = NULL,
    y = "Pooled AUC"
  ) +
  theme_minimal(base_size = 12)

save_med_plot(
  auc_comparison_plot,
  "03F_auc_comparison_medication_sensitivity.png",
  width = 8,
  height = 5
)

# Delta AUC relative to full candidate
full_candidate_auc <- pooled_modality_performance %>%
  filter(model == "candidate_full_with_medications") %>%
  pull(pooled_auc)

delta_auc_plot <- pooled_modality_performance %>%
  filter(model != "candidate_full_with_medications") %>%
  left_join(model_label_table, by = "model") %>%
  mutate(
    delta_auc_vs_full = pooled_auc - full_candidate_auc,
    short_label = factor(short_label, levels = short_label[order(delta_auc_vs_full)])
  ) %>%
  ggplot(aes(x = short_label, y = delta_auc_vs_full)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Change in AUC after medication-variable removal",
    subtitle = "Negative values indicate lower AUC than the full candidate model",
    x = NULL,
    y = "Delta pooled AUC versus full candidate"
  ) +
  theme_minimal(base_size = 12)

save_med_plot(
  delta_auc_plot,
  "03F_delta_auc_vs_full_candidate.png",
  width = 7,
  height = 4
)

# Fold-level validation AUC
fold_auc_plot <- modality_fold_performance %>%
  left_join(model_label_table, by = "model") %>%
  mutate(model_label = factor(model_label, levels = model_label_table$model_label)) %>%
  ggplot(aes(x = model_label, y = validation_auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.75, size = 2) +
  labs(
    title = "Validation-fold AUC stability",
    subtitle = "Each point represents one outer validation fold",
    x = NULL,
    y = "Validation-fold AUC"
  ) +
  theme_minimal(base_size = 12)

save_med_plot(
  fold_auc_plot,
  "03F_fold_auc_stability.png",
  width = 8,
  height = 5
)

# Sensitivity and specificity at Youden threshold
sensitivity_specificity_plot <- pooled_modality_performance %>%
  left_join(model_label_table, by = "model") %>%
  select(model, short_label, sensitivity, specificity) %>%
  pivot_longer(
    cols = c(sensitivity, specificity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric, sensitivity = "Sensitivity", specificity = "Specificity"),
    short_label = factor(short_label, levels = model_label_table$short_label)
  ) %>%
  ggplot(aes(x = short_label, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75)) +
  coord_flip() +
  labs(
    title = "Sensitivity and specificity after medication removal",
    subtitle = "Calculated at the pooled Youden threshold",
    x = NULL,
    y = "Metric value",
    fill = "Metric"
  ) +
  theme_minimal(base_size = 12)

save_med_plot(
  sensitivity_specificity_plot,
  "03F_sensitivity_specificity_medication_sensitivity.png",
  width = 8,
  height = 5
)

# ROC curves
roc_curve_data <- modality_predictions %>%
  group_by(model) %>%
  group_modify(~ {
    roc_obj <- pROC::roc(.x$outcome, .x$pred, quiet = TRUE)
    tibble(
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      false_positive_rate = 1 - specificity
    )
  }) %>%
  ungroup() %>%
  left_join(
    pooled_modality_performance %>% select(model, pooled_auc),
    by = "model"
  ) %>%
  left_join(model_label_table, by = "model") %>%
  mutate(curve_label = paste0(short_label, " | AUC ", round(pooled_auc, 3)))

write.csv(
  roc_curve_data,
  "outputs/tables/medication_sensitivity_script03/03F_roc_curve_data.csv",
  row.names = FALSE
)

roc_curve_plot <- roc_curve_data %>%
  ggplot(aes(x = false_positive_rate, y = sensitivity, colour = curve_label)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "ROC curves for medication sensitivity models",
    x = "False positive rate",
    y = "Sensitivity",
    colour = "Model"
  ) +
  theme_minimal(base_size = 11)

save_med_plot(
  roc_curve_plot,
  "03F_roc_curves_medication_sensitivity.png",
  width = 8,
  height = 6
)

# Feature stability and modality labels
modality_feature_votes_with_labels <- modality_feature_votes %>%
  mutate(feature_modality = feature_modality_label(feature)) %>%
  left_join(model_label_table, by = "model")

write.csv(
  modality_feature_votes_with_labels,
  "outputs/tables/medication_sensitivity_script03/03F_feature_votes_with_modality_labels.csv",
  row.names = FALSE
)

top_feature_stability_plot <- modality_feature_votes_with_labels %>%
  group_by(model) %>%
  slice_max(order_by = selection_frequency, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    feature = factor(feature, levels = rev(unique(feature))),
    model_label = factor(model_label, levels = model_label_table$model_label)
  ) %>%
  ggplot(aes(x = feature, y = selection_frequency, fill = feature_modality)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ model_label, scales = "free_y") +
  labs(
    title = "Top selected features after medication removal",
    subtitle = "Selection frequency across outer validation folds",
    x = NULL,
    y = "Selection frequency",
    fill = "Feature type"
  ) +
  theme_minimal(base_size = 10)

save_med_plot(
  top_feature_stability_plot,
  "03F_top_feature_stability_medication_sensitivity.png",
  width = 12,
  height = 8
)

# Coefficient magnitude
top_coefficient_plot <- modality_feature_votes_with_labels %>%
  group_by(model) %>%
  slice_max(order_by = mean_abs_coefficient, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    feature = factor(feature, levels = rev(unique(feature))),
    model_label = factor(model_label, levels = model_label_table$model_label)
  ) %>%
  ggplot(aes(x = feature, y = mean_abs_coefficient, fill = direction)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ model_label, scales = "free_y") +
  labs(
    title = "Top features by mean coefficient magnitude",
    subtitle = "Elastic Net coefficients averaged across outer validation folds",
    x = NULL,
    y = "Mean absolute coefficient",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 10)

save_med_plot(
  top_coefficient_plot,
  "03F_top_coefficient_magnitude_medication_sensitivity.png",
  width = 12,
  height = 8
)

# Stable feature counts by modality
feature_modality_summary <- modality_feature_votes_with_labels %>%
  mutate(stable_feature = selection_frequency >= 0.6) %>%
  group_by(model, model_label, feature_modality) %>%
  summarise(
    n_selected_at_least_once = n(),
    n_stable_features = sum(stable_feature, na.rm = TRUE),
    mean_selection_frequency = mean(selection_frequency, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  feature_modality_summary,
  "outputs/tables/medication_sensitivity_script03/03F_feature_modality_summary.csv",
  row.names = FALSE
)

feature_modality_summary_plot <- feature_modality_summary %>%
  mutate(model_label = factor(model_label, levels = model_label_table$model_label)) %>%
  ggplot(aes(x = feature_modality, y = n_stable_features, fill = feature_modality)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ model_label) +
  labs(
    title = "Stable selected features by modality",
    subtitle = "Stable = selected in at least 3 of 5 validation folds",
    x = "Feature type",
    y = "Number of stable features"
  ) +
  theme_minimal(base_size = 11)

save_med_plot(
  feature_modality_summary_plot,
  "03F_stable_features_by_modality.png",
  width = 10,
  height = 5
)

# Confusion matrices
confusion_matrix_data <- modality_predictions %>%
  left_join(
    pooled_modality_performance %>% select(model, threshold),
    by = "model"
  ) %>%
  left_join(model_label_table, by = "model") %>%
  mutate(
    observed = if_else(outcome == 1, "MACE", "No MACE"),
    predicted = if_else(pred >= threshold, "MACE", "No MACE"),
    cell = case_when(
      observed == "No MACE" & predicted == "No MACE" ~ "TN",
      observed == "No MACE" & predicted == "MACE" ~ "FP",
      observed == "MACE" & predicted == "No MACE" ~ "FN",
      observed == "MACE" & predicted == "MACE" ~ "TP",
      TRUE ~ NA_character_
    )
  ) %>%
  count(model, model_label, short_label, observed, predicted, cell) %>%
  group_by(model, observed) %>%
  mutate(
    row_total = sum(n),
    row_percent = n / row_total
  ) %>%
  ungroup()

write.csv(
  confusion_matrix_data,
  "outputs/tables/medication_sensitivity_script03/03F_confusion_matrix_data.csv",
  row.names = FALSE
)

confusion_matrix_plot <- confusion_matrix_data %>%
  mutate(
    observed = factor(observed, levels = c("MACE", "No MACE")),
    predicted = factor(predicted, levels = c("No MACE", "MACE")),
    model_label = factor(model_label, levels = model_label_table$model_label),
    label = paste0(cell, "\n", n, "\n(", sprintf("%.1f", 100 * row_percent), "% row)"),
    cell_type = if_else(cell %in% c("TN", "TP"), "Correct", "Incorrect")
  ) %>%
  ggplot(aes(x = predicted, y = observed, fill = cell_type)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = label), fontface = "bold", size = 4.2, lineheight = 0.9) +
  coord_equal() +
  facet_wrap(~ model_label) +
  scale_fill_manual(values = c("Correct" = "#4CAF50", "Incorrect" = "#E64B35")) +
  labs(
    title = "Confusion matrices for medication sensitivity models",
    subtitle = "Percentages are row percentages within observed outcome class",
    x = "Predicted outcome",
    y = "Observed outcome",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

save_med_plot(
  confusion_matrix_plot,
  "03F_confusion_matrices_medication_sensitivity.png",
  width = 10,
  height = 4.5
)

message("Figures saved to outputs/figures/medication_sensitivity_script03.")

# Print key result tables to console for immediate inspection.
print(modality_performance_summary)
print(pooled_modality_performance)
print(clinical_family_comparison)

message("Medication sensitivity comparison complete.")
message("Tables saved to outputs/tables/medication_sensitivity_script03.")
message("Model objects saved to outputs/models/medication_sensitivity_script03.")
message("Per-model checkpoints saved to outputs/models/checkpoints_medication_sensitivity_script03.")
