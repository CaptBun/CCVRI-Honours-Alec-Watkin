# Repository stage: no-aspirin
# Source supplied as: 03_aucfilter_expanded_clinical_no_aspirin_comparison(1).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# Aspirin-ablation modality comparison using two-layer nested cross-validation
#
# Project:
#   Cardiovascular Risk Stratification in Liver Cirrhosis: Multiomics integration
#
# Purpose:
#   This script evaluates expanded-clinical models after removing aspirin only, using the same omics combinations and nested validation framework as
#   the original modality comparison. The remaining medication indicators are retained. The new models are then
#   compared against all successfully completed models from the original run.
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
dir.create("outputs/models/checkpoints_v8_no_aspirin", showWarnings = FALSE, recursive = TRUE)

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

# Define aspirin exclusion explicitly
#
# All expanded clinical predictors are retained except med_aspirin. This creates
# a targeted sensitivity analysis that addresses the clinician's concern while
# preserving metformin, statin, PPI, GLP-1/GIP, and no-medication indicators.
aspirin_var <- "med_aspirin"

clinical_expanded_no_aspirin_vars <- setdiff(
  clinical_expanded_vars,
  aspirin_var
)

removed_aspirin_vars <- intersect(
  aspirin_var,
  clinical_expanded_vars
)

missing_expected_aspirin_vars <- setdiff(
  aspirin_var,
  names(model_source)
)

message(
  "Expanded clinical aspirin ablation: removed ",
  length(removed_aspirin_vars),
  " aspirin variable; retained ",
  length(clinical_expanded_no_aspirin_vars),
  " clinical variables."
)

if (length(missing_expected_aspirin_vars) > 0) {
  stop(
    "The expected aspirin variable was not present in the compiled cohort: ",
    paste(missing_expected_aspirin_vars, collapse = ", ")
  )
}

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
    "clinical_expanded_no_aspirin",
    names(omics_list)
  ),
  n_variables = c(
    length(clinical_expanded_no_aspirin_vars),
    map_int(omics_list, length)
  )
)

write.csv(
  variable_count_table,
  "outputs/tables/model_comparison_no_aspirin_variable_counts.csv",
  row.names = FALSE
)

# Generate model sets automatically
#
# The previous version manually listed every model combination. That becomes
# error-prone once bile acids and pathology are added. Here, combinations are
# generated programmatically, but the default is intentionally not exhaustive.
#
# With six non-clinical modalities, this focused exhaustive run produces:
#   one aspirin-excluded clinical model + (2^6 - 1) aspirin-excluded clinical-plus-omics
#   models = 64 model specifications.
#
# Each specification is then evaluated with two-layer nested cross-validation.
# This can become computationally expensive and increases the chance of
# over-interpreting small, fold-dependent differences. The default "screened"
# mode therefore runs a defensible first-pass comparison:
#   1. expanded clinical without medications
#   2. aspirin-excluded expanded clinical + each single modality
#   3. aspirin-excluded expanded clinical + all non-clinical modalities
#
# Set model_run_mode <- "exhaustive" only if you want every possible modality
# combination. This is better used after the screened run identifies a small
# subset of plausible modality groups.
model_run_mode <- ccrvi_config$model_run_mode

if (model_run_mode == "screened") {
  omics_combinations <- c(
    as.list(names(omics_list)),
    list(names(omics_list))
  )
} else if (model_run_mode == "pairs") {
  omics_combinations <- c(
    as.list(names(omics_list)),
    combn(names(omics_list), 2, simplify = FALSE),
    list(names(omics_list))
  )
} else if (model_run_mode == "exhaustive") {
  max_omics_combination_size <- length(omics_list)

  omics_combinations <- map(
    seq_len(max_omics_combination_size),
    ~ combn(names(omics_list), .x, simplify = FALSE)
  ) %>%
    flatten()
} else {
  stop("model_run_mode must be one of: screened, pairs, exhaustive")
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
  clinical_expanded_no_aspirin = make_model_spec(
    model_name = "clinical_expanded_no_aspirin",
    clinical_set = clinical_expanded_no_aspirin_vars,
    omics_modalities = character()
  )
)

for (combo in omics_combinations) {
  combo_name <- paste(combo, collapse = "_")

  model_name <- paste0("clinical_expanded_no_aspirin_", combo_name)

  model_sets[[model_name]] <- make_model_spec(
    model_name = model_name,
    clinical_set = clinical_expanded_no_aspirin_vars,
    omics_modalities = combo
  )
}

# Optional model filtering
#
# Running every model can be computationally expensive. This regular expression
# can be used to restrict the run without editing the model generation logic.
#
# Examples:
#   model_name_filter <- "bile_acids"
#   model_name_filter <- "clinical_baseline|clinical_expanded"
#   model_name_filter <- NULL
model_name_filter <- NULL

if (!is.null(model_name_filter)) {
  model_sets <- model_sets[str_detect(names(model_sets), model_name_filter)]
}

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
        identical(model_spec$clinical, clinical_expanded_no_aspirin_vars) ~
          "clinical_expanded_no_aspirin_plus_or_only",
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
checkpoint_dir <- ccrvi_checkpoint_directory("outputs/models/checkpoints_v8_no_aspirin")
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
  "outputs/models/aucfilter_model_result_list_checkpointed_v7_no_aspirin.rds"
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
      "outputs/tables/aucfilter_no_aspirin_failed_model_runs.csv",
      row.names = FALSE
    )

  saveRDS(
    failed_model_runs,
    "outputs/models/aucfilter_no_aspirin_failed_model_runs.rds"
  )
}

if (!"predictions" %in% names(successful_model_results)) {
  stop(
    "No successful model results contained a predictions column. ",
    "Check outputs/tables/aucfilter_no_aspirin_failed_model_runs.csv and checkpoint files."
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

# Compare aspirin-excluded models against the complete original run
#
# The original run is not recomputed here. Its saved output object is loaded and
# combined with the new aspirin-excluded results. This preserves the original model
# estimates and makes the ablation script substantially faster.
original_outputs_path <-
  ccrvi_config$files$original_results

if (!file.exists(original_outputs_path)) {
  stop(
    "Original model output file was not found at: ",
    original_outputs_path,
    ". Run the original script first or update original_outputs_path."
  )
}

original_outputs <- readRDS(original_outputs_path)

required_original_objects <- c(
  "modality_performance_summary",
  "pooled_modality_performance",
  "modality_fold_performance",
  "modality_predictions"
)

missing_original_objects <- setdiff(
  required_original_objects,
  names(original_outputs)
)

if (length(missing_original_objects) > 0) {
  stop(
    "The original output object is missing: ",
    paste(missing_original_objects, collapse = ", ")
  )
}

# Add an analysis-source flag before combining tables.
original_performance_summary <- original_outputs$modality_performance_summary %>%
  mutate(analysis_source = "original_run")

no_aspirin_performance_summary <- modality_performance_summary %>%
  mutate(analysis_source = "aspirin_excluded_run")

combined_performance_summary <- bind_rows(
  original_performance_summary,
  no_aspirin_performance_summary
) %>%
  arrange(desc(mean_auc))

original_pooled_performance <- original_outputs$pooled_modality_performance %>%
  mutate(analysis_source = "original_run")

no_aspirin_pooled_performance <- pooled_modality_performance %>%
  mutate(analysis_source = "aspirin_excluded_run")

combined_pooled_performance <- bind_rows(
  original_pooled_performance,
  no_aspirin_pooled_performance
) %>%
  arrange(desc(pooled_auc))

# Construct matched comparisons between each expanded model and its corresponding
# aspirin-excluded version. Omics-only and baseline-clinical models remain in the
# all-model comparison but do not have a medication-ablation counterpart.
no_aspirin_matched <- no_aspirin_pooled_performance %>%
  mutate(
    matched_original_model = case_when(
      model == "clinical_expanded_no_aspirin" ~ "clinical_expanded",
      str_starts(model, "clinical_expanded_no_aspirin_") ~
        str_replace(model, "^clinical_expanded_no_aspirin_", "clinical_expanded_"),
      TRUE ~ NA_character_
    )
  )

matched_aspirin_exclusion_comparison <- no_aspirin_matched %>%
  select(
    no_aspirin_model = model,
    matched_original_model,
    omics_modalities,
    no_aspirin_pooled_auc = pooled_auc,
    no_aspirin_auc_ci_low = auc_ci_low,
    no_aspirin_auc_ci_high = auc_ci_high,
    no_aspirin_pooled_auprc = pooled_auprc,
    no_aspirin_calibration_intercept = calibration_intercept,
    no_aspirin_calibration_slope = calibration_slope,
    no_aspirin_brier_score = brier_score
  ) %>%
  left_join(
    original_pooled_performance %>%
      select(
        matched_original_model = model,
        original_pooled_auc = pooled_auc,
        original_auc_ci_low = auc_ci_low,
        original_auc_ci_high = auc_ci_high,
        original_pooled_auprc = pooled_auprc,
        original_calibration_intercept = calibration_intercept,
        original_calibration_slope = calibration_slope,
        original_brier_score = brier_score
      ),
    by = "matched_original_model"
  ) %>%
  mutate(
    delta_auc_no_aspirin_minus_original =
      no_aspirin_pooled_auc - original_pooled_auc,
    delta_auprc_no_aspirin_minus_original =
      no_aspirin_pooled_auprc - original_pooled_auprc,
    delta_brier_no_aspirin_minus_original =
      no_aspirin_brier_score - original_brier_score
  ) %>%
  arrange(desc(no_aspirin_pooled_auc))

# Paired outer-fold differences use identical held-out patients because the global
# folds and random seed are unchanged. These fold-level differences are descriptive;
# five folds provide too few independent units for reliable formal inference.
matched_fold_comparison <- modality_fold_performance %>%
  mutate(
    matched_original_model = case_when(
      model == "clinical_expanded_no_aspirin" ~ "clinical_expanded",
      str_starts(model, "clinical_expanded_no_aspirin_") ~
        str_replace(model, "^clinical_expanded_no_aspirin_", "clinical_expanded_"),
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    no_aspirin_model = model,
    matched_original_model,
    validation_fold,
    no_aspirin_validation_auc = validation_auc,
    no_aspirin_validation_auprc = validation_auprc
  ) %>%
  left_join(
    original_outputs$modality_fold_performance %>%
      select(
        matched_original_model = model,
        validation_fold,
        original_validation_auc = validation_auc,
        original_validation_auprc = validation_auprc
      ),
    by = c("matched_original_model", "validation_fold")
  ) %>%
  mutate(
    delta_validation_auc =
      no_aspirin_validation_auc - original_validation_auc,
    delta_validation_auprc =
      no_aspirin_validation_auprc - original_validation_auprc
  )


# Save all outputs
aucfilter_aspirin_exclusion_outputs <- list(
  model_run_mode = model_run_mode,
  original_outputs_path = original_outputs_path,
  model_sets = model_sets,
  omics_definitions = omics_definitions,
  aspirin_var = aspirin_var,
  removed_aspirin_vars = removed_aspirin_vars,
  missing_expected_aspirin_vars = missing_expected_aspirin_vars,
  clinical_expanded_vars = clinical_expanded_vars,
  clinical_expanded_no_aspirin_vars = clinical_expanded_no_aspirin_vars,
  variable_count_table = variable_count_table,
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
  combined_performance_summary = combined_performance_summary,
  combined_pooled_performance = combined_pooled_performance,
  matched_aspirin_exclusion_comparison = matched_aspirin_exclusion_comparison,
  matched_fold_comparison = matched_fold_comparison
)

saveRDS(
  aucfilter_aspirin_exclusion_outputs,
  ccrvi_config$files$no_aspirin_results
)

write.csv(
  modality_performance_summary,
  "outputs/tables/aucfilter_no_aspirin_performance_summary.csv",
  row.names = FALSE
)

write.csv(
  pooled_modality_performance,
  "outputs/tables/aucfilter_no_aspirin_pooled_performance.csv",
  row.names = FALSE
)

write.csv(
  modality_fold_performance,
  "outputs/tables/aucfilter_no_aspirin_fold_performance.csv",
  row.names = FALSE
)

write.csv(
  modality_feature_votes,
  "outputs/tables/aucfilter_no_aspirin_feature_votes.csv",
  row.names = FALSE
)

write.csv(
  modality_nested_results,
  "outputs/tables/aucfilter_no_aspirin_nested_results.csv",
  row.names = FALSE
)

write.csv(
  combined_performance_summary,
  "outputs/tables/aucfilter_all_models_plus_no_aspirin_summary.csv",
  row.names = FALSE
)

write.csv(
  combined_pooled_performance,
  "outputs/tables/aucfilter_all_models_plus_no_aspirin_pooled.csv",
  row.names = FALSE
)

write.csv(
  matched_aspirin_exclusion_comparison,
  "outputs/tables/aucfilter_matched_aspirin_exclusion_comparison.csv",
  row.names = FALSE
)

write.csv(
  matched_fold_comparison,
  "outputs/tables/aucfilter_matched_aspirin_exclusion_fold_comparison.csv",
  row.names = FALSE
)

print(matched_aspirin_exclusion_comparison)
print(combined_pooled_performance)

message("Aspirin-exclusion model comparison complete.")
message("New model checkpoints: outputs/models/checkpoints_v8_no_aspirin")
message("Matched comparison: outputs/tables/aucfilter_matched_aspirin_exclusion_comparison.csv")
