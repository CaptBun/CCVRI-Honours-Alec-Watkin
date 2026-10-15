# Repository stage: eda
# Source supplied as: 02_exploratory_analysis_v2_pathology(5).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# Exploratory analysis
# Purpose:
# This script characterises the complete multimodal cohort used for downstream
# cardiovascular risk prediction modelling in cirrhosis.
#
# The script is designed to be descriptive rather than predictive. It summarises
# cohort structure, missingness, liver disease aetiology, MASLD stratification,
# metabolic syndrome burden, univariate clinical associations, univariate omics
# associations, within-modality correlation, cross-modality correlation, and PCA
# structure.

# Load packages
library(tidyverse)
library(janitor)
library(broom)
library(pROC)
library(pheatmap)
library(recipes)
library(factoextra)
library(naniar)
library(patchwork)

set.seed(123)

# Load data
# This file should be created by the compilation script after clinical and omics
# data have been cleaned, joined, restricted to the complete multimodal cohort,
# and filtered for excessive missingness if appropriate.
eda_input_file <- ifelse(
  file.exists(ccrvi_config$files$unfiltered_cohort),
  ccrvi_config$files$unfiltered_cohort,
  ccrvi_config$files$compiled_cohort
)

combined_data <- read.csv(eda_input_file)

message("EDA input file: ", eda_input_file)

# Create output folders
# Tables and figures are separated so that outputs can be easily located for
# thesis writing, presentation preparation, and later manuscript drafting.
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures/eda", showWarnings = FALSE, recursive = TRUE)

# Helper function for saving tables
# All EDA tables are written as CSV files so they can be inspected in R, Excel,
# or imported into a thesis/manuscript table workflow.
save_table <- function(x, filename) {
  write.csv(
    x,
    file = file.path("outputs/tables", filename),
    row.names = FALSE
  )
}

# Helper function for saving ggplot figures
# Figures are exported at 300 dpi to make them usable in presentations and
# thesis drafts without needing to regenerate them manually.
save_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path("outputs/figures/eda", filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

# Safe numeric conversion
# Some omics values may be imported as character values because of formatting,
# missing-value codes, or spreadsheet artefacts. This function converts values
# to numeric while suppressing warnings caused by expected missing or non-numeric
# entries.
safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

# Identify likely metadata variables by name
#
# The compiled dataset should already remove technical identifiers, but the EDA
# script applies this additional safeguard so that source IDs, visit labels,
# plate labels, sample IDs, and other non-biological columns are not accidentally
# analysed as molecular or pathology features.
is_metadata_variable <- function(variable_names) {
  str_detect(
    variable_names,
    regex(
      paste(
        c(
          "sourceID",
          "source.ID",
          "Study.ID.1",
          "sampleID",
          "sample.ID",
          "plateID",
          "plate.ID",
          "tubeLabel",
          "tube.Label",
          "Visit",
          "Timepoint",
          "Label",
          "Batch",
          "Location",
          "Sample$"
        ),
        collapse = "|"
      ),
      ignore_case = TRUE
    )
  )
}

# Mean absolute correlation from an upper-triangle correlation matrix
# This is used as a compact summary of redundancy within each modality.
mean_abs_cor <- function(cor_mat) {
  cor_mat[upper.tri(cor_mat)] %>%
    abs() %>%
    mean(na.rm = TRUE)
}

# Extract highly correlated feature pairs
# Highly correlated pairs are useful for understanding redundancy within omics
# modalities. This does not automatically remove features; it only documents
# correlation structure.
get_high_cor_pairs <- function(cor_mat, cutoff = 0.80) {
  
  if (ncol(cor_mat) < 2) {
    return(tibble())
  }
  
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  
  pair_df <- expand.grid(
    feature_1 = rownames(cor_mat),
    feature_2 = colnames(cor_mat),
    stringsAsFactors = FALSE
  )
  
  pair_df$correlation <- as.vector(cor_mat)
  
  pair_df %>%
    as_tibble() %>%
    filter(
      !is.na(correlation),
      abs(correlation) >= cutoff
    ) %>%
    arrange(desc(abs(correlation)))
}

# Prepare numeric data for omics correlation and PCA analyses
# This function performs conservative cleaning before unsupervised analysis:
#   1. keeps requested variables that are present in the data
#   2. converts them to numeric
#   3. removes variables with fewer than two observed values
#   4. removes variables with zero or undefined variance
#
# This avoids PCA and correlation errors caused by constant or empty columns.
clean_numeric_matrix <- function(data, vars) {
  
  mat <- data %>%
    select(any_of(vars)) %>%
    mutate(across(everything(), safe_numeric))
  
  mat <- mat %>%
    select(where(~ sum(!is.na(.x)) > 1))
  
  mat <- mat %>%
    select(where(~ {
      v <- var(.x, na.rm = TRUE)
      !is.na(v) && is.finite(v) && v > 0
    }))
  
  mat
}

# Create a within-modality correlation heatmap
# Spearman correlation is used because omics data are commonly non-normal and may
# contain outliers. Pairwise complete observations are used because missingness is
# expected in multimodal biomedical data.
make_correlation_heatmap <- function(mat, title, filename) {
  
  if (ncol(mat) < 2) {
    message("Skipping ", title, ": fewer than two valid variables.")
    return(NULL)
  }
  
  cor_mat <- cor(
    mat,
    use = "pairwise.complete.obs",
    method = "spearman"
  )
  
  keep <- apply(
    cor_mat,
    1,
    function(x) all(is.finite(x))
  )
  
  cor_mat <- cor_mat[keep, keep, drop = FALSE]
  
  if (ncol(cor_mat) < 2) {
    message("Skipping ", title, ": correlation matrix became too small after cleaning.")
    return(NULL)
  }
  
  png(
    filename = file.path("outputs/figures/eda", filename),
    width = 2200,
    height = 2000,
    res = 300
  )
  
  pheatmap(
    cor_mat,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = title
  )
  
  dev.off()
  
  cor_mat
}

# Continuous clinical group comparison
# This function compares one continuous variable between outcome groups using:
#   • Wilcoxon rank-sum test
#   • descriptive medians and IQRs
#   • univariate ROC AUC
#
# The AUC is descriptive only. It is not a validated model performance estimate.
run_continuous_group_test <- function(data, variable, group_var = "mace") {
  
  test_df <- data %>%
    select(
      group = all_of(group_var),
      value = all_of(variable)
    ) %>%
    mutate(
      group = factor(group),
      value = safe_numeric(value)
    ) %>%
    filter(
      !is.na(group),
      !is.na(value)
    )
  
  if (n_distinct(test_df$group) != 2 || nrow(test_df) < 5) {
    return(tibble())
  }
  
  group_summary <- test_df %>%
    group_by(group) %>%
    summarise(
      n = n(),
      median = median(value, na.rm = TRUE),
      iqr = IQR(value, na.rm = TRUE),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  wilcox_out <- tryCatch(
    wilcox.test(value ~ group, data = test_df),
    error = function(e) NULL
  )
  
  auc_out <- tryCatch(
    as.numeric(
      pROC::auc(
        pROC::roc(
          response = test_df$group,
          predictor = test_df$value,
          quiet = TRUE
        )
      )
    ),
    error = function(e) NA_real_
  )
  
  tibble(
    variable = variable,
    test = "wilcoxon_rank_sum",
    n_total = nrow(test_df),
    p_value = ifelse(is.null(wilcox_out), NA_real_, wilcox_out$p.value),
    auc = auc_out,
    auc_strength = ifelse(is.na(auc_out), NA_real_, max(auc_out, 1 - auc_out)),
    group_0_n = group_summary$n[group_summary$group == levels(test_df$group)[1]],
    group_1_n = group_summary$n[group_summary$group == levels(test_df$group)[2]],
    group_0_median = group_summary$median[group_summary$group == levels(test_df$group)[1]],
    group_1_median = group_summary$median[group_summary$group == levels(test_df$group)[2]],
    group_0_iqr = group_summary$iqr[group_summary$group == levels(test_df$group)[1]],
    group_1_iqr = group_summary$iqr[group_summary$group == levels(test_df$group)[2]]
  )
}

# Categorical clinical group comparison
# Fisher's exact test is used rather than chi-square because some clinical
# categories may be sparse in a relatively small complete multimodal cohort.
run_categorical_group_test <- function(data, variable, group_var = "mace") {
  
  test_df <- data %>%
    select(
      group = all_of(group_var),
      variable_value = all_of(variable)
    ) %>%
    mutate(
      group = factor(group),
      variable_value = factor(variable_value)
    ) %>%
    filter(
      !is.na(group),
      !is.na(variable_value)
    )
  
  if (
    n_distinct(test_df$group) != 2 ||
    n_distinct(test_df$variable_value) < 2 ||
    nrow(test_df) < 5
  ) {
    return(tibble())
  }
  
  tab <- table(test_df$variable_value, test_df$group)
  
  test_out <- tryCatch(
    fisher.test(tab),
    error = function(e) NULL
  )
  
  tibble(
    variable = variable,
    test = "fisher_exact",
    n_total = nrow(test_df),
    n_levels = n_distinct(test_df$variable_value),
    p_value = ifelse(is.null(test_out), NA_real_, test_out$p.value)
  )
}

# Detect sparse categorical tables before logistic regression
# Sparse cells can produce unstable odds ratios or complete/quasi-complete
# separation. Instead of forcing an unstable model, these variables are skipped
# and labelled in the output table.
has_sparse_or_separated_table <- function(data, variable, outcome = "mace", min_cell_n = 3) {
  
  tab_df <- data %>%
    select(
      outcome_value = all_of(outcome),
      predictor = all_of(variable)
    ) %>%
    filter(
      !is.na(outcome_value),
      !is.na(predictor)
    ) %>%
    mutate(
      outcome_value = factor(outcome_value),
      predictor = factor(predictor)
    )
  
  if (nrow(tab_df) < 10) {
    return(TRUE)
  }
  
  if (n_distinct(tab_df$outcome_value) != 2 || n_distinct(tab_df$predictor) < 2) {
    return(TRUE)
  }
  
  tab <- table(tab_df$predictor, tab_df$outcome_value)
  
  any(tab < min_cell_n)
}

# Univariate logistic regression
# This estimates crude associations between individual clinical variables and
# MACE. Odds ratios are descriptive and are not adjusted unless explicitly done
# in a separate model.
run_logistic_or <- function(data, variable, outcome = "mace") {
  
  model_df <- data %>%
    select(
      outcome_value = all_of(outcome),
      predictor = all_of(variable)
    ) %>%
    filter(
      !is.na(outcome_value),
      !is.na(predictor)
    )
  
  if (nrow(model_df) < 10 || n_distinct(model_df$outcome_value) != 2) {
    return(tibble())
  }
  
  if (is.numeric(model_df$predictor)) {
    if (var(model_df$predictor, na.rm = TRUE) == 0) {
      return(tibble())
    }
  } else {
    model_df <- model_df %>%
      mutate(predictor = factor(predictor))
    
    if (n_distinct(model_df$predictor) < 2) {
      return(tibble())
    }
  }
  
  if (!is.numeric(model_df$predictor)) {
    
    if (has_sparse_or_separated_table(data, variable, outcome, min_cell_n = 3)) {
      return(
        tibble(
          variable = variable,
          term = NA_character_,
          odds_ratio = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_,
          p.value = NA_real_,
          note = "Skipped: sparse categorical predictor or possible separation"
        )
      )
    }
  }
  
  fit <- tryCatch(
    glm(
      outcome_value ~ predictor,
      data = model_df,
      family = binomial()
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble())
  }
  
  broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      variable = variable,
      note = NA_character_
    ) %>%
    select(
      variable,
      term,
      odds_ratio = estimate,
      conf.low,
      conf.high,
      p.value,
      note
    )
}

# Descriptive omics univariate analysis
# Each molecular feature is analysed independently against MACE using:
#   • Wilcoxon rank-sum test
#   • univariate ROC AUC
#   • median difference between MACE and non-MACE groups
#
# AUC strength is defined as max(AUC, 1 - AUC). This treats variables with
# inverse associations as potentially informative while preserving the original
# AUC direction in the auc column.
#
# These results should be interpreted as exploratory biology, not final feature
# selection.
run_omics_univariate <- function(data, vars, modality_name, outcome = "mace") {
  
  map_dfr(vars, function(feature_name) {
    
    test_df <- data %>%
      select(
        outcome_value = all_of(outcome),
        feature_value = all_of(feature_name)
      ) %>%
      mutate(
        outcome_value = as.numeric(outcome_value),
        feature_value = safe_numeric(feature_value)
      ) %>%
      filter(
        !is.na(outcome_value),
        !is.na(feature_value)
      )
    
    if (
      nrow(test_df) < 10 ||
      n_distinct(test_df$outcome_value) != 2 ||
      var(test_df$feature_value, na.rm = TRUE) == 0
    ) {
      return(
        tibble(
          modality = modality_name,
          feature = feature_name,
          n = nrow(test_df),
          p_value = NA_real_,
          auc = NA_real_,
          auc_strength = NA_real_,
          median_no_mace = NA_real_,
          median_mace = NA_real_,
          median_difference = NA_real_
        )
      )
    }
    
    wilcox_out <- tryCatch(
      wilcox.test(feature_value ~ outcome_value, data = test_df),
      error = function(e) NULL
    )
    
    auc_out <- tryCatch(
      as.numeric(
        pROC::auc(
          pROC::roc(
            response = test_df$outcome_value,
            predictor = test_df$feature_value,
            quiet = TRUE
          )
        )
      ),
      error = function(e) NA_real_
    )
    
    median_no_mace <- median(
      test_df$feature_value[test_df$outcome_value == 0],
      na.rm = TRUE
    )
    
    median_mace <- median(
      test_df$feature_value[test_df$outcome_value == 1],
      na.rm = TRUE
    )
    
    tibble(
      modality = modality_name,
      feature = feature_name,
      n = nrow(test_df),
      p_value = ifelse(is.null(wilcox_out), NA_real_, wilcox_out$p.value),
      auc = auc_out,
      auc_strength = ifelse(is.na(auc_out), NA_real_, max(auc_out, 1 - auc_out)),
      median_no_mace = median_no_mace,
      median_mace = median_mace,
      median_difference = median_mace - median_no_mace
    )
  }) %>%
    group_by(modality) %>%
    mutate(
      p_fdr = p.adjust(p_value, method = "BH")
    ) %>%
    ungroup() %>%
    arrange(modality, desc(auc_strength))
}

# PCA by modality
# PCA is performed independently for each omics modality.
#
# The purpose is exploratory:
#   • visualise dominant unsupervised structure
#   • assess whether samples cluster by MACE outcome
#   • assess whether samples cluster by MASLD status
#   • estimate the number of PCs needed to explain modality variance
#
# KNN imputation is used only to allow PCA visualisation. These imputed values
# are not written back to the modelling dataset and are not used for final model
# training.
run_pca_modality <- function(data, vars, modality_name, outcome_var = "mace", strat_var = "masld") {
  
  message("Running PCA for ", modality_name)
  
  numeric_mat <- clean_numeric_matrix(data, vars)
  
  if (ncol(numeric_mat) < 2) {
    message("Skipping PCA for ", modality_name, ": fewer than two valid variables.")
    return(NULL)
  }
  
  row_missing <- rowMeans(is.na(numeric_mat))
  
  numeric_mat <- numeric_mat[row_missing < 0.50, , drop = FALSE]
  
  if (nrow(numeric_mat) < 10) {
    message("Skipping PCA for ", modality_name, ": fewer than ten samples after row missingness filtering.")
    return(NULL)
  }
  
  pca_meta <- data[row_missing < 0.50, ] %>%
    select(
      Study.ID,
      outcome_value = all_of(outcome_var),
      strat_value = all_of(strat_var)
    )
  
  pca_recipe <- recipe(~ ., data = numeric_mat) %>%
    step_impute_knn(all_predictors(), neighbors = 5) %>%
    prep(training = numeric_mat)
  
  imputed_mat <- bake(pca_recipe, new_data = numeric_mat)
  
  pca_fit <- prcomp(
    imputed_mat,
    center = TRUE,
    scale. = TRUE
  )
  
  variance_tbl <- tibble(
    modality = modality_name,
    pc = paste0("PC", seq_along(pca_fit$sdev)),
    variance_explained = pca_fit$sdev^2 / sum(pca_fit$sdev^2),
    cumulative_variance = cumsum(pca_fit$sdev^2 / sum(pca_fit$sdev^2))
  )
  
  save_table(
    variance_tbl,
    paste0("eda_pca_variance_", modality_name, ".csv")
  )
  
  scree_plot <- variance_tbl %>%
    slice_head(n = min(20, nrow(variance_tbl))) %>%
    ggplot(
      aes(
        x = factor(pc, levels = pc),
        y = variance_explained
      )
    ) +
    geom_col() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = paste0(modality_name, " PCA scree plot"),
      x = "Principal component",
      y = "Variance explained"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot(
    scree_plot,
    paste0("eda_pca_scree_", modality_name, ".png"),
    width = 8,
    height = 5
  )
  
  score_tbl <- as.data.frame(pca_fit$x[, 1:2, drop = FALSE]) %>%
    as_tibble() %>%
    bind_cols(pca_meta) %>%
    mutate(
      outcome_label = factor(
        outcome_value,
        levels = c(0, 1),
        labels = c("No MACE", "MACE")
      )
    )
  
  outcome_plot <- ggplot(
    score_tbl,
    aes(
      x = PC1,
      y = PC2,
      colour = outcome_label
    )
  ) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(level = 0.80, na.rm = TRUE) +
    labs(
      title = paste0(modality_name, " PCA by MACE status"),
      x = paste0(
        "PC1 (",
        round(100 * variance_tbl$variance_explained[1], 1),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(100 * variance_tbl$variance_explained[2], 1),
        "%)"
      ),
      colour = NULL
    ) +
    theme_minimal(base_size = 13)
  
  save_plot(
    outcome_plot,
    paste0("eda_pca_mace_", modality_name, ".png"),
    width = 7,
    height = 6
  )
  
  masld_plot <- ggplot(
    score_tbl,
    aes(
      x = PC1,
      y = PC2,
      colour = strat_value
    )
  ) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(level = 0.80, na.rm = TRUE) +
    labs(
      title = paste0(modality_name, " PCA by MASLD status"),
      x = paste0(
        "PC1 (",
        round(100 * variance_tbl$variance_explained[1], 1),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(100 * variance_tbl$variance_explained[2], 1),
        "%)"
      ),
      colour = "MASLD"
    ) +
    theme_minimal(base_size = 13)
  
  save_plot(
    masld_plot,
    paste0("eda_pca_masld_", modality_name, ".png"),
    width = 7,
    height = 6
  )
  
  list(
    pca_fit = pca_fit,
    variance_tbl = variance_tbl,
    scores = score_tbl
  )
}

# Standardise key outcome and stratification variables
# These variables are used repeatedly across analyses. Coercing them early avoids
# inconsistent behaviour between plotting, summary tables, and statistical tests.
combined_data <- combined_data %>%
  mutate(
    mace = as.numeric(mace),
    cvd = as.numeric(cvd),
    masld = factor(masld),
    primary_aetiology = factor(primary_aetiology)
  )

# Define clinical variable domains
# Clinical variables are grouped by biological or clinical role to make the EDA
# easier to read and to support later thesis tables.
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
  "mets_count",
  "family_history_cardiac"
)

liver_disease_vars <- c(
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

clinical_vars <- c(
  demographic_vars,
  lifestyle_vars,
  cardiometabolic_vars,
  liver_disease_vars,
  medication_vars
)

clinical_vars <- intersect(clinical_vars, names(combined_data))

clinical_variable_domains <- tibble(
  variable = clinical_vars,
  domain = case_when(
    variable %in% demographic_vars ~ "demographic",
    variable %in% lifestyle_vars ~ "lifestyle",
    variable %in% cardiometabolic_vars ~ "cardiometabolic",
    variable %in% liver_disease_vars ~ "liver_disease",
    variable %in% medication_vars ~ "medication",
    TRUE ~ "other"
  )
)

save_table(
  clinical_variable_domains,
  "eda_clinical_variable_domains.csv"
)

# Define omics modalities
# To add a new omics platform, add one row here with the modality name and the
# variable prefix created by the compilation script. All downstream EDA sections
# will then include the modality automatically.
omics_definitions <- tribble(
  ~modality,          ~prefix,
  "lipoproteins",     "lipo_",
  "lipids",           "lipid_",
  "cytokines",        "cyto_",
  "small_molecules",  "sm_",
  "bile_acids",       "ba_",
  "pathology",        "path_"
)

# Identify available variables for each omics modality
# Metadata columns such as timepoint and label are excluded because they are not
# biological measurements and should not be included in molecular EDA.
omics_list <- omics_definitions %>%
  mutate(
    variables = map(
      prefix,
      ~ names(combined_data) %>%
        str_subset(paste0("^", .x)) %>%
        discard(is_metadata_variable)
    )
  ) %>%
  select(modality, variables) %>%
  deframe()

omics_variable_counts <- imap_dfr(
  omics_list,
  ~ tibble(
    modality = .y,
    n_variables = length(.x)
  )
)

save_table(
  omics_variable_counts,
  "eda_omics_variable_counts.csv"
)

# Cohort overview
# This table provides a compact description of the modelling cohort, including
# outcome prevalence and selected clinical characteristics.
cohort_overview <- tibble(
  n_patients = nrow(combined_data),
  n_mace = sum(combined_data$mace == 1, na.rm = TRUE),
  mace_rate = mean(combined_data$mace == 1, na.rm = TRUE),
  n_cvd = sum(combined_data$cvd == 1, na.rm = TRUE),
  cvd_rate = mean(combined_data$cvd == 1, na.rm = TRUE),
  n_masld_yes = sum(combined_data$masld == "Yes", na.rm = TRUE),
  masld_rate = mean(combined_data$masld == "Yes", na.rm = TRUE),
  age_median = median(combined_data$age, na.rm = TRUE),
  age_iqr = IQR(combined_data$age, na.rm = TRUE),
  bmi_median = median(combined_data$bmi, na.rm = TRUE),
  bmi_iqr = IQR(combined_data$bmi, na.rm = TRUE),
  mets_count_median = median(combined_data$mets_count, na.rm = TRUE),
  mets_count_iqr = IQR(combined_data$mets_count, na.rm = TRUE)
)

save_table(cohort_overview, "eda_cohort_overview.csv")

# Outcome balance
# Outcome prevalence is important because MACE is expected to be a minority class.
# This affects interpretation of AUC, AUPRC, PPV, threshold selection, and model
# stability in later predictive analyses.
outcome_balance <- combined_data %>%
  count(mace, name = "n") %>%
  mutate(
    proportion = n / sum(n)
  )

save_table(outcome_balance, "eda_outcome_balance.csv")

masld_mace_balance <- combined_data %>%
  count(masld, mace, name = "n") %>%
  group_by(masld) %>%
  mutate(
    within_masld_proportion = n / sum(n)
  ) %>%
  ungroup()

save_table(masld_mace_balance, "eda_masld_mace_balance.csv")

masld_mace_plot <- ggplot(
  masld_mace_balance,
  aes(
    x = masld,
    y = within_masld_proportion,
    fill = factor(mace)
  )
) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "MACE prevalence by MASLD status",
    x = "MASLD status",
    y = "Within-group proportion",
    fill = "MACE"
  ) +
  theme_minimal(base_size = 13)

save_plot(
  masld_mace_plot,
  "eda_mace_by_masld.png",
  width = 6,
  height = 5
)

# Missingness analysis
# Missingness is assessed before modelling to describe data quality, identify
# platforms with substantial missingness, and document variables that may be
# unsuitable for downstream predictive modelling.
#
# This section is descriptive. Variable exclusion for modelling should be handled
# explicitly in the modelling scripts or compilation script.
all_omics_prefixes <- omics_definitions$prefix

missingness_table <- combined_data %>%
  summarise(across(everything(), ~ mean(is.na(.x)))) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "missing_prop"
  ) %>%
  mutate(
    modality = case_when(
      str_detect(variable, "^lipo_") ~ "lipoproteins",
      str_detect(variable, "^lipid_") ~ "lipids",
      str_detect(variable, "^cyto_") ~ "cytokines",
      str_detect(variable, "^sm_") ~ "small_molecules",
      str_detect(variable, "^ba_") ~ "bile_acids",
      str_detect(variable, "^path_") ~ "pathology",
      TRUE ~ "clinical_or_outcome"
    )
  ) %>%
  arrange(desc(missing_prop))

save_table(missingness_table, "eda_missingness_by_variable.csv")

missingness_modality_summary <- missingness_table %>%
  group_by(modality) %>%
  summarise(
    n_variables = n(),
    median_missing = median(missing_prop),
    mean_missing = mean(missing_prop),
    n_missing_gt_20 = sum(missing_prop > 0.20),
    n_missing_gt_50 = sum(missing_prop > 0.50),
    n_missing_gt_80 = sum(missing_prop > 0.80),
    .groups = "drop"
  )

save_table(
  missingness_modality_summary,
  "eda_missingness_by_modality.csv"
)

missingness_plot <- ggplot(
  missingness_table,
  aes(
    x = missing_prop,
    fill = modality
  )
) +
  geom_histogram(bins = 40, alpha = 0.8) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Feature missingness distribution",
    x = "Proportion missing",
    y = "Number of variables",
    fill = "Variable group"
  ) +
  theme_minimal(base_size = 13)

save_plot(
  missingness_plot,
  "eda_missingness_distribution.png",
  width = 8,
  height = 5
)

# Missingness visualisation for the first 100 variables
# This provides a quick visual check of whether missingness is scattered or
# block-structured. Block missingness may suggest modality-specific assay absence
# or data integration issues.
png(
  filename = file.path("outputs/figures/eda", "eda_missingness_map_first_100_variables.png"),
  width = 2200,
  height = 1600,
  res = 300
)

naniar::vis_miss(
  combined_data %>%
    select(1:min(100, ncol(combined_data)))
)

dev.off()

# Liver disease aetiology analysis
# This section describes how MACE outcomes are distributed across the primary
# liver disease aetiology groups. Both total-cohort and within-aetiology
# proportions are useful:
#   • total proportions show how much each group contributes to the complete cohort
#   • within-aetiology proportions show outcome prevalence within each group

aetiology_total <- combined_data %>%
  count(primary_aetiology, mace, name = "n") %>%
  mutate(
    proportion_total = n / nrow(combined_data)
  )

save_table(
  aetiology_total,
  "eda_aetiology_mace_total_proportions.csv"
)

aetiology_within <- combined_data %>%
  count(primary_aetiology, mace, name = "n") %>%
  group_by(primary_aetiology) %>%
  mutate(
    within_aetiology_proportion = n / sum(n)
  ) %>%
  ungroup()

save_table(
  aetiology_within,
  "eda_aetiology_mace_within_aetiology_proportions.csv"
)

aetiology_total_plot <- ggplot(
  aetiology_total,
  aes(
    x = reorder(primary_aetiology, -proportion_total),
    y = proportion_total,
    fill = factor(mace)
  )
) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Primary liver disease aetiology and MACE outcomes",
    subtitle = paste0("Proportion of complete cohort, n = ", nrow(combined_data)),
    x = "Primary liver disease aetiology",
    y = "Percentage of total cohort",
    fill = "MACE"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_plot(
  aetiology_total_plot,
  "eda_aetiology_mace_total_proportions.png",
  width = 9,
  height = 5
)

aetiology_within_plot <- ggplot(
  aetiology_within,
  aes(
    x = primary_aetiology,
    y = within_aetiology_proportion,
    fill = factor(mace)
  )
) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "MACE prevalence by primary liver disease aetiology",
    x = "Primary liver disease aetiology",
    y = "Within-aetiology proportion",
    fill = "MACE"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_plot(
  aetiology_within_plot,
  "eda_aetiology_mace_within_aetiology_proportions.png",
  width = 9,
  height = 5
)

masld_aetiology_summary <- combined_data %>%
  summarise(
    primary_masld_pct = mean(primary_aetiology == "MASLD", na.rm = TRUE) * 100,
    any_masld_pct = mean(masld == "Yes", na.rm = TRUE) * 100,
    masld_not_primary_pct = mean(
      masld == "Yes" & primary_aetiology != "MASLD",
      na.rm = TRUE
    ) * 100
  )

save_table(
  masld_aetiology_summary,
  "eda_primary_vs_any_masld_summary.csv"
)

# Metabolic syndrome analysis
# Metabolic syndrome burden is central to the MASLD-CVD relationship and is also
# a clinically interpretable risk construct. These summaries assess whether MACE
# and MASLD groups differ in metabolic syndrome component count.
mets_by_mace <- combined_data %>%
  group_by(mace) %>%
  summarise(
    n = n(),
    mean_mets = mean(mets_count, na.rm = TRUE),
    sd_mets = sd(mets_count, na.rm = TRUE),
    median_mets = median(mets_count, na.rm = TRUE),
    iqr_mets = IQR(mets_count, na.rm = TRUE),
    .groups = "drop"
  )

save_table(mets_by_mace, "eda_mets_count_by_mace.csv")

mets_by_masld <- combined_data %>%
  group_by(masld) %>%
  summarise(
    n = n(),
    mean_mets = mean(mets_count, na.rm = TRUE),
    sd_mets = sd(mets_count, na.rm = TRUE),
    median_mets = median(mets_count, na.rm = TRUE),
    iqr_mets = IQR(mets_count, na.rm = TRUE),
    .groups = "drop"
  )

save_table(mets_by_masld, "eda_mets_count_by_masld.csv")

mets_mace_plot <- ggplot(
  combined_data,
  aes(
    x = factor(mets_count),
    fill = factor(mace)
  )
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Metabolic syndrome component count by MACE status",
    x = "Metabolic syndrome component count",
    y = "Proportion",
    fill = "MACE"
  ) +
  theme_minimal(base_size = 13)

save_plot(
  mets_mace_plot,
  "eda_mets_count_by_mace.png",
  width = 7,
  height = 5
)

mets_masld_plot <- ggplot(
  combined_data,
  aes(
    x = factor(mets_count),
    fill = masld
  )
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Metabolic syndrome component count by MASLD status",
    x = "Metabolic syndrome component count",
    y = "Proportion",
    fill = "MASLD"
  ) +
  theme_minimal(base_size = 13)

save_plot(
  mets_masld_plot,
  "eda_mets_count_by_masld.png",
  width = 7,
  height = 5
)

# Metabolic syndrome logistic models
# The first model estimates the crude association between metabolic syndrome
# component count and MACE. The second includes MASLD status to assess whether
# the metabolic syndrome signal persists after accounting for MASLD.
mets_models <- bind_rows(
  run_logistic_or(combined_data, "mets_count", outcome = "mace") %>%
    mutate(model = "mace_unadjusted"),
  broom::tidy(
    glm(
      mace ~ masld + mets_count,
      data = combined_data,
      family = binomial()
    ),
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      variable = term,
      model = "mace_masld_adjusted"
    ) %>%
    select(
      variable,
      term,
      odds_ratio = estimate,
      conf.low,
      conf.high,
      p.value,
      model
    )
)

save_table(mets_models, "eda_mets_count_logistic_models.csv")

# Clinical univariate association analysis
# Clinical variables are summarised separately as continuous or categorical
# predictors. This supports interpretation of baseline clinical differences
# between MACE and non-MACE participants.
clinical_data <- combined_data %>%
  select(mace, all_of(clinical_vars))

clinical_numeric_vars <- clinical_data %>%
  select(-mace) %>%
  select(where(is.numeric)) %>%
  names()

clinical_categorical_vars <- setdiff(
  names(clinical_data %>% select(-mace)),
  clinical_numeric_vars
)

clinical_continuous_tests <- map_dfr(
  clinical_numeric_vars,
  ~ run_continuous_group_test(
    data = clinical_data,
    variable = .x,
    group_var = "mace"
  )
) %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_fdr)

save_table(
  clinical_continuous_tests,
  "eda_clinical_continuous_univariate_vs_mace.csv"
)

# Clinical univariate AUC table
#
# This table puts numeric clinical predictors on the same descriptive AUC scale
# used for omics features. For categorical clinical predictors, AUC is not
# calculated here because multi-level categorical coding can make a single
# univariate AUC difficult to interpret without specifying a reference contrast.
clinical_auc_summary <- clinical_continuous_tests %>%
  select(
    variable,
    n_total,
    auc,
    auc_strength,
    p_value,
    p_fdr,
    group_0_median,
    group_1_median,
    group_0_iqr,
    group_1_iqr
  ) %>%
  arrange(desc(auc_strength))

save_table(
  clinical_auc_summary,
  "eda_clinical_numeric_auc_summary_vs_mace.csv"
)

clinical_categorical_tests <- map_dfr(
  clinical_categorical_vars,
  ~ run_categorical_group_test(
    data = clinical_data,
    variable = .x,
    group_var = "mace"
  )
) %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_fdr)

save_table(
  clinical_categorical_tests,
  "eda_clinical_categorical_univariate_vs_mace.csv"
)

clinical_logistic_univariate <- map_dfr(
  clinical_vars,
  ~ run_logistic_or(
    data = combined_data,
    variable = .x,
    outcome = "mace"
  )
) %>%
  mutate(
    p_fdr = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(p_fdr)

save_table(
  clinical_logistic_univariate,
  "eda_clinical_univariate_logistic_or_vs_mace.csv"
)

# Omics univariate association analysis against MACE
# This runs the same descriptive univariate analysis for every modality defined
# in omics_definitions, including bile acids if ba_ variables are present.
omics_univariate_results <- imap_dfr(
  omics_list,
  ~ run_omics_univariate(
    data = combined_data,
    vars = .x,
    modality_name = .y,
    outcome = "mace"
  )
)

save_table(
  omics_univariate_results,
  "eda_omics_univariate_vs_mace.csv"
)

# Dedicated pathology univariate table
#
# Pathology variables are routine clinical laboratory measures rather than
# high-throughput molecular assays. They are included in the modality framework
# for consistency, but a dedicated export makes them easier to inspect and cite
# separately in the thesis.
pathology_univariate_results <- omics_univariate_results %>%
  filter(modality == "pathology") %>%
  arrange(desc(auc_strength))

save_table(
  pathology_univariate_results,
  "eda_pathology_univariate_vs_mace.csv"
)

# Plot the strongest univariate omics features by modality
# This figure is useful for exploratory interpretation but should not be treated
# as a feature-selection result for the final prediction model.
top_omics_plot <- omics_univariate_results %>%
  filter(!is.na(auc_strength)) %>%
  group_by(modality) %>%
  slice_max(auc_strength, n = 15) %>%
  ungroup() %>%
  mutate(
    feature = reorder(feature, auc_strength)
  ) %>%
  ggplot(
    aes(
      x = feature,
      y = auc_strength
    )
  ) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ modality, scales = "free_y") +
  labs(
    title = "Top univariate omics features by MACE AUC strength",
    subtitle = "AUC strength is max(AUC, 1 - AUC); descriptive only",
    x = NULL,
    y = "Univariate AUC strength"
  ) +
  theme_minimal(base_size = 12)

save_plot(
  top_omics_plot,
  "eda_top_omics_univariate_auc_strength.png",
  width = 11,
  height = 8
)

# MASLD-stratified omics univariate analysis
# This evaluates whether feature-outcome relationships appear different within
# MASLD strata. Small subgroup sizes mean these results should be interpreted as
# exploratory rather than confirmatory.
omics_univariate_masld_stratified <- combined_data %>%
  filter(!is.na(masld)) %>%
  group_split(masld) %>%
  map_dfr(function(df) {
    
    masld_label <- unique(df$masld)
    
    imap_dfr(
      omics_list,
      ~ run_omics_univariate(
        data = df,
        vars = .x,
        modality_name = .y,
        outcome = "mace"
      )
    ) %>%
      mutate(masld_group = as.character(masld_label))
  })

save_table(
  omics_univariate_masld_stratified,
  "eda_omics_univariate_vs_mace_by_masld.csv"
)

# Within-modality correlation analysis
# This section quantifies redundancy within each omics modality and creates
# heatmaps to visualise correlation structure. This is particularly relevant for
# elastic net modelling because highly correlated predictors can affect feature
# stability across validation folds.
modality_matrices <- imap(
  omics_list,
  ~ clean_numeric_matrix(combined_data, .x)
)

# Clinical model matrix for clinical correlation structure
# Categorical variables are expanded into dummy variables so clinical correlation
# can be assessed on a numeric design-matrix scale.
clinical_matrix <- combined_data %>%
  select(any_of(clinical_vars)) %>%
  mutate(
    across(where(is.character), factor),
    across(where(is.logical), factor)
  )

clinical_matrix <- clinical_matrix %>%
  select(where(~ n_distinct(.x, na.rm = TRUE) > 1))

clinical_model_matrix <- model.matrix(
  ~ .,
  data = clinical_matrix
)[, -1, drop = FALSE]

clinical_model_matrix <- clinical_model_matrix[
  ,
  apply(clinical_model_matrix, 2, var, na.rm = TRUE) > 0,
  drop = FALSE
]

cor_clinical <- make_correlation_heatmap(
  as.data.frame(clinical_model_matrix),
  "Expanded clinical variable Spearman correlation structure",
  "eda_correlation_heatmap_clinical.png"
)

cor_mats <- imap(
  modality_matrices,
  ~ make_correlation_heatmap(
    .x,
    paste0(.y, " Spearman correlation structure"),
    paste0("eda_correlation_heatmap_", .y, ".png")
  )
)

within_modality_correlation_summary <- imap_dfr(
  cor_mats,
  function(cor_mat, modality_name) {
    
    if (is.null(cor_mat)) {
      return(
        tibble(
          modality = modality_name,
          n_variables = NA_integer_,
          mean_abs_correlation = NA_real_,
          n_pairs_abs_cor_ge_0_80 = NA_integer_,
          n_pairs_abs_cor_ge_0_90 = NA_integer_
        )
      )
    }
    
    upper_vals <- cor_mat[upper.tri(cor_mat)]
    
    tibble(
      modality = modality_name,
      n_variables = ncol(cor_mat),
      mean_abs_correlation = mean(abs(upper_vals), na.rm = TRUE),
      n_pairs_abs_cor_ge_0_80 = sum(abs(upper_vals) >= 0.80, na.rm = TRUE),
      n_pairs_abs_cor_ge_0_90 = sum(abs(upper_vals) >= 0.90, na.rm = TRUE)
    )
  }
)

if (!is.null(cor_clinical)) {
  
  clinical_upper <- cor_clinical[upper.tri(cor_clinical)]
  
  within_modality_correlation_summary <- bind_rows(
    tibble(
      modality = "clinical",
      n_variables = ncol(cor_clinical),
      mean_abs_correlation = mean(abs(clinical_upper), na.rm = TRUE),
      n_pairs_abs_cor_ge_0_80 = sum(abs(clinical_upper) >= 0.80, na.rm = TRUE),
      n_pairs_abs_cor_ge_0_90 = sum(abs(clinical_upper) >= 0.90, na.rm = TRUE)
    ),
    within_modality_correlation_summary
  )
}

save_table(
  within_modality_correlation_summary,
  "eda_within_modality_correlation_summary.csv"
)

high_cor_pairs <- imap_dfr(
  cor_mats,
  function(cor_mat, modality_name) {
    
    if (is.null(cor_mat)) {
      return(tibble())
    }
    
    get_high_cor_pairs(cor_mat, cutoff = 0.80) %>%
      mutate(modality = modality_name)
  }
)

save_table(
  high_cor_pairs,
  "eda_highly_correlated_feature_pairs_abs_r_ge_0_80.csv"
)

# Cross-modality correlations
# To keep this output readable, only the top 100 strongest feature-feature
# correlations are exported for each pair of omics modalities.
cross_modality_correlations <- crossing(
  modality_1 = names(modality_matrices),
  modality_2 = names(modality_matrices)
) %>%
  filter(modality_1 < modality_2) %>%
  pmap_dfr(function(modality_1, modality_2) {
    
    mat_1 <- modality_matrices[[modality_1]]
    mat_2 <- modality_matrices[[modality_2]]
    
    if (ncol(mat_1) < 1 || ncol(mat_2) < 1) {
      return(tibble())
    }
    
    cor_mat <- cor(
      mat_1,
      mat_2,
      use = "pairwise.complete.obs",
      method = "spearman"
    )
    
    cor_df <- expand.grid(
      feature_1 = rownames(cor_mat),
      feature_2 = colnames(cor_mat),
      stringsAsFactors = FALSE
    )
    
    cor_df$correlation <- as.vector(cor_mat)
    
    cor_df %>%
      as_tibble() %>%
      filter(is.finite(correlation)) %>%
      mutate(
        modality_1 = modality_1,
        modality_2 = modality_2,
        abs_correlation = abs(correlation)
      ) %>%
      arrange(desc(abs_correlation)) %>%
      slice_head(n = 100)
  }) %>%
  arrange(desc(abs_correlation))

save_table(
  cross_modality_correlations,
  "eda_top_cross_modality_correlations.csv"
)

cross_modality_summary <- cross_modality_correlations %>%
  group_by(modality_1, modality_2) %>%
  summarise(
    top_abs_correlation = max(abs_correlation, na.rm = TRUE),
    median_top_100_abs_correlation = median(abs_correlation, na.rm = TRUE),
    n_top_pairs_abs_cor_ge_0_70 = sum(abs_correlation >= 0.70, na.rm = TRUE),
    n_top_pairs_abs_cor_ge_0_80 = sum(abs_correlation >= 0.80, na.rm = TRUE),
    .groups = "drop"
  )

save_table(
  cross_modality_summary,
  "eda_cross_modality_correlation_summary.csv"
)

# PCA analysis by modality
# PCA outputs include variance tables, scree plots, MACE-coloured score plots,
# and MASLD-coloured score plots for each available modality.
pca_outputs <- imap(
  omics_list,
  ~ run_pca_modality(
    data = combined_data,
    vars = .x,
    modality_name = .y,
    outcome_var = "mace",
    strat_var = "masld"
  )
)

# Cross-modality PC correlations
# Correlating leading PCs across modalities gives a compact view of whether
# dominant unsupervised structure is shared across platforms.
pc_score_tables <- imap(
  pca_outputs,
  function(pca_output, modality_name) {
    
    if (is.null(pca_output)) {
      return(NULL)
    }
    
    pca_output$scores %>%
      select(Study.ID, starts_with("PC")) %>%
      select(Study.ID, any_of(paste0("PC", 1:5))) %>%
      rename_with(
        ~ paste0(modality_name, "_", .x),
        starts_with("PC")
      )
  }
)

pc_score_tables <- pc_score_tables[!map_lgl(pc_score_tables, is.null)]

if (length(pc_score_tables) >= 2) {
  
  pc_combined <- reduce(
    pc_score_tables,
    inner_join,
    by = "Study.ID"
  )
  
  pc_cor <- pc_combined %>%
    select(-Study.ID) %>%
    cor(use = "pairwise.complete.obs", method = "spearman")
  
  png(
    filename = file.path("outputs/figures/eda", "eda_cross_modality_pc_correlation_heatmap.png"),
    width = 2200,
    height = 2000,
    res = 300
  )
  
  pheatmap(
    pc_cor,
    main = "Cross-modality principal component correlations"
  )
  
  dev.off()
  
  pc_cor_table <- expand.grid(
    pc_1 = rownames(pc_cor),
    pc_2 = colnames(pc_cor),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      correlation = as.vector(pc_cor)
    ) %>%
    filter(pc_1 < pc_2) %>%
    mutate(
      abs_correlation = abs(correlation)
    ) %>%
    arrange(desc(abs_correlation))
  
  save_table(
    pc_cor_table,
    "eda_cross_modality_pc_correlations.csv"
  )
}

# Modality-level summary table
# This table is intended as a thesis/manuscript-friendly overview of each data
# layer. It combines variable counts, missingness, univariate signal summaries,
# correlation redundancy, and PCA complexity.
pca_90_summary <- imap_dfr(
  pca_outputs,
  function(pca_output, modality_name) {
    
    if (is.null(pca_output)) {
      return(
        tibble(
          modality = modality_name,
          n_pcs_90_variance = NA_integer_,
          n_pcs_95_variance = NA_integer_
        )
      )
    }
    
    tibble(
      modality = modality_name,
      n_pcs_90_variance = which(pca_output$variance_tbl$cumulative_variance >= 0.90)[1],
      n_pcs_95_variance = which(pca_output$variance_tbl$cumulative_variance >= 0.95)[1]
    )
  }
)

omics_signal_summary <- omics_univariate_results %>%
  group_by(modality) %>%
  summarise(
    median_auc_strength = median(auc_strength, na.rm = TRUE),
    max_auc_strength = max(auc_strength, na.rm = TRUE),
    n_features_auc_strength_ge_0_60 = sum(auc_strength >= 0.60, na.rm = TRUE),
    n_features_fdr_lt_0_05 = sum(p_fdr < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

omics_missingness_summary <- missingness_table %>%
  filter(modality != "clinical_or_outcome") %>%
  group_by(modality) %>%
  summarise(
    n_features = n(),
    median_missing = median(missing_prop, na.rm = TRUE),
    mean_missing = mean(missing_prop, na.rm = TRUE),
    n_missing_gt_50 = sum(missing_prop > 0.50, na.rm = TRUE),
    .groups = "drop"
  )

modality_summary_table <- omics_missingness_summary %>%
  left_join(omics_signal_summary, by = "modality") %>%
  left_join(within_modality_correlation_summary, by = "modality") %>%
  left_join(pca_90_summary, by = "modality") %>%
  arrange(modality)

save_table(
  modality_summary_table,
  "eda_modality_summary_table.csv"
)

# Distribution diagnostics by modality
#
# These summaries describe the numerical shape of each data layer. Omics and
# pathology variables are often skewed, so documenting skewness, coefficient of
# variation, and median absolute deviation helps justify the use of Spearman
# correlation, non-parametric tests, and regularised models.
calc_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || sd(x) == 0) return(NA_real_)
  mean((x - mean(x))^3) / sd(x)^3
}

feature_distribution_summary <- imap_dfr(
  omics_list,
  function(vars, modality_name) {
    clean_numeric_matrix(combined_data, vars) %>%
      summarise(
        across(
          everything(),
          list(
            n_observed = ~ sum(!is.na(.x)),
            mean = ~ mean(.x, na.rm = TRUE),
            sd = ~ sd(.x, na.rm = TRUE),
            median = ~ median(.x, na.rm = TRUE),
            iqr = ~ IQR(.x, na.rm = TRUE),
            mad = ~ mad(.x, na.rm = TRUE),
            cv = ~ ifelse(abs(mean(.x, na.rm = TRUE)) > 0, sd(.x, na.rm = TRUE) / abs(mean(.x, na.rm = TRUE)), NA_real_),
            skewness = ~ calc_skewness(.x)
          ),
          .names = "{.col}__{.fn}"
        )
      ) %>%
      pivot_longer(
        everything(),
        names_to = c("feature", ".value"),
        names_sep = "__"
      ) %>%
      mutate(modality = modality_name) %>%
      select(modality, feature, everything())
  }
)

save_table(
  feature_distribution_summary,
  "eda_feature_distribution_summary.csv"
)

modality_distribution_summary <- feature_distribution_summary %>%
  group_by(modality) %>%
  summarise(
    median_cv = median(cv, na.rm = TRUE),
    median_iqr = median(iqr, na.rm = TRUE),
    median_mad = median(mad, na.rm = TRUE),
    median_abs_skewness = median(abs(skewness), na.rm = TRUE),
    .groups = "drop"
  )

save_table(
  modality_distribution_summary,
  "eda_modality_distribution_summary.csv"
)

# Extend modality summary with distribution diagnostics
modality_summary_table <- modality_summary_table %>%
  left_join(modality_distribution_summary, by = "modality")

save_table(
  modality_summary_table,
  "eda_modality_summary_table.csv"
)

# Final script completion message
message("Exploratory analysis complete. Tables saved to outputs/tables and figures saved to outputs/figures/eda.")
