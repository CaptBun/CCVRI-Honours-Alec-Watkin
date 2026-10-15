# Repository stage: medications
# Source supplied as: 03E_medication_confounding_investigation(5).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# 03E Medication confounding and clinical phenotype investigation
#
# Purpose:
#   Investigate whether medication predictors, especially statin and aspirin,
#   are acting as independent predictors, treatment proxies, or surrogate markers
#   for existing cardiovascular/cardiometabolic disease.
#
# Inputs:
#   combined_multimodal_complete_cohort.csv
#
# Outputs:
#   outputs/tables/medication_investigation/03E_*.csv
#   outputs/figures/medication_investigation/03E_*.png
#   outputs/models/03E_medication_confounding_investigation_outputs.rds

library(tidyverse)
library(broom)

set.seed(123)

dir.create("outputs/tables/medication_investigation", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures/medication_investigation", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)

data_file <- ccrvi_config$files$compiled_cohort

medication_vars <- c(
  "med_metformin",
  "med_statin",
  "med_aspirin",
  "med_ppi",
  "med_glp1_gip",
  "med_none"
)

primary_medication_vars <- c(
  "med_statin",
  "med_aspirin"
)

expanded_clinical_vars <- c(
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
  medication_vars
)

cardiometabolic_context_vars <- c(
  "age",
  "gender",
  "bmi",
  "smoking_bin",
  "masld",
  "diabetes_bin",
  "hypertension",
  "dyslipidemia",
  "family_history_cardiac"
)

liver_context_vars <- c(
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

save_table <- function(x, filename) {
  write.csv(
    x,
    file.path("outputs/tables/medication_investigation", filename),
    row.names = FALSE
  )
}

save_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path("outputs/figures/medication_investigation", filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

clean_binary <- function(x) {
  if (is.numeric(x)) {
    return(ifelse(is.na(x), NA_real_, ifelse(x > 0, 1, 0)))
  }

  x_chr <- str_to_lower(str_trim(as.character(x)))

  case_when(
    x_chr %in% c("1", "yes", "y", "true", "t", "present", "positive") ~ 1,
    x_chr %in% c("0", "no", "n", "false", "f", "absent", "negative") ~ 0,
    is.na(x_chr) | x_chr == "" ~ NA_real_,
    TRUE ~ suppressWarnings(as.numeric(x_chr))
  )
}

clean_factor_label <- function(x) {
  replace_na(as.character(x), "Missing")
}

safe_fisher_or_chisq <- function(tbl) {
  if (any(dim(tbl) < 2) || sum(tbl) == 0) {
    return(tibble(test = NA_character_, p_value = NA_real_))
  }

  expected <- suppressWarnings(chisq.test(tbl)$expected)

  if (any(expected < 5)) {
    test_out <- tryCatch(
      fisher.test(tbl),
      error = function(e) NULL
    )

    return(
      tibble(
        test = "Fisher exact",
        p_value = ifelse(is.null(test_out), NA_real_, test_out$p.value)
      )
    )
  }

  test_out <- suppressWarnings(chisq.test(tbl))

  tibble(
    test = "Chi-square",
    p_value = test_out$p.value
  )
}

binary_or_table <- function(data, exposure, outcome = "mace") {
  df <- data %>%
    select(all_of(c(exposure, outcome))) %>%
    filter(!is.na(.data[[exposure]]), !is.na(.data[[outcome]]))

  if (nrow(df) == 0 || n_distinct(df[[exposure]]) < 2 || n_distinct(df[[outcome]]) < 2) {
    return(
      tibble(
        medication = exposure,
        n = nrow(df),
        exposed_n = sum(df[[exposure]] == 1, na.rm = TRUE),
        unexposed_n = sum(df[[exposure]] == 0, na.rm = TRUE),
        exposed_mace_n = NA_integer_,
        exposed_mace_rate = NA_real_,
        unexposed_mace_n = NA_integer_,
        unexposed_mace_rate = NA_real_,
        odds_ratio = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        test = NA_character_,
        p_value = NA_real_
      )
    )
  }

  tbl <- table(exposed = df[[exposure]], mace = df[[outcome]])

  complete_tbl <- matrix(
    0,
    nrow = 2,
    ncol = 2,
    dimnames = list(exposed = c("0", "1"), mace = c("0", "1"))
  )

  complete_tbl[rownames(tbl), colnames(tbl)] <- tbl

  exposed_mace_n <- complete_tbl["1", "1"]
  exposed_nonmace_n <- complete_tbl["1", "0"]
  unexposed_mace_n <- complete_tbl["0", "1"]
  unexposed_nonmace_n <- complete_tbl["0", "0"]

  fisher_out <- tryCatch(
    fisher.test(complete_tbl),
    error = function(e) NULL
  )

  test_tbl <- safe_fisher_or_chisq(complete_tbl)

  tibble(
    medication = exposure,
    n = nrow(df),
    exposed_n = exposed_mace_n + exposed_nonmace_n,
    unexposed_n = unexposed_mace_n + unexposed_nonmace_n,
    exposed_mace_n = exposed_mace_n,
    exposed_mace_rate = exposed_mace_n / (exposed_mace_n + exposed_nonmace_n),
    unexposed_mace_n = unexposed_mace_n,
    unexposed_mace_rate = unexposed_mace_n / (unexposed_mace_n + unexposed_nonmace_n),
    odds_ratio = ifelse(is.null(fisher_out), NA_real_, as.numeric(fisher_out$estimate)),
    conf_low = ifelse(is.null(fisher_out), NA_real_, fisher_out$conf.int[1]),
    conf_high = ifelse(is.null(fisher_out), NA_real_, fisher_out$conf.int[2])
  ) %>%
    bind_cols(test_tbl)
}

cramers_v <- function(x, y) {
  tbl <- table(x, y, useNA = "no")

  if (any(dim(tbl) < 2) || sum(tbl) == 0) {
    return(NA_real_)
  }

  chi <- suppressWarnings(chisq.test(tbl)$statistic)
  n <- sum(tbl)
  k <- min(nrow(tbl), ncol(tbl))

  as.numeric(sqrt(chi / (n * (k - 1))))
}

safe_spearman <- function(x, y) {
  n_complete <- sum(!is.na(x) & !is.na(y))

  if (n_complete < 5) {
    return(tibble(r = NA_real_, p_value = NA_real_, n = n_complete))
  }

  if (length(unique(x[!is.na(x)])) < 2 || length(unique(y[!is.na(y)])) < 2) {
    return(tibble(r = NA_real_, p_value = NA_real_, n = n_complete))
  }

  test_out <- suppressWarnings(
    cor.test(x, y, method = "spearman", exact = FALSE)
  )

  tibble(
    r = as.numeric(test_out$estimate),
    p_value = test_out$p.value,
    n = n_complete
  )
}

variable_type <- function(x) {
  if (is.numeric(x)) {
    non_missing <- x[!is.na(x)]
    unique_vals <- unique(non_missing)

    if (length(unique_vals) <= 2 && all(unique_vals %in% c(0, 1))) {
      return("binary")
    }

    return("numeric")
  }

  non_missing <- unique(x[!is.na(x)])

  if (length(non_missing) <= 2) {
    return("binary_or_two_level")
  }

  "categorical"
}

medication_association <- function(data, med_var, comparison_var) {
  med_numeric <- clean_binary(data[[med_var]])
  x <- data[[comparison_var]]
  type_x <- variable_type(x)

  if (type_x == "numeric") {
    spearman <- safe_spearman(med_numeric, as.numeric(x))

    return(
      tibble(
        medication = med_var,
        variable = comparison_var,
        variable_type = "numeric",
        association_metric = "Spearman rho",
        association_value = spearman$r,
        p_value = spearman$p_value,
        n = spearman$n
      )
    )
  }

  x_factor <- clean_factor_label(x)
  med_factor <- factor(med_numeric, levels = c(0, 1), labels = c("No", "Yes"))
  tbl <- table(med_factor, x_factor, useNA = "no")
  test_tbl <- safe_fisher_or_chisq(tbl)

  tibble(
    medication = med_var,
    variable = comparison_var,
    variable_type = type_x,
    association_metric = "Cramers V",
    association_value = cramers_v(med_factor, x_factor),
    p_value = test_tbl$p_value,
    n = sum(!is.na(med_numeric) & !is.na(x_factor))
  )
}

summarise_numeric_by_med <- function(data, med_var, variable) {
  df <- data %>%
    select(all_of(c(med_var, variable))) %>%
    filter(!is.na(.data[[med_var]]), !is.na(.data[[variable]]))

  if (nrow(df) < 5) {
    return(tibble())
  }

  test_out <- tryCatch(
    wilcox.test(df[[variable]] ~ df[[med_var]]),
    error = function(e) NULL
  )

  df %>%
    group_by(.data[[med_var]]) %>%
    summarise(
      n = n(),
      median = median(.data[[variable]], na.rm = TRUE),
      q1 = quantile(.data[[variable]], 0.25, na.rm = TRUE),
      q3 = quantile(.data[[variable]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      medication = med_var,
      variable = variable,
      variable_type = "numeric",
      medication_present = .data[[med_var]],
      p_value = ifelse(is.null(test_out), NA_real_, test_out$p.value)
    ) %>%
    select(medication, variable, variable_type, medication_present, everything(), -all_of(med_var))
}

summarise_categorical_by_med <- function(data, med_var, variable) {
  df <- data %>%
    select(all_of(c(med_var, variable))) %>%
    mutate(variable_value = clean_factor_label(.data[[variable]])) %>%
    filter(!is.na(.data[[med_var]]))

  tbl <- table(df[[med_var]], df$variable_value)
  test_tbl <- safe_fisher_or_chisq(tbl)

  df %>%
    count(.data[[med_var]], variable_value) %>%
    group_by(.data[[med_var]]) %>%
    mutate(
      group_total = sum(n),
      proportion = n / group_total
    ) %>%
    ungroup() %>%
    mutate(
      medication = med_var,
      variable = variable,
      variable_type = "categorical",
      medication_present = .data[[med_var]],
      p_value = test_tbl$p_value,
      test = test_tbl$test
    ) %>%
    select(
      medication,
      variable,
      variable_type,
      medication_present,
      variable_value,
      n,
      group_total,
      proportion,
      test,
      p_value,
      -all_of(med_var)
    )
}

# Load and prepare data
model_source <- read.csv(data_file) %>%
  mutate(mace = as.numeric(mace)) %>%
  filter(!is.na(mace))

medication_vars <- intersect(medication_vars, names(model_source))
primary_medication_vars <- intersect(primary_medication_vars, names(model_source))
expanded_clinical_vars <- intersect(expanded_clinical_vars, names(model_source))
cardiometabolic_context_vars <- intersect(cardiometabolic_context_vars, names(model_source))
liver_context_vars <- intersect(liver_context_vars, names(model_source))

pathology_vars <- names(model_source) %>%
  str_subset("^path_") %>%
  str_subset("sample|tube|plate|source|id|label|timepoint|visit", negate = TRUE)

if (length(medication_vars) == 0) {
  stop("No medication variables were found in the dataset.")
}

model_source <- model_source %>%
  mutate(across(all_of(medication_vars), clean_binary))

medication_label_table <- tibble(
  medication = medication_vars,
  medication_label = medication %>%
    str_remove("^med_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
)

# Medication prevalence and MACE enrichment
medication_prevalence <- model_source %>%
  summarise(
    n_total = n(),
    n_mace = sum(mace == 1, na.rm = TRUE),
    mace_rate = mean(mace == 1, na.rm = TRUE),
    across(
      all_of(medication_vars),
      list(
        n = ~ sum(.x == 1, na.rm = TRUE),
        prevalence = ~ mean(.x == 1, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

save_table(medication_prevalence, "03E_medication_prevalence_overall.csv")

medication_mace_summary <- map_dfr(
  medication_vars,
  ~ binary_or_table(model_source, exposure = .x, outcome = "mace")
) %>%
  left_join(medication_label_table, by = "medication") %>%
  arrange(desc(exposed_mace_rate), desc(odds_ratio))

save_table(medication_mace_summary, "03E_medication_mace_enrichment.csv")

medication_mace_long <- model_source %>%
  select(mace, all_of(medication_vars)) %>%
  pivot_longer(
    cols = all_of(medication_vars),
    names_to = "medication",
    values_to = "medication_present"
  ) %>%
  filter(!is.na(medication_present)) %>%
  mutate(
    mace_status = if_else(mace == 1, "MACE", "No MACE"),
    medication_present = factor(
      medication_present,
      levels = c(0, 1),
      labels = c("Not prescribed", "Prescribed")
    )
  ) %>%
  count(medication, mace_status, medication_present) %>%
  group_by(medication, mace_status) %>%
  mutate(
    group_total = sum(n),
    proportion = n / group_total
  ) %>%
  ungroup() %>%
  left_join(medication_label_table, by = "medication")

save_table(medication_mace_long, "03E_medication_mace_status_counts_long.csv")

medication_mace_plot <- medication_mace_long %>%
  filter(medication_present == "Prescribed") %>%
  ggplot(aes(x = medication_label, y = proportion, fill = mace_status)) +
  geom_col(position = position_dodge(width = 0.75)) +
  coord_flip() +
  labs(
    title = "Medication prevalence by MACE status",
    subtitle = "Proportion prescribed within MACE and non-MACE groups",
    x = NULL,
    y = "Proportion prescribed",
    fill = "Outcome"
  ) +
  theme_minimal(base_size = 12)

save_plot(medication_mace_plot, "03E_medication_prevalence_by_mace_status.png", width = 8, height = 5)

# Medication overlap
medication_overlap_patient_table <- model_source %>%
  mutate(
    medication_pattern = pmap_chr(
      select(., all_of(medication_vars)),
      function(...) {
        values <- c(...)
        present <- medication_vars[as.numeric(values) == 1]

        if (length(present) == 0) {
          return("No listed medications")
        }

        present %>%
          str_remove("^med_") %>%
          str_replace_all("_", " ") %>%
          str_to_title() %>%
          paste(collapse = " + ")
      }
    )
  ) %>%
  count(medication_pattern, mace, sort = TRUE) %>%
  group_by(medication_pattern) %>%
  mutate(
    pattern_total = sum(n),
    proportion_within_pattern = n / pattern_total
  ) %>%
  ungroup() %>%
  arrange(desc(pattern_total), desc(mace))

save_table(medication_overlap_patient_table, "03E_medication_overlap_patterns_by_mace.csv")

medication_pairwise_overlap <- expand_grid(
  medication_1 = medication_vars,
  medication_2 = medication_vars
) %>%
  mutate(
    both_present_n = map2_int(
      medication_1,
      medication_2,
      ~ sum(model_source[[.x]] == 1 & model_source[[.y]] == 1, na.rm = TRUE)
    ),
    med1_present_n = map_int(
      medication_1,
      ~ sum(model_source[[.x]] == 1, na.rm = TRUE)
    ),
    med2_present_n = map_int(
      medication_2,
      ~ sum(model_source[[.x]] == 1, na.rm = TRUE)
    ),
    overlap_proportion_of_med1 = ifelse(
      med1_present_n > 0,
      both_present_n / med1_present_n,
      NA_real_
    )
  ) %>%
  left_join(
    medication_label_table %>% rename(medication_1_label = medication_label),
    by = c("medication_1" = "medication")
  ) %>%
  left_join(
    medication_label_table %>% rename(medication_2_label = medication_label),
    by = c("medication_2" = "medication")
  )

save_table(medication_pairwise_overlap, "03E_medication_pairwise_overlap.csv")

medication_overlap_plot <- medication_pairwise_overlap %>%
  ggplot(aes(x = medication_1_label, y = medication_2_label, fill = overlap_proportion_of_med1)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * overlap_proportion_of_med1)), size = 3) +
  coord_equal() +
  labs(
    title = "Medication overlap",
    subtitle = "Cell value = proportion of row medication users also using column medication",
    x = "Row medication",
    y = "Column medication",
    fill = "Overlap"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(medication_overlap_plot, "03E_medication_pairwise_overlap_heatmap.png", width = 7, height = 6)

# Medication by liver disease aetiology and MACE
aetiology_var <- if ("primary_aetiology" %in% names(model_source)) {
  "primary_aetiology"
} else if ("aetiology" %in% names(model_source)) {
  "aetiology"
} else {
  NA_character_
}

if (!is.na(aetiology_var)) {
  medication_aetiology_mace <- model_source %>%
    select(mace, all_of(aetiology_var), all_of(medication_vars)) %>%
    mutate(
      aetiology = clean_factor_label(.data[[aetiology_var]]),
      mace_status = if_else(mace == 1, "MACE", "No MACE")
    ) %>%
    pivot_longer(
      cols = all_of(medication_vars),
      names_to = "medication",
      values_to = "medication_present"
    ) %>%
    filter(!is.na(medication_present)) %>%
    mutate(
      medication_status = if_else(medication_present == 1, "Prescribed", "Not prescribed")
    ) %>%
    count(medication, aetiology, mace_status, medication_status) %>%
    group_by(medication, aetiology, medication_status) %>%
    mutate(
      total = sum(n),
      mace_rate = sum(n[mace_status == "MACE"]) / total
    ) %>%
    ungroup() %>%
    left_join(medication_label_table, by = "medication")

  save_table(medication_aetiology_mace, "03E_medication_by_aetiology_and_mace.csv")

  medication_aetiology_mace_plot <- medication_aetiology_mace %>%
    filter(medication %in% primary_medication_vars, medication_status == "Prescribed") %>%
    ggplot(aes(x = aetiology, y = n, fill = mace_status)) +
    geom_col(position = "stack") +
    facet_wrap(~ medication_label) +
    coord_flip() +
    labs(
      title = "MACE events among statin and aspirin users by liver disease aetiology",
      x = "Primary liver disease aetiology",
      y = "Patient count",
      fill = "Outcome"
    ) +
    theme_minimal(base_size = 12)

  save_plot(medication_aetiology_mace_plot, "03E_primary_medications_mace_by_aetiology.png", width = 9, height = 5)
} else {
  medication_aetiology_mace <- tibble()
}

# Medication associations with expanded clinical variables
clinical_association_vars <- expanded_clinical_vars %>%
  setdiff(medication_vars) %>%
  intersect(names(model_source))

medication_clinical_associations <- expand_grid(
  medication = medication_vars,
  variable = clinical_association_vars
) %>%
  pmap_dfr(function(medication, variable) {
    medication_association(model_source, medication, variable)
  }) %>%
  left_join(medication_label_table, by = "medication") %>%
  arrange(medication, desc(abs(association_value)))

save_table(medication_clinical_associations, "03E_medication_associations_expanded_clinical.csv")

primary_medication_clinical_association_plot <- medication_clinical_associations %>%
  filter(medication %in% primary_medication_vars) %>%
  group_by(medication) %>%
  slice_max(order_by = abs(association_value), n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(variable = factor(variable, levels = rev(unique(variable)))) %>%
  ggplot(aes(x = variable, y = association_value, fill = variable_type)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ medication_label, scales = "free_y") +
  labs(
    title = "Clinical variables most associated with statin and aspirin use",
    subtitle = "Spearman rho for numeric variables; Cramer's V for categorical variables",
    x = NULL,
    y = "Association strength",
    fill = "Variable type"
  ) +
  theme_minimal(base_size = 11)

save_plot(primary_medication_clinical_association_plot, "03E_primary_medications_clinical_associations.png", width = 10, height = 7)

# Medication associations with pathology variables
if (length(pathology_vars) > 0) {
  medication_pathology_associations <- expand_grid(
    medication = medication_vars,
    variable = pathology_vars
  ) %>%
    pmap_dfr(function(medication, variable) {
      medication_association(model_source, medication, variable)
    }) %>%
    left_join(medication_label_table, by = "medication") %>%
    arrange(medication, desc(abs(association_value)))

  save_table(medication_pathology_associations, "03E_medication_associations_pathology.csv")

  primary_medication_pathology_plot <- medication_pathology_associations %>%
    filter(medication %in% primary_medication_vars) %>%
    group_by(medication) %>%
    slice_max(order_by = abs(association_value), n = 20, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      variable_clean = variable %>%
        str_remove("^path_") %>%
        str_replace_all("_", " "),
      variable_clean = factor(variable_clean, levels = rev(unique(variable_clean)))
    ) %>%
    ggplot(aes(x = variable_clean, y = association_value, fill = medication_label)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    facet_wrap(~ medication_label, scales = "free_y") +
    labs(
      title = "Pathology variables most associated with statin and aspirin use",
      subtitle = "Spearman rho for numeric variables; Cramer's V for categorical variables",
      x = NULL,
      y = "Association strength"
    ) +
    theme_minimal(base_size = 11)

  save_plot(primary_medication_pathology_plot, "03E_primary_medications_pathology_associations.png", width = 10, height = 7)
} else {
  medication_pathology_associations <- tibble()
}

# Focused medication phenotype tables
phenotype_vars <- unique(c(
  cardiometabolic_context_vars,
  liver_context_vars,
  "mace",
  pathology_vars
)) %>%
  intersect(names(model_source))

medication_phenotype_numeric <- expand_grid(
  medication = primary_medication_vars,
  variable = phenotype_vars
) %>%
  filter(map_lgl(variable, ~ is.numeric(model_source[[.x]]) && !.x %in% c("mace", medication_vars))) %>%
  pmap_dfr(function(medication, variable) {
    summarise_numeric_by_med(model_source, medication, variable)
  }) %>%
  left_join(medication_label_table, by = "medication")

medication_phenotype_categorical <- expand_grid(
  medication = primary_medication_vars,
  variable = phenotype_vars
) %>%
  filter(map_lgl(variable, ~ !is.numeric(model_source[[.x]]) || .x %in% c("mace", medication_vars))) %>%
  pmap_dfr(function(medication, variable) {
    summarise_categorical_by_med(model_source, medication, variable)
  }) %>%
  left_join(medication_label_table, by = "medication")

save_table(medication_phenotype_numeric, "03E_primary_medication_numeric_phenotype_table.csv")
save_table(medication_phenotype_categorical, "03E_primary_medication_categorical_phenotype_table.csv")

# Key binary context summary for statin/aspirin
key_binary_context_vars <- c(
  "mace",
  "masld",
  "diabetes_bin",
  "hypertension",
  "dyslipidemia",
  "family_history_cardiac",
  "portal_hypertension",
  "thrombocytopenia",
  "tips",
  "decomp_ascites",
  "decomp_encephalopathy",
  "decomp_variceal_bleed",
  "decomp_child_pugh7"
) %>%
  intersect(names(model_source))

key_medication_context_summary <- expand_grid(
  medication = primary_medication_vars,
  variable = key_binary_context_vars
) %>%
  pmap_dfr(function(medication, variable) {
    df <- model_source %>%
      select(all_of(c(medication, variable))) %>%
      mutate(variable_binary = clean_binary(.data[[variable]])) %>%
      filter(!is.na(.data[[medication]]), !is.na(variable_binary))

    if (nrow(df) == 0) {
      return(tibble())
    }

    tbl <- table(df[[medication]], df$variable_binary)
    test_tbl <- safe_fisher_or_chisq(tbl)

    df %>%
      group_by(.data[[medication]]) %>%
      summarise(
        n = n(),
        variable_positive_n = sum(variable_binary == 1, na.rm = TRUE),
        variable_positive_rate = mean(variable_binary == 1, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        medication = medication,
        variable = variable,
        medication_present = .data[[medication]],
        test = test_tbl$test,
        p_value = test_tbl$p_value
      ) %>%
      select(
        medication,
        variable,
        medication_present,
        n,
        variable_positive_n,
        variable_positive_rate,
        test,
        p_value,
        -all_of(medication)
      )
  }) %>%
  left_join(medication_label_table, by = "medication") %>%
  arrange(medication, variable, desc(medication_present))

save_table(key_medication_context_summary, "03E_key_context_variables_by_primary_medication.csv")

# Interpretation summary
strongest_clinical_associations <- medication_clinical_associations %>%
  filter(medication %in% primary_medication_vars) %>%
  group_by(medication) %>%
  slice_max(order_by = abs(association_value), n = 8, with_ties = FALSE) %>%
  summarise(
    strongest_clinical_associations = paste(
      paste0(variable, " (", round(association_value, 2), ")"),
      collapse = "; "
    ),
    .groups = "drop"
  )

strongest_pathology_associations <- if (nrow(medication_pathology_associations) > 0) {
  medication_pathology_associations %>%
    filter(medication %in% primary_medication_vars) %>%
    group_by(medication) %>%
    slice_max(order_by = abs(association_value), n = 8, with_ties = FALSE) %>%
    summarise(
      strongest_pathology_associations = paste(
        paste0(variable, " (", round(association_value, 2), ")"),
        collapse = "; "
      ),
      .groups = "drop"
    )
} else {
  tibble(
    medication = primary_medication_vars,
    strongest_pathology_associations = NA_character_
  )
}

medication_interpretation_summary <- medication_mace_summary %>%
  filter(medication %in% primary_medication_vars) %>%
  select(
    medication,
    medication_label,
    exposed_n,
    exposed_mace_rate,
    unexposed_mace_rate,
    odds_ratio,
    conf_low,
    conf_high,
    p_value
  ) %>%
  left_join(strongest_clinical_associations, by = "medication") %>%
  left_join(strongest_pathology_associations, by = "medication") %>%
  mutate(
    likely_interpretation = case_when(
      medication == "med_statin" ~ "Likely cardiometabolic/CVD treatment proxy; interpret as prescribed risk-marker rather than causal protective effect.",
      medication == "med_aspirin" ~ "Likely secondary-prevention or prior vascular disease proxy; investigate pre-existing CVD indication if available.",
      TRUE ~ "Medication variable may reflect treatment indication and underlying clinical phenotype."
    )
  )

save_table(medication_interpretation_summary, "03E_primary_medication_interpretation_summary.csv")

# Presentation-focused figures
medication_or_plot <- medication_mace_summary %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(
    medication_label = factor(
      medication_label,
      levels = medication_label[order(odds_ratio)]
    )
  ) %>%
  ggplot(
    aes(
      x = medication_label,
      y = odds_ratio,
      ymin = conf_low,
      ymax = conf_high
    )
  ) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  scale_y_log10() +
  labs(
    title = "Medication enrichment among patients with MACE",
    subtitle = "Odds ratios compare medication users with non-users",
    x = NULL,
    y = "Odds ratio for MACE (log scale)"
  ) +
  theme_minimal(base_size = 12)

save_plot(medication_or_plot, "03E_medication_mace_odds_ratios.png", width = 8, height = 5)

key_context_heatmap <- key_medication_context_summary %>%
  filter(medication_present == 1) %>%
  mutate(
    variable_label = variable %>%
      str_replace_all("_", " ") %>%
      str_to_title()
  ) %>%
  ggplot(
    aes(
      x = medication_label,
      y = variable_label,
      fill = variable_positive_rate
    )
  ) +
  geom_tile(colour = "white") +
  geom_text(
    aes(label = sprintf("%.0f%%", 100 * variable_positive_rate)),
    size = 3
  ) +
  labs(
    title = "Clinical phenotype of statin and aspirin users",
    subtitle = "Cell value = prevalence of clinical feature among medication users",
    x = "Medication",
    y = "Clinical feature",
    fill = "Prevalence"
  ) +
  theme_minimal(base_size = 11)

save_plot(key_context_heatmap, "03E_primary_medication_clinical_context_heatmap.png", width = 7, height = 7)

# Save all outputs
medication_investigation_outputs <- list(
  settings = list(
    medication_vars = medication_vars,
    primary_medication_vars = primary_medication_vars,
    expanded_clinical_vars = expanded_clinical_vars,
    pathology_vars = pathology_vars
  ),
  medication_prevalence = medication_prevalence,
  medication_mace_summary = medication_mace_summary,
  medication_mace_long = medication_mace_long,
  medication_overlap_patient_table = medication_overlap_patient_table,
  medication_pairwise_overlap = medication_pairwise_overlap,
  medication_aetiology_mace = medication_aetiology_mace,
  medication_clinical_associations = medication_clinical_associations,
  medication_pathology_associations = medication_pathology_associations,
  medication_phenotype_numeric = medication_phenotype_numeric,
  medication_phenotype_categorical = medication_phenotype_categorical,
  key_medication_context_summary = key_medication_context_summary,
  medication_interpretation_summary = medication_interpretation_summary
)

saveRDS(
  medication_investigation_outputs,
  "outputs/models/03E_medication_confounding_investigation_outputs.rds"
)

message("Medication investigation complete.")
message("Tables saved to outputs/tables/medication_investigation/")
message("Figures saved to outputs/figures/medication_investigation/")
message("Primary medication interpretation summary:")
print(medication_interpretation_summary)






#the following is a secondary script section that includes summary tables and descriptions about medication data 
#that has been imported from the former 3G script
# Inputs:
#   combined_multimodal_complete_cohort.csv
#
# Outputs:
#   outputs/tables/medication_summary/03G_*.csv
#   outputs/figures/medication_summary/03G_*.png

library(tidyverse)

set.seed(123)

# Create output folders
dir.create(
  "outputs/tables/medication_summary",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/figures/medication_summary",
  recursive = TRUE,
  showWarnings = FALSE
)

# Settings
data_file <- ccrvi_config$files$compiled_cohort

medication_vars <- c(
  "med_metformin",
  "med_statin",
  "med_aspirin",
  "med_ppi",
  "med_glp1_gip",
  "med_none"
)

primary_medications <- c(
  "med_statin",
  "med_aspirin"
)

# Helper functions
save_table <- function(x, filename) {
  write.csv(
    x,
    file.path("outputs/tables/medication_summary", filename),
    row.names = FALSE
  )
}

save_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path("outputs/figures/medication_summary", filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

clean_binary <- function(x) {
  if (is.numeric(x)) {
    return(ifelse(is.na(x), NA_real_, ifelse(x > 0, 1, 0)))
  }
  
  x_chr <- str_to_lower(str_trim(as.character(x)))
  
  case_when(
    x_chr %in% c("1", "yes", "y", "true", "t", "present", "positive") ~ 1,
    x_chr %in% c("0", "no", "n", "false", "f", "absent", "negative") ~ 0,
    is.na(x_chr) | x_chr == "" ~ NA_real_,
    TRUE ~ suppressWarnings(as.numeric(x_chr))
  )
}

safe_test <- function(tbl) {
  if (any(dim(tbl) < 2) || sum(tbl) == 0) {
    return(tibble(test = NA_character_, p_value = NA_real_))
  }
  
  expected <- suppressWarnings(chisq.test(tbl)$expected)
  
  if (any(expected < 5)) {
    result <- tryCatch(fisher.test(tbl), error = function(e) NULL)
    
    return(
      tibble(
        test = "Fisher exact",
        p_value = ifelse(is.null(result), NA_real_, result$p.value)
      )
    )
  }
  
  result <- suppressWarnings(chisq.test(tbl))
  
  tibble(
    test = "Chi-square",
    p_value = result$p.value
  )
}

summarise_medication_mace <- function(data, medication) {
  df <- data %>%
    select(mace, all_of(medication)) %>%
    filter(
      !is.na(mace),
      !is.na(.data[[medication]])
    )
  
  users <- df %>% filter(.data[[medication]] == 1)
  non_users <- df %>% filter(.data[[medication]] == 0)
  
  contingency <- table(
    medication = df[[medication]],
    mace = df$mace
  )
  
  complete_contingency <- matrix(
    0,
    nrow = 2,
    ncol = 2,
    dimnames = list(
      medication = c("0", "1"),
      mace = c("0", "1")
    )
  )
  
  complete_contingency[
    rownames(contingency),
    colnames(contingency)
  ] <- contingency
  
  fisher_result <- tryCatch(
    fisher.test(complete_contingency),
    error = function(e) NULL
  )
  
  test_result <- safe_test(complete_contingency)
  
  user_mace_rate <- if (nrow(users) > 0) {
    mean(users$mace == 1, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  non_user_mace_rate <- if (nrow(non_users) > 0) {
    mean(non_users$mace == 1, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  tibble(
    medication = medication,
    total_cohort_n = nrow(df),
    
    users_n = nrow(users),
    users_percent = 100 * nrow(users) / nrow(df),
    users_mace_n = sum(users$mace == 1, na.rm = TRUE),
    users_nonmace_n = sum(users$mace == 0, na.rm = TRUE),
    users_mace_rate = user_mace_rate,
    
    non_users_n = nrow(non_users),
    non_users_percent = 100 * nrow(non_users) / nrow(df),
    non_users_mace_n = sum(non_users$mace == 1, na.rm = TRUE),
    non_users_nonmace_n = sum(non_users$mace == 0, na.rm = TRUE),
    non_users_mace_rate = non_user_mace_rate,
    
    absolute_mace_rate_difference =
      user_mace_rate - non_user_mace_rate,
    
    risk_ratio = ifelse(
      is.na(non_user_mace_rate) || non_user_mace_rate == 0,
      NA_real_,
      user_mace_rate / non_user_mace_rate
    ),
    
    odds_ratio = ifelse(
      is.null(fisher_result),
      NA_real_,
      as.numeric(fisher_result$estimate)
    ),
    
    odds_ratio_ci_low = ifelse(
      is.null(fisher_result),
      NA_real_,
      fisher_result$conf.int[1]
    ),
    
    odds_ratio_ci_high = ifelse(
      is.null(fisher_result),
      NA_real_,
      fisher_result$conf.int[2]
    )
  ) %>%
    bind_cols(test_result)
}

# Load and prepare data
model_source <- read.csv(
  data_file,
  check.names = FALSE
) %>%
  mutate(mace = as.numeric(mace)) %>%
  filter(!is.na(mace))

medication_vars <- intersect(
  medication_vars,
  names(model_source)
)

primary_medications <- intersect(
  primary_medications,
  names(model_source)
)

if (length(medication_vars) == 0) {
  stop("No medication variables were found in the cohort.")
}

model_source <- model_source %>%
  mutate(
    across(
      all_of(medication_vars),
      clean_binary
    )
  )

medication_label_table <- tibble(
  medication = medication_vars,
  medication_label = medication %>%
    str_remove("^med_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
)

# Cohort summary
cohort_summary <- tibble(
  total_cohort_n = nrow(model_source),
  mace_n = sum(model_source$mace == 1, na.rm = TRUE),
  nonmace_n = sum(model_source$mace == 0, na.rm = TRUE),
  mace_rate = mean(model_source$mace == 1, na.rm = TRUE),
  mace_percent = 100 * mace_rate
)

save_table(
  cohort_summary,
  "03G_total_cohort_summary.csv"
)

# Medication prevalence in the total cohort
medication_prevalence_table <- map_dfr(
  medication_vars,
  function(medication) {
    non_missing_n <- sum(
      !is.na(model_source[[medication]])
    )
    
    users_n <- sum(
      model_source[[medication]] == 1,
      na.rm = TRUE
    )
    
    tibble(
      medication = medication,
      total_cohort_n = nrow(model_source),
      non_missing_n = non_missing_n,
      users_n = users_n,
      users_percent_of_total = 100 * users_n / nrow(model_source),
      users_percent_of_non_missing = 100 * users_n / non_missing_n,
      non_users_n = sum(
        model_source[[medication]] == 0,
        na.rm = TRUE
      ),
      missing_n = sum(
        is.na(model_source[[medication]])
      )
    )
  }
) %>%
  left_join(
    medication_label_table,
    by = "medication"
  ) %>%
  select(
    medication,
    medication_label,
    total_cohort_n,
    non_missing_n,
    users_n,
    users_percent_of_total,
    users_percent_of_non_missing,
    non_users_n,
    missing_n
  ) %>%
  arrange(desc(users_percent_of_total))

save_table(
  medication_prevalence_table,
  "03G_medication_prevalence_total_cohort.csv"
)

# MACE rate among users and non-users
medication_mace_summary <- map_dfr(
  medication_vars,
  ~ summarise_medication_mace(
    data = model_source,
    medication = .x
  )
) %>%
  left_join(
    medication_label_table,
    by = "medication"
  ) %>%
  arrange(desc(users_mace_rate))

save_table(
  medication_mace_summary,
  "03G_mace_rate_by_medication.csv"
)

# Presentation-ready medication table
medication_presentation_table <- medication_mace_summary %>%
  transmute(
    Medication = medication_label,
    `Total cohort, n` = total_cohort_n,
    `Users, n` = users_n,
    `Users, % of cohort` = round(users_percent, 1),
    `MACE among users, n` = users_mace_n,
    `MACE rate among users, %` = round(
      100 * users_mace_rate,
      1
    ),
    `Non-users, n` = non_users_n,
    `MACE among non-users, n` = non_users_mace_n,
    `MACE rate among non-users, %` = round(
      100 * non_users_mace_rate,
      1
    ),
    `Absolute difference, percentage points` = round(
      100 * absolute_mace_rate_difference,
      1
    ),
    `Risk ratio` = round(risk_ratio, 2),
    `Odds ratio` = round(odds_ratio, 2),
    `Odds ratio 95% CI` = paste0(
      round(odds_ratio_ci_low, 2),
      "–",
      round(odds_ratio_ci_high, 2)
    ),
    `p-value` = signif(p_value, 3)
  )

save_table(
  medication_presentation_table,
  "03G_medication_mace_presentation_table.csv"
)

# Focused statin and aspirin summary
statin_aspirin_summary <- medication_mace_summary %>%
  filter(
    medication %in% primary_medications
  ) %>%
  arrange(
    match(
      medication,
      primary_medications
    )
  )

save_table(
  statin_aspirin_summary,
  "03G_statin_aspirin_mace_summary.csv"
)

statin_aspirin_presentation_table <- statin_aspirin_summary %>%
  transmute(
    Medication = medication_label,
    `Total cohort, n` = total_cohort_n,
    `Users, n` = users_n,
    `Users, % of cohort` = round(users_percent, 1),
    `MACE among users, n` = users_mace_n,
    `MACE rate among users, %` = round(
      100 * users_mace_rate,
      1
    ),
    `Non-users, n` = non_users_n,
    `MACE among non-users, n` = non_users_mace_n,
    `MACE rate among non-users, %` = round(
      100 * non_users_mace_rate,
      1
    ),
    `Risk ratio` = round(risk_ratio, 2),
    `Odds ratio` = round(odds_ratio, 2),
    `Odds ratio 95% CI` = paste0(
      round(odds_ratio_ci_low, 2),
      "–",
      round(odds_ratio_ci_high, 2)
    ),
    `p-value` = signif(p_value, 3)
  )

save_table(
  statin_aspirin_presentation_table,
  "03G_statin_aspirin_presentation_table.csv"
)

# Statin and aspirin overlap
if (all(c("med_statin", "med_aspirin") %in% names(model_source))) {
  statin_aspirin_overlap <- model_source %>%
    mutate(
      statin_status = if_else(
        med_statin == 1,
        "Statin",
        "No statin"
      ),
      aspirin_status = if_else(
        med_aspirin == 1,
        "Aspirin",
        "No aspirin"
      )
    ) %>%
    group_by(
      statin_status,
      aspirin_status
    ) %>%
    summarise(
      total_n = n(),
      mace_n = sum(mace == 1, na.rm = TRUE),
      nonmace_n = sum(mace == 0, na.rm = TRUE),
      mace_rate = mean(mace == 1, na.rm = TRUE),
      mace_percent = 100 * mace_rate,
      .groups = "drop"
    ) %>%
    arrange(desc(mace_rate))
  
  save_table(
    statin_aspirin_overlap,
    "03G_statin_aspirin_overlap_summary.csv"
  )
} else {
  statin_aspirin_overlap <- tibble()
}

# Medication use by liver disease aetiology
aetiology_var <- if ("primary_aetiology" %in% names(model_source)) {
  "primary_aetiology"
} else if ("aetiology" %in% names(model_source)) {
  "aetiology"
} else {
  NA_character_
}

if (!is.na(aetiology_var)) {
  medication_by_aetiology <- model_source %>%
    mutate(
      aetiology = replace_na(
        as.character(.data[[aetiology_var]]),
        "Missing"
      )
    ) %>%
    select(
      aetiology,
      mace,
      all_of(medication_vars)
    ) %>%
    pivot_longer(
      cols = all_of(medication_vars),
      names_to = "medication",
      values_to = "medication_present"
    ) %>%
    filter(!is.na(medication_present)) %>%
    group_by(
      aetiology,
      medication
    ) %>%
    summarise(
      aetiology_total_n = n(),
      users_n = sum(
        medication_present == 1,
        na.rm = TRUE
      ),
      users_percent = 100 * mean(
        medication_present == 1,
        na.rm = TRUE
      ),
      users_mace_n = sum(
        medication_present == 1 &
          mace == 1,
        na.rm = TRUE
      ),
      users_mace_rate = ifelse(
        users_n > 0,
        users_mace_n / users_n,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    left_join(
      medication_label_table,
      by = "medication"
    ) %>%
    arrange(
      medication,
      desc(users_percent)
    )
  
  save_table(
    medication_by_aetiology,
    "03G_medication_by_liver_aetiology.csv"
  )
  
  save_table(
    medication_by_aetiology %>%
      filter(
        medication %in% primary_medications
      ),
    "03G_statin_aspirin_by_liver_aetiology.csv"
  )
} else {
  medication_by_aetiology <- tibble()
}

# Figures
medication_prevalence_plot <- medication_prevalence_table %>%
  mutate(
    medication_label = factor(
      medication_label,
      levels = rev(
        medication_label[
          order(users_percent_of_total)
        ]
      )
    )
  ) %>%
  ggplot(
    aes(
      x = medication_label,
      y = users_percent_of_total
    )
  ) +
  geom_col() +
  geom_text(
    aes(
      label = paste0(
        users_n,
        " (",
        round(users_percent_of_total, 1),
        "%)"
      )
    ),
    hjust = -0.05,
    size = 3.5
  ) +
  coord_flip() +
  expand_limits(
    y = max(
      medication_prevalence_table$users_percent_of_total,
      na.rm = TRUE
    ) * 1.2
  ) +
  labs(
    title = "Medication use in the complete cohort",
    x = NULL,
    y = "Percentage of cohort"
  ) +
  theme_minimal(base_size = 12)

save_plot(
  medication_prevalence_plot,
  "03G_medication_prevalence_total_cohort.png",
  width = 8,
  height = 5
)

mace_rate_plot <- medication_mace_summary %>%
  select(
    medication_label,
    users_mace_rate,
    non_users_mace_rate
  ) %>%
  pivot_longer(
    cols = c(
      users_mace_rate,
      non_users_mace_rate
    ),
    names_to = "group",
    values_to = "mace_rate"
  ) %>%
  mutate(
    group = recode(
      group,
      users_mace_rate = "Medication users",
      non_users_mace_rate = "Non-users"
    )
  ) %>%
  ggplot(
    aes(
      x = medication_label,
      y = 100 * mace_rate,
      fill = group
    )
  ) +
  geom_col(
    position = position_dodge(width = 0.75)
  ) +
  coord_flip() +
  labs(
    title = "MACE rate by medication status",
    x = NULL,
    y = "MACE rate (%)",
    fill = "Group"
  ) +
  theme_minimal(base_size = 12)

save_plot(
  mace_rate_plot,
  "03G_mace_rate_by_medication.png",
  width = 9,
  height = 6
)

if (nrow(statin_aspirin_overlap) > 0) {
  statin_aspirin_overlap_plot <- statin_aspirin_overlap %>%
    unite(
      "combination",
      statin_status,
      aspirin_status,
      sep = " + ",
      remove = FALSE
    ) %>%
    mutate(
      combination = factor(
        combination,
        levels = combination[order(mace_rate)]
      )
    ) %>%
    ggplot(
      aes(
        x = combination,
        y = mace_percent
      )
    ) +
    geom_col() +
    geom_text(
      aes(
        label = paste0(
          mace_n,
          "/",
          total_n,
          "\n",
          round(mace_percent, 1),
          "%"
        )
      ),
      hjust = -0.05,
      size = 3.5
    ) +
    coord_flip() +
    expand_limits(
      y = max(
        statin_aspirin_overlap$mace_percent,
        na.rm = TRUE
      ) * 1.2
    ) +
    labs(
      title = "MACE rate by statin and aspirin combination",
      x = NULL,
      y = "MACE rate (%)"
    ) +
    theme_minimal(base_size = 12)
  
  save_plot(
    statin_aspirin_overlap_plot,
    "03G_statin_aspirin_overlap_mace_rate.png",
    width = 8,
    height = 5
  )
}

# Save combined outputs
medication_summary_outputs <- list(
  cohort_summary = cohort_summary,
  medication_prevalence_table = medication_prevalence_table,
  medication_mace_summary = medication_mace_summary,
  medication_presentation_table = medication_presentation_table,
  statin_aspirin_summary = statin_aspirin_summary,
  statin_aspirin_presentation_table = statin_aspirin_presentation_table,
  statin_aspirin_overlap = statin_aspirin_overlap,
  medication_by_aetiology = medication_by_aetiology
)

saveRDS(
  medication_summary_outputs,
  "outputs/tables/medication_summary/03G_medication_summary_outputs.rds"
)

message("Medication descriptive summary complete.")
message("Tables saved to outputs/tables/medication_summary/")
message("Figures saved to outputs/figures/medication_summary/")
message("Statin and aspirin summary:")
print(statin_aspirin_presentation_table)
