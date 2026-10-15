# Repository stage: compile
# Source supplied as: 01_compile_multimodal_from_new_workbook_v2_bileacid_fix(8).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# Compile clinical, outcome, omics, bile acid, and pathology data
#
# Purpose
# This script builds the patient-level modelling cohort for the cirrhosis CVD
# risk project using:
#   1. Clinical REDCap-style CSV data
#   2. Outcome CSV data
#   3. One Excel workbook containing all omics, bile acids, and pathology data
#
# The new workbook contains several technical identifier columns. This script
# deliberately avoids automatic ID detection for omics/pathology sheets and
# instead uses an explicit sheet-by-sheet ID map. This is important because
# columns such as sampleID, plateID, tubeLabel, Sample, and lipoprotein variables
# containing "ID" are not patient identifiers.
#
# Expected files in the working directory:
#   clindata.csv
#   outcomedata.csv
#   Alec Honours.xlsx
#
# Load packages
library(tidyverse)
library(readxl)

set.seed(123)

# File paths
clinical_file <- ccrvi_config$files$clinical
outcome_file <- ccrvi_config$files$outcome
omics_workbook <- ccrvi_config$files$workbook

# Output folders
# Audit files are saved
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/audit", recursive = TRUE, showWarnings = FALSE)

# Helper functions

clean_id <- function(x) {
  # Standardise patient identifiers across files.
  #
  # Several omics sheets use IDs such as "ALF001", "alf001", or "alf001-1".
  # These represent the same patient.
  x %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower() %>%
    str_remove("-1$")
}

first_nonmissing <- function(x) {
  # Return the first observed non-empty value.
  #
  # This is used for clinical fields after restricting to enrolment records,
  # where repeated rows can occur but patient-level baseline values are desired.
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA else x[1]
}

clean_column_names <- function(df) {
  # Reproduce the column-name style used in the earlier scripts.
  #
  # make.names() converts spaces and punctuation to syntactically valid names.
  # unique = TRUE prevents duplicate names after cleaning.
  names(df) <- make.names(names(df), unique = TRUE)
  df
}

save_table <- function(x, filename, subdir = "outputs/tables") {
  write.csv(
    x,
    file = file.path(subdir, filename),
    row.names = FALSE
  )
}

find_required_col <- function(df, candidates, context) {
  # Locate a required column using a set of possible cleaned names.
  #
  # The clinical/outcome exports may change slightly between REDCap versions.
  # This helper allows the script to tolerate minor differences in punctuation
  # after make.names() has standardised the column names.
  hit <- intersect(candidates, names(df))[1]

  if (is.na(hit)) {
    stop(
      paste0(
        "Could not find required column for ", context, ".\n",
        "Tried: ", paste(candidates, collapse = ", "), "\n",
        "Available columns include: ", paste(head(names(df), 30), collapse = ", ")
      )
    )
  }

  hit
}

has_med <- function(x, pattern) {
  # Convert medication free-text/multi-select fields into binary indicators.
  # Missing medication fields remain missing, rather than being assumed absent.
  case_when(
    is.na(x) | x == "" ~ NA_real_,
    str_detect(x, regex(pattern, ignore_case = TRUE)) ~ 1,
    TRUE ~ 0
  )
}

parse_numeric_safely <- function(x) {
  # Convert values to numeric while tolerating common pathology notation.
  #
  # Pathology files can include values such as "<2.0". parse_number() extracts
  # the numeric component. Non-result text such as "clotted" or "no result" is
  # treated as missing. Warnings are suppressed because these values are expected
  # laboratory non-results rather than coding errors.
  if (is.numeric(x)) {
    return(x)
  }

  suppressWarnings(
    readr::parse_number(
      as.character(x),
      na = c(
        "",
        "NA",
        "NaN",
        "N/A",
        "clotted",
        "Clotted",
        "no result",
        "No result",
        "not done",
        "Not done"
      )
    )
  )
}

collapse_replicates <- function(df, id_col = "Study.ID") {
  # Collapse technical or duplicate rows to one row per patient.
  #
  # For numeric assay variables, duplicate rows are averaged using available
  # values. This is preferable to arbitrarily keeping the first row when the
  # duplicate rows appear to be technical repeats or batch re-runs.
  #
  # For non-numeric metadata variables, the first non-missing value is retained.
  # Most metadata variables are removed before modelling, but keeping them in
  # audit objects can help diagnose unexpected behaviour.
  df %>%
    group_by(.data[[id_col]]) %>%
    summarise(
      across(
        everything(),
        ~ if (is.numeric(.x)) {
          if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
        } else {
          first_nonmissing(.x)
        }
      ),
      .groups = "drop"
    )
}

# Read study-ID reference sheet
#
# The workbook includes a sheet listing expected project IDs. This is used as an
# audit and optional safety filter. It helps prevent technical sample IDs or
# unexpected non-project rows entering the final dataset.
study_id_reference <- read_excel(
  omics_workbook,
  sheet = "Study IDs for Alec Project",
  col_names = FALSE
) %>%
  rename(raw_study_id = 1) %>%
  transmute(Study.ID = clean_id(raw_study_id)) %>%
  filter(!is.na(Study.ID), Study.ID != "") %>%
  distinct()

# Import clinical baseline data

clinical_raw <- read.csv(
  clinical_file,
  fileEncoding = "latin1",
  check.names = FALSE
) %>%
  clean_column_names()

clinical_id_col <- find_required_col(
  clinical_raw,
  candidates = c("Study.ID", "Study.ID."),
  context = "clinical patient ID"
)

event_col <- find_required_col(
  clinical_raw,
  candidates = c("Event.Name", "Event.Name."),
  context = "clinical event name"
)

repeat_col <- find_required_col(
  clinical_raw,
  candidates = c("Repeat.Instrument", "Repeat.Instrument."),
  context = "clinical repeat instrument"
)

clinical_source <- clinical_raw %>%
  rename(Study.ID = all_of(clinical_id_col)) %>%
  mutate(Study.ID = clean_id(Study.ID))

baseline_source <- clinical_source %>%
  filter(
    .data[[event_col]] == "Enrolment",
    is.na(.data[[repeat_col]]) | .data[[repeat_col]] == ""
  )

# Helper for clinical columns after make.names()
get_clin <- function(name) {
  find_required_col(baseline_source, name, context = paste("clinical variable", name[1]))
}

baseline_data <- baseline_source %>%
  arrange(Study.ID) %>%
  group_by(Study.ID) %>%
  summarise(
    # Demographics
    age = first_nonmissing(.data[[get_clin(c("Age.at.Visit"))]]),
    gender = first_nonmissing(.data[[get_clin(c("Gender"))]]),
    ethnicity = first_nonmissing(.data[[get_clin(c("Ethnicity"))]]),
    bmi = first_nonmissing(.data[[get_clin(c("BMI"))]]),
    smoking = first_nonmissing(.data[[get_clin(c("Smoking"))]]),

    # Alcohol exposure
    alcohol_frequency = first_nonmissing(.data[[get_clin(c(
      "How.often.has.the.participant.had.a.drink.containing.alcohol.in.the.past.year."
    ))]]),
    alcohol_drinks = first_nonmissing(.data[[get_clin(c(
      "How.many.drinks.did.the.participant.have.on.a.typical.day.when.drinking.in.the.past.year."
    ))]]),
    alcohol_dependency = first_nonmissing(.data[[get_clin(c(
      "Does.the.participant.have.a.history.of.alcohol.dependency."
    ))]]),

    # Cardiometabolic factors
    masld = first_nonmissing(.data[[get_clin(c("NAFLD.NASH", "MAFLD", "MASLD"))]]),
    diabetes = first_nonmissing(.data[[get_clin(c("Diabetes"))]]),
    hypertension = first_nonmissing(.data[[get_clin(c("Hypertension"))]]),
    dyslipidemia = first_nonmissing(.data[[get_clin(c("Dyslipidemia"))]]),
    mets_count = first_nonmissing(.data[[get_clin(c("Number.of.metabolic.syndrome.risk.factors"))]]),
    family_history_cardiac = first_nonmissing(.data[[get_clin(c(
      "Does.the.participant.have.a.family.history.of.cardiac.disease."
    ))]]),

    # Medication and recent exposure fields
    medication_use = first_nonmissing(.data[[get_clin(c("Medication.use"))]]),
    rifaximin_current = first_nonmissing(.data[[get_clin(c("Rifaximin.use..current."))]]),
    antivirals_3mo = first_nonmissing(.data[[get_clin(c(
      "Has.the.participant.used.antivirals.in.the.last.three.months."
    ))]]),
    antibiotics_3mo = first_nonmissing(.data[[get_clin(c(
      "Has.the.participant.used.antibiotics.in.the.last.three.months..excluding.rifaximin.."
    ))]]),
    supplement_use = first_nonmissing(.data[[get_clin(c("Supplement.use"))]]),

    # Liver disease severity and complications
    portal_hypertension = first_nonmissing(.data[[get_clin(c("Portal.hypertension"))]]),
    cirrhosis_stage = first_nonmissing(.data[[get_clin(c("Stage.of.cirrhosis"))]]),
    decompensated_cirrhosis_raw = first_nonmissing(.data[[get_clin(c("Decompensated.cirrhosis"))]]),
    encephalopathy = first_nonmissing(.data[[get_clin(c("History.of.encephalopathy"))]]),
    ascites = first_nonmissing(.data[[get_clin(c("Ascites"))]]),
    thrombocytopenia = first_nonmissing(.data[[get_clin(c("Thrombocytopenia"))]]),
    tips = first_nonmissing(.data[[get_clin(c(
      "TIPS..Transjugular.intrahepatic.portosystemic.shunt."
    ))]]),

    # Primary liver disease aetiology
    # A single mutually-exclusive factor is constructed from the source columns.
    # This can be collinear with MASLD status and should be interpreted carefully.
    primary_aetiology = case_when(
      first_nonmissing(.data[[get_clin(c("Is.MASLD.the.primary.aetiology."))]]) == "Yes" ~ "MASLD",
      first_nonmissing(.data[[get_clin(c("Is.alcoholic.liver.disease.the.primary.aetiology."))]]) == "Yes" ~ "Alcohol",
      first_nonmissing(.data[[get_clin(c("Is.Hepatitis.C.the.primary.aetiology."))]]) == "Yes" ~ "HCV",
      first_nonmissing(.data[[get_clin(c("Is.Hepatitis.B.the.primary.aetiology."))]]) == "Yes" ~ "HBV",
      first_nonmissing(.data[[get_clin(c("Is.autoimmune.hepatitis.the.primary.aetiology."))]]) == "Yes" ~ "AIH",
      first_nonmissing(.data[[get_clin(c("Is.haemochromatosis.the.primary.aetiology."))]]) == "Yes" ~ "Haemochromatosis",
      first_nonmissing(.data[[get_clin(c("Is.primary.biliary.cholangitis.the.primary.aetiology."))]]) == "Yes" ~ "PBC",
      TRUE ~ "Other/Unknown"
    ),

    .groups = "drop"
  ) %>%
  mutate(
    # Medication indicators
    med_metformin = has_med(medication_use, "Metformin"),
    med_statin = has_med(medication_use, "Statin"),
    med_aspirin = has_med(medication_use, "Aspirin"),
    med_ppi = has_med(medication_use, "Proton-pump inhibitors"),
    med_glp1_gip = has_med(medication_use, "GLP-1|GIP"),
    med_none = if_else(medication_use == "None", 1, 0),

    # Binary simplifications
    smoking_bin = case_when(
      smoking == "Non-smoker" ~ 0,
      smoking %in% c("Ex-smoker", "Current smoker") ~ 1,
      TRUE ~ NA_real_
    ),
    diabetes_bin = case_when(
      diabetes == "No" ~ 0,
      diabetes %in% c("Type 1", "Type 2") ~ 1,
      TRUE ~ NA_real_
    ),

    # Numeric conversions
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    mets_count = as.numeric(mets_count),

    # Decompensation indicators from the multi-select raw field
    decomp_ascites = if_else(
      !is.na(decompensated_cirrhosis_raw) & str_detect(decompensated_cirrhosis_raw, "Ascites"),
      1, 0
    ),
    decomp_encephalopathy = if_else(
      !is.na(decompensated_cirrhosis_raw) & str_detect(decompensated_cirrhosis_raw, "Encephalopathy"),
      1, 0
    ),
    decomp_variceal_bleed = if_else(
      !is.na(decompensated_cirrhosis_raw) & str_detect(decompensated_cirrhosis_raw, "Variceal"),
      1, 0
    ),
    decomp_child_pugh7 = if_else(
      !is.na(decompensated_cirrhosis_raw) & str_detect(decompensated_cirrhosis_raw, "Child-Pugh"),
      1, 0
    ),

    # Modelling factors
    decomp_ascites = factor(decomp_ascites),
    decomp_encephalopathy = factor(decomp_encephalopathy),
    decomp_variceal_bleed = factor(decomp_variceal_bleed),
    decomp_child_pugh7 = factor(decomp_child_pugh7),
    primary_aetiology = factor(primary_aetiology),
    masld = relevel(factor(masld), ref = "No"),
    gender = relevel(factor(gender), ref = "Female"),
    diabetes_bin = relevel(factor(diabetes_bin), ref = "0"),
    smoking_bin = relevel(factor(smoking_bin), ref = "0")
  )

# Import and collapse outcome data

outcome_raw <- read.csv(
  outcome_file,
  fileEncoding = "latin1",
  check.names = FALSE
) %>%
  clean_column_names()

outcome_id_col <- find_required_col(
  outcome_raw,
  candidates = c("Study.ID", "Study.ID."),
  context = "outcome patient ID"
)

mace_hf_col <- find_required_col(
  outcome_raw,
  candidates = c("MACE.5..HF.", "MACE.5.HF", "MACE.5..HF.."),
  context = "MACE heart failure outcome"
)

mace_revasc_col <- find_required_col(
  outcome_raw,
  candidates = c("MACE.5..revasc.", "MACE.5.revasc", "MACE.5..revasc.."),
  context = "MACE revascularisation outcome"
)

any_cvd_col <- find_required_col(
  outcome_raw,
  candidates = c("Any.CVD", "Any.CVD."),
  context = "Any CVD outcome"
)

outcome_patient <- outcome_raw %>%
  rename(Study.ID = all_of(outcome_id_col)) %>%
  mutate(Study.ID = clean_id(Study.ID)) %>%
  group_by(Study.ID) %>%
  summarise(
    mace = as.numeric(any(
      parse_numeric_safely(.data[[mace_hf_col]]) > 0 |
        parse_numeric_safely(.data[[mace_revasc_col]]) > 0,
      na.rm = TRUE
    )),
    cvd = as.numeric(any(
      parse_numeric_safely(.data[[any_cvd_col]]) > 0,
      na.rm = TRUE
    )),
    .groups = "drop"
  )

# Construct clinical model data

baseline_model_data <- baseline_data %>%
  select(
    -decompensated_cirrhosis_raw,
    -medication_use,
    -smoking,
    -diabetes
  )

clinical_model_data <- baseline_model_data %>%
  left_join(outcome_patient, by = "Study.ID") %>%
  filter(!is.na(mace), !is.na(cvd))

save_table(
  clinical_model_data %>%
    summarise(
      n_patients = n(),
      n_mace = sum(mace == 1, na.rm = TRUE),
      n_no_mace = sum(mace == 0, na.rm = TRUE),
      mace_rate = mean(mace == 1, na.rm = TRUE),
      n_cvd = sum(cvd == 1, na.rm = TRUE),
      cvd_rate = mean(cvd == 1, na.rm = TRUE)
    ),
  "audit_clinical_outcome_summary.csv",
  subdir = "outputs/audit"
)

# Import workbook modalities

modality_specs <- tribble(
  ~modality,          ~sheet,             ~id_col_raw, ~prefix,  ~metadata_cols_raw,
  "lipids",          "Lipids",           "sourceID",  "lipid_", list(c("sourceID", "Timepoint", "Batch", "Label2")),
  "small_molecules", "Small Molecules",  "Study ID",  "sm_",    list(c("Study ID", "Visit #")),
  "lipoproteins",    "Lipoproteins",     "sourceID",  "lipo_",  list(c("sourceID", "Visit #")),
  "bile_acids",      "Bile Acids",       "sourceID",  "ba_",    list(c("sampleID", "plateID", "tubeLabel", "sourceID")),
  "pathology",       "Pathology",        "Study ID",  "path_",  list(c("Study ID", "Timepoint", "HBV Viral Load Date", "HCV Viral Load Date")),
  "cytokines",       "Cytokines",        "sourceID",  "cyto_",  list(c("sampleID", "sourceID", "tubeLabel", "plateID", "Location", "Sample"))
)

read_workbook_modality <- function(sheet, id_col_raw, prefix, metadata_cols_raw, modality) {
  # Read one workbook sheet and return one row per patient.
  #
  # The ID column is explicit. Do not replace this with automatic ID detection:
  # several sheets contain sample-level or technical ID columns that are not the
  # patient identifier required for merging.
  raw_df <- read_excel(
    omics_workbook,
    sheet = sheet,
    .name_repair = "minimal"
  )

  raw_names <- names(raw_df)

  if (!id_col_raw %in% raw_names) {
    stop(
      paste0(
        "Expected ID column '", id_col_raw, "' not found in sheet '", sheet, "'.\n",
        "Available columns include: ", paste(head(raw_names, 25), collapse = ", ")
      )
    )
  }

  # Audit all ID-like columns before cleaning so anomalous columns can be checked.
  id_like_cols <- raw_names[str_detect(raw_names, regex("id|sample|tube|source|study", ignore_case = TRUE))]

  id_audit <- map_dfr(id_like_cols, function(col) {
    tibble(
      modality = modality,
      sheet = sheet,
      candidate_id_column = col,
      used_as_patient_id = col == id_col_raw,
      n_non_missing = sum(!is.na(raw_df[[col]])),
      n_unique_raw = n_distinct(raw_df[[col]], na.rm = TRUE),
      example_values = paste(head(na.omit(as.character(raw_df[[col]])), 5), collapse = "; ")
    )
  })

  cleaned <- raw_df %>%
    mutate(Study.ID = clean_id(.data[[id_col_raw]]))

  duplicate_audit <- cleaned %>%
    count(Study.ID, name = "n_rows") %>%
    filter(!is.na(Study.ID), Study.ID != "", n_rows > 1) %>%
    mutate(
      modality = modality,
      sheet = sheet,
      id_col_used = id_col_raw
    ) %>%
    select(modality, sheet, id_col_used, Study.ID, n_rows)

  # Keep only assay/lab columns. Technical metadata are dropped before modelling.
  #
  # This is deliberately based on explicit sheet-specific metadata columns rather
  # than a broad pattern such as "id". Broad matching is unsafe here because bile
  # acid variables literally end in "acid" and would be incorrectly removed.
  data_cols_raw <- setdiff(raw_names, metadata_cols_raw)

  if (length(data_cols_raw) == 0) {
    stop(
      paste0(
        "No assay/lab columns were retained for sheet '", sheet, "'. ",
        "Check metadata_cols_raw; it may be too broad."
      )
    )
  }

  message(
    "Retained ", length(data_cols_raw), " assay/lab columns for ", modality,
    ": ", paste(head(data_cols_raw, 8), collapse = ", "),
    ifelse(length(data_cols_raw) > 8, ", ...", "")
  )

  cleaned_data <- cleaned %>%
    select(Study.ID, all_of(data_cols_raw))

  # Clean assay variable names after metadata removal.
  names(cleaned_data) <- make.names(names(cleaned_data), unique = TRUE)

  # Convert all assay/lab columns to numeric. This is appropriate for omics and
  # pathology sheets; any non-numeric text becomes NA and is captured by
  # missingness audit tables.
  cleaned_data <- cleaned_data %>%
    mutate(
      across(
        -Study.ID,
        parse_numeric_safely
      )
    )

  collapsed <- cleaned_data %>%
    collapse_replicates(id_col = "Study.ID")

  # Prefix retained assay/lab variables. The guard makes failure modes clearer if
  # a sheet unexpectedly contains only Study.ID after metadata removal.
  if (ncol(collapsed) > 1) {
    collapsed <- collapsed %>%
      rename_with(~ paste0(prefix, .x), -Study.ID)
  }

  list(
    data = collapsed,
    id_audit = id_audit,
    duplicate_audit = duplicate_audit
  )
}

modality_imports <- pmap(
  modality_specs,
  function(modality, sheet, id_col_raw, prefix, metadata_cols_raw) {
    read_workbook_modality(
      sheet = sheet,
      id_col_raw = id_col_raw,
      prefix = prefix,
      metadata_cols_raw = metadata_cols_raw,
      modality = modality
    )
  }
)

names(modality_imports) <- modality_specs$modality

omics_data <- map(modality_imports, "data")

id_column_audit <- map_dfr(modality_imports, "id_audit")
duplicate_id_audit <- map_dfr(modality_imports, "duplicate_audit")

save_table(id_column_audit, "audit_workbook_id_columns.csv", subdir = "outputs/audit")
save_table(duplicate_id_audit, "audit_workbook_duplicate_patient_ids.csv", subdir = "outputs/audit")

# Patient overlap audit

id_sets <- c(
  list(
    study_reference = study_id_reference$Study.ID,
    clinical = clinical_model_data$Study.ID
  ),
  map(omics_data, ~ .x$Study.ID)
)

id_overlap_summary <- tibble(
  dataset = names(id_sets),
  n_unique_ids = map_int(id_sets, ~ n_distinct(.x)),
  n_overlap_with_reference = map_int(id_sets, ~ length(intersect(.x, study_id_reference$Study.ID))),
  n_extra_vs_reference = map_int(id_sets, ~ length(setdiff(.x, study_id_reference$Study.ID))),
  n_missing_from_reference = map_int(id_sets, ~ length(setdiff(study_id_reference$Study.ID, .x)))
)

save_table(id_overlap_summary, "audit_id_overlap_summary.csv", subdir = "outputs/audit")

id_overlap_details <- imap_dfr(id_sets, function(ids, dataset_name) {
  tibble(
    dataset = dataset_name,
    id_status = c("extra_vs_reference", "missing_from_dataset"),
    ids = c(
      paste(sort(setdiff(ids, study_id_reference$Study.ID)), collapse = "; "),
      paste(sort(setdiff(study_id_reference$Study.ID, ids)), collapse = "; ")
    )
  )
})

save_table(id_overlap_details, "audit_id_overlap_details.csv", subdir = "outputs/audit")

# Restrict to IDs present in all required datasets
#
# This maintains the previous complete-case-by-modality design. If you later want
# modality-specific models to use larger available subsets, this is the section
# to change. For now, all generated modelling files are based on the complete
# intersection across clinical, outcome, omics, bile acids, and pathology.
common_ids <- Reduce(intersect, id_sets)

clinical_common <- clinical_model_data %>%
  filter(Study.ID %in% common_ids)

modalities_common <- map(
  omics_data,
  ~ .x %>% filter(Study.ID %in% common_ids)
)

combined_data <- reduce(
  modalities_common,
  left_join,
  by = "Study.ID",
  .init = clinical_common
)

combined_data <- combined_data %>%
  select(
    -any_of(c(
      "ba_sampleID",
      "ba_plateID",
      "ba_tubeLabel",
      "ba_sourceID"
    ))
  )

metadata_columns_to_remove <- c(
  "lipid_sourceID",
  "lipid_Timepoint",
  "lipid_Batch",
  "lipid_Label2",
  "sm_Study.ID.1",
  "sm_Visit..",
  "lipo_sourceID",
  "lipo_Visit..",
  "cyto_sampleID",
  "cyto_sourceID",
  "cyto_tubeLabel",
  "cyto_plateID",
  "cyto_Location",
  "cyto_Sample",
  "path_Study.ID.1",
  "path_Timepoint"
)

combined_data <- combined_data %>%
  select(-any_of(metadata_columns_to_remove))

# Final sanity checks and audit summaries

cohort_summary <- combined_data %>%
  summarise(
    n_patients = n(),
    n_mace = sum(mace == 1, na.rm = TRUE),
    n_no_mace = sum(mace == 0, na.rm = TRUE),
    mace_rate = mean(mace == 1, na.rm = TRUE),
    n_cvd = sum(cvd == 1, na.rm = TRUE),
    cvd_rate = mean(cvd == 1, na.rm = TRUE),
    n_duplicate_study_ids = sum(duplicated(Study.ID))
  )

save_table(cohort_summary, "combined_cohort_summary.csv")
save_table(cohort_summary, "audit_final_cohort_summary.csv", subdir = "outputs/audit")

# Missingness by variable

missing_prop <- combined_data %>%
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

save_table(missing_prop, "combined_missingness_by_variable.csv")

removed_missing <- missing_prop %>%
  filter(missing_prop > 0.50) %>%
  arrange(desc(missing_prop))

save_table(removed_missing, "combined_removed_gt50_missing.csv")

# Missingness and feature count by modality

modality_prefixes <- tribble(
  ~modality,              ~prefix,
  "Clinical/outcome",     NA_character_,
  "Lipoproteins",         "lipo_",
  "Lipids",               "lipid_",
  "Cytokines",            "cyto_",
  "Small molecules",      "sm_",
  "Bile acids",           "ba_",
  "Pathology",            "path_"
)

omics_prefixes <- modality_prefixes %>%
  filter(!is.na(prefix)) %>%
  pull(prefix)

clinical_columns <- names(combined_data)[
  !str_detect(
    names(combined_data),
    paste0("^(" , paste(omics_prefixes, collapse = "|"), ")")
  )
]

clinical_columns <- setdiff(clinical_columns, c("Study.ID", "mace", "cvd"))

missing_by_modality <- modality_prefixes %>%
  mutate(
    feature_names = map(prefix, function(pfx) {
      if (is.na(pfx)) {
        clinical_columns
      } else {
        names(combined_data)[str_detect(names(combined_data), paste0("^", pfx))]
      }
    }),
    n_features = map_int(feature_names, length),
    mean_missing = map_dbl(feature_names, function(vars) {
      if (length(vars) == 0) {
        NA_real_
      } else {
        mean(colMeans(is.na(combined_data[, vars, drop = FALSE])))
      }
    }),
    median_missing = map_dbl(feature_names, function(vars) {
      if (length(vars) == 0) {
        NA_real_
      } else {
        median(colMeans(is.na(combined_data[, vars, drop = FALSE])))
      }
    }),
    n_missing_gt_20 = map_int(feature_names, function(vars) {
      if (length(vars) == 0) 0 else sum(colMeans(is.na(combined_data[, vars, drop = FALSE])) > 0.20)
    }),
    n_missing_gt_50 = map_int(feature_names, function(vars) {
      if (length(vars) == 0) 0 else sum(colMeans(is.na(combined_data[, vars, drop = FALSE])) > 0.50)
    }),
    n_missing_gt_80 = map_int(feature_names, function(vars) {
      if (length(vars) == 0) 0 else sum(colMeans(is.na(combined_data[, vars, drop = FALSE])) > 0.80)
    })
  ) %>%
  select(-feature_names)

save_table(missing_by_modality, "combined_missingness_by_modality.csv")

# Apply missingness filter
#
# This keeps the historical >50% missingness threshold. Outcomes and Study.ID are
# retained even if their missingness is zero by construction.
keep_cols <- missing_prop %>%
  filter(missing_prop <= 0.50) %>%
  pull(variable)

combined_filtered <- combined_data %>%
  select(all_of(keep_cols))

# Save final modelling dataset

write.csv(
  combined_filtered,
  ccrvi_config$files$compiled_cohort,
  row.names = FALSE
)

# Optional: save an unfiltered version for audit and EDA.
write.csv(
  combined_data,
  ccrvi_config$files$unfiltered_cohort,
  row.names = FALSE
)

message("Compilation complete.")
message("Final complete-intersection cohort: ", nrow(combined_filtered), " patients.")
message("Filtered dataset saved as combined_multimodal_complete_cohort.csv")
message("Audit outputs saved in outputs/audit")
