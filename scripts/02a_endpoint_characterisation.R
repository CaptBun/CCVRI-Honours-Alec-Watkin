# Repository stage: endpoint
# Source supplied as: 02a_mace_clinical_event_characterisation_v2(5).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# Clinical and event-level characterisation of the MACE endpoint
#
# Project:
#   Cardiovascular Risk Stratification in Liver Cirrhosis: Multiomics integration
#
# Purpose:
#   This script characterises the definitive binary MACE endpoint used in the
#   complete multimodal modelling cohort and examines the cardiovascular events
#   represented by the two clinical-team-defined MACE-5 composites.
#
# Definitive modelling endpoint:
#   mace = MACE-5 (with HF) OR MACE-5 (with revascularisation)
#
# Clinical-team definitions:
#
#   MACE-5 with heart failure:
#     - Myocardial infarction, baseline or follow-up
#     - Cerebrovascular accident or stroke, baseline or follow-up
#     - Decompensated heart failure requiring hospitalisation,
#       baseline or follow-up
#     - Unstable angina, baseline or follow-up
#     - Cardiovascular cause of death
#
#   MACE-5 with revascularisation:
#     - Myocardial infarction, baseline or follow-up
#     - Cerebrovascular accident or stroke, baseline or follow-up
#     - Coronary revascularisation by stent or CABG,
#       baseline or follow-up
#     - Unstable angina, baseline or follow-up
#     - Cardiovascular cause of death
#
# Important interpretation:
#   The compiled cohort variable 'mace' remains the definitive endpoint.
#   Raw cardiovascular fields are used only to characterise recorded events.
#   They are not used to replace or redefine the adjudicated endpoint.
#
#   Failure to reconstruct every MACE case from raw questionnaire fields does
#   not invalidate the endpoint. It may indicate that adjudication incorporated
#   information not completely represented in the exported item-level fields.
#
# Main outputs:
#   1. Cohort and linkage audit
#   2. Formal endpoint-definition table
#   3. MACE-5 composite frequency and overlap
#   4. Recorded cardiovascular-event frequencies among MACE cases
#   5. Patient-level event matrix
#   6. Event co-occurrence table and heatmap
#   7. Raw-event coverage of the adjudicated MACE endpoint
#   8. MACE versus non-MACE clinical description
#   9. Medication-specific review
#  10. Variable audit for pre-modelling decisions

# Load packages
library(tidyverse)

set.seed(123)

# Create output folders
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/models", showWarnings = FALSE, recursive = TRUE)

# File paths
compiled_cohort_file <- ccrvi_config$files$compiled_cohort
raw_clinical_file <- ccrvi_config$files$clinical

# Load compiled modelling cohort
model_source <- read.csv(
  compiled_cohort_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Load raw clinical data
clinical_source <- read.csv(
  raw_clinical_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "Windows-1252"
)

# Remove unnamed columns from the raw clinical export
valid_clinical_names <- !is.na(names(clinical_source)) &
  str_trim(names(clinical_source)) != ""

clinical_source <- clinical_source[
  ,
  valid_clinical_names,
  drop = FALSE
]

names(clinical_source) <- str_trim(names(clinical_source))

# Stop if duplicate named columns remain after trimming
duplicate_clinical_names <- unique(
  names(clinical_source)[duplicated(names(clinical_source))]
)

if (length(duplicate_clinical_names) > 0) {
  stop(
    paste0(
      "Duplicate named columns remain in the raw clinical file: ",
      paste(duplicate_clinical_names, collapse = ", ")
    )
  )
}

# Helper: standardise participant identifiers
standardise_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower() %>%
    str_remove("-1$")
}

# Helper: convert common binary codings to 0/1
safe_binary <- function(x) {
  x_chr <- x %>%
    as.character() %>%
    str_trim() %>%
    str_to_lower()

  case_when(
    x_chr %in% c(
      "1",
      "yes",
      "y",
      "true",
      "positive"
    ) ~ 1,
    x_chr %in% c(
      "0",
      "no",
      "n",
      "false",
      "negative"
    ) ~ 0,
    TRUE ~ suppressWarnings(
      as.numeric(x_chr)
    )
  )
}

# Helper: combine binary fields across rows or visits
#
# Returns 1 if any source field is positive, 0 if at least one source field is
# observed and all observed fields are negative, and NA if every source field is
# missing.
combine_binary_fields <- function(...) {
  values <- c(...)

  if (all(is.na(values))) {
    return(NA_real_)
  }

  if (any(values == 1, na.rm = TRUE)) {
    return(1)
  }

  0
}

# Identify participant-ID column in the compiled cohort
possible_model_id_columns <- c(
  "Study.ID",
  "Study ID",
  "study_id",
  "StudyID"
)

model_id_column <- possible_model_id_columns[
  possible_model_id_columns %in% names(model_source)
][1]

if (is.na(model_id_column)) {
  stop(
    paste0(
      "No recognised participant-ID column was found in the compiled cohort. ",
      "Expected one of: ",
      paste(possible_model_id_columns, collapse = ", ")
    )
  )
}

# Confirm definitive endpoint
if (!"mace" %in% names(model_source)) {
  stop(
    "The compiled modelling cohort does not contain a column named 'mace'."
  )
}

# Standardise the modelling cohort
model_source <- model_source %>%
  mutate(
    study_id = standardise_id(
      .data[[model_id_column]]
    ),
    mace = safe_binary(mace)
  ) %>%
  filter(
    !is.na(study_id),
    study_id != "",
    mace %in% c(0, 1)
  )

if (anyDuplicated(model_source$study_id) > 0) {
  warning(
    "Duplicate participant IDs were found in the compiled cohort. ",
    "Only the first row per participant is retained."
  )

  model_source <- model_source %>%
    distinct(
      study_id,
      .keep_all = TRUE
    )
}

# Define exact clinical columns
#
# These names are explicit to prevent unrelated historical cardiovascular
# variables from being selected by broad regular-expression matching.
event_column_map <- list(
  myocardial_infarction = c(
    "Have they had a myocardial infarction?",
    "Have they had a myocardial infarction? .1"
  ),
  stroke = c(
    "Cerebrovascular Accident (either hemorrhagic or ischaemic)?",
    paste0(
      "Has the participant experienced a cerebrovascular accident ",
      "(either hemorrhagic or ischemic) since the last study visit?"
    )
  ),
  decompensated_heart_failure = c(
    paste0(
      "Have they had decompensated heart failure ",
      "(requiring hospital admission)?"
    ),
    paste0(
      "Have they had decompensated heart failure ",
      "(requiring hospital admission)? .1"
    )
  ),
  unstable_angina = c(
    "Have they had any episodes of unstable angina?",
    "Have they had any episodes of unstable angina? .1"
  ),
  coronary_stent = c(
    "Do they have a coronary artery stent?",
    "Do they have a coronary artery stent?.1"
  ),
  coronary_bypass = c(
    "Have they had a coronary artery bypass?",
    "Have they had a coronary artery bypass?.1"
  ),
  cardiovascular_death = c(
    "Cause of death"
  )
)

required_composite_columns <- c(
  "Study ID",
  "MACE-5 (HF)",
  "MACE-5 (revasc)"
)

missing_composite_columns <- setdiff(
  required_composite_columns,
  names(clinical_source)
)

if (length(missing_composite_columns) > 0) {
  stop(
    paste0(
      "Required endpoint columns are missing from the raw clinical file: ",
      paste(missing_composite_columns, collapse = ", ")
    )
  )
}

# Audit event columns actually available
event_column_audit <- imap_dfr(
  event_column_map,
  function(column_names, event_name) {
    tibble(
      event = event_name,
      requested_column = column_names,
      column_present = column_names %in%
        names(clinical_source)
    )
  }
)

write.csv(
  event_column_audit,
  "outputs/tables/mace_event_column_audit.csv",
  row.names = FALSE
)

# Retrieve a column safely
#
# A missing source column returns an all-NA vector so the script remains usable
# while documenting incomplete raw-field availability.
get_clinical_column <- function(data, column_name) {
  if (!column_name %in% names(data)) {
    return(
      rep(
        NA_character_,
        nrow(data)
      )
    )
  }

  data[[column_name]]
}

# Prepare raw visit-level cardiovascular events
clinical_events_visit <- clinical_source %>%
  transmute(
    study_id = standardise_id(`Study ID`),

    mace5_hf_count = suppressWarnings(
      as.numeric(
        as.character(`MACE-5 (HF)`)
      )
    ),

    mace5_revascularisation_count = suppressWarnings(
      as.numeric(
        as.character(`MACE-5 (revasc)`)
      )
    ),

    mi_baseline = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$myocardial_infarction[1]
      )
    ),

    mi_followup = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$myocardial_infarction[2]
      )
    ),

    stroke_baseline = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$stroke[1]
      )
    ),

    stroke_followup = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$stroke[2]
      )
    ),

    decomp_hf_baseline = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$decompensated_heart_failure[1]
      )
    ),

    decomp_hf_followup = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$decompensated_heart_failure[2]
      )
    ),

    unstable_angina_baseline = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$unstable_angina[1]
      )
    ),

    unstable_angina_followup = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$unstable_angina[2]
      )
    ),

    stent_baseline = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$coronary_stent[1]
      )
    ),

    stent_followup = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$coronary_stent[2]
      )
    ),

    cabg_baseline = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$coronary_bypass[1]
      )
    ),

    cabg_followup = safe_binary(
      get_clinical_column(
        clinical_source,
        event_column_map$coronary_bypass[2]
      )
    ),

    cause_of_death = as.character(
      get_clinical_column(
        clinical_source,
        event_column_map$cardiovascular_death[1]
      )
    )
  ) %>%
  filter(
    !is.na(study_id),
    study_id != ""
  ) %>%
  mutate(
    cardiovascular_death = case_when(
      is.na(cause_of_death) |
        str_trim(cause_of_death) == "" ~
        NA_real_,
      str_detect(
        cause_of_death,
        regex(
          "cardio|cardiac|vascular|myocard|heart",
          ignore_case = TRUE
        )
      ) ~ 1,
      TRUE ~ 0
    )
  )

# Collapse repeated visits to one participant-level event record
#
# A positive event at any recorded visit classifies the participant as positive.
# Event counts use the maximum value recorded across visits.
clinical_events <- clinical_events_visit %>%
  group_by(study_id) %>%
  summarise(
    mace5_hf_count = if (
      all(is.na(mace5_hf_count))
    ) {
      NA_real_
    } else {
      max(
        mace5_hf_count,
        na.rm = TRUE
      )
    },

    mace5_revascularisation_count = if (
      all(is.na(mace5_revascularisation_count))
    ) {
      NA_real_
    } else {
      max(
        mace5_revascularisation_count,
        na.rm = TRUE
      )
    },

    myocardial_infarction = combine_binary_fields(
      mi_baseline,
      mi_followup
    ),

    stroke = combine_binary_fields(
      stroke_baseline,
      stroke_followup
    ),

    decompensated_heart_failure =
      combine_binary_fields(
        decomp_hf_baseline,
        decomp_hf_followup
      ),

    unstable_angina = combine_binary_fields(
      unstable_angina_baseline,
      unstable_angina_followup
    ),

    coronary_stent = combine_binary_fields(
      stent_baseline,
      stent_followup
    ),

    coronary_bypass = combine_binary_fields(
      cabg_baseline,
      cabg_followup
    ),

    cardiovascular_death = combine_binary_fields(
      cardiovascular_death
    ),

    .groups = "drop"
  ) %>%
  mutate(
    mace5_hf = case_when(
      is.na(mace5_hf_count) ~ NA_real_,
      mace5_hf_count > 0 ~ 1,
      TRUE ~ 0
    ),

    mace5_revascularisation = case_when(
      is.na(mace5_revascularisation_count) ~
        NA_real_,
      mace5_revascularisation_count > 0 ~ 1,
      TRUE ~ 0
    ),

    coronary_revascularisation = case_when(
      coronary_stent == 1 |
        coronary_bypass == 1 ~ 1,
      coronary_stent == 0 &
        coronary_bypass == 0 ~ 0,
      TRUE ~ NA_real_
    ),

    any_recorded_qualifying_event = case_when(
      myocardial_infarction == 1 |
        stroke == 1 |
        decompensated_heart_failure == 1 |
        unstable_angina == 1 |
        coronary_revascularisation == 1 |
        cardiovascular_death == 1 ~ 1,

      myocardial_infarction == 0 &
        stroke == 0 &
        decompensated_heart_failure == 0 &
        unstable_angina == 0 &
        coronary_revascularisation == 0 &
        cardiovascular_death == 0 ~ 0,

      TRUE ~ NA_real_
    ),

    composite_pattern = case_when(
      mace5_hf == 1 &
        mace5_revascularisation == 1 ~
        "Both MACE-5 composites",

      mace5_hf == 1 &
        mace5_revascularisation == 0 ~
        "MACE-5 with HF only",

      mace5_hf == 0 &
        mace5_revascularisation == 1 ~
        "MACE-5 with revascularisation only",

      mace5_hf == 0 &
        mace5_revascularisation == 0 ~
        "Neither MACE-5 composite",

      TRUE ~
        "Incomplete composite data"
    )
  )

# Link raw clinical events to the exact complete modelling cohort
endpoint_data <- model_source %>%
  select(
    study_id,
    mace
  ) %>%
  left_join(
    clinical_events,
    by = "study_id"
  )

# Formal endpoint-definition table
endpoint_definition <- tribble(
  ~event, ~mace5_with_hf, ~mace5_with_revascularisation,
  "Myocardial infarction", TRUE, TRUE,
  "Cerebrovascular accident / stroke", TRUE, TRUE,
  "Decompensated heart failure requiring hospitalisation", TRUE, FALSE,
  "Coronary revascularisation: stent or CABG", FALSE, TRUE,
  "Unstable angina", TRUE, TRUE,
  "Cardiovascular cause of death", TRUE, TRUE
)

# Cohort and linkage audit
cohort_audit <- tibble(
  n_complete_cohort = nrow(model_source),
  n_mace = sum(
    model_source$mace == 1
  ),
  n_non_mace = sum(
    model_source$mace == 0
  ),
  mace_prevalence = mean(
    model_source$mace == 1
  ),
  n_linked_to_raw_clinical = sum(
    endpoint_data$study_id %in%
      clinical_events$study_id
  ),
  n_unlinked = sum(
    !endpoint_data$study_id %in%
      clinical_events$study_id
  ),
  n_mace_with_composite_data = sum(
    endpoint_data$mace == 1 &
      (
        !is.na(endpoint_data$mace5_hf) |
          !is.na(
            endpoint_data$mace5_revascularisation
          )
      )
  )
)

# MACE-5 composite frequency among definitive MACE cases
composite_frequency <- tibble(
  composite = c(
    "MACE-5 with heart failure",
    "MACE-5 with revascularisation"
  ),

  n_positive_full_cohort = c(
    sum(
      endpoint_data$mace5_hf == 1,
      na.rm = TRUE
    ),
    sum(
      endpoint_data$mace5_revascularisation == 1,
      na.rm = TRUE
    )
  ),

  n_positive_among_mace = c(
    sum(
      endpoint_data$mace == 1 &
        endpoint_data$mace5_hf == 1,
      na.rm = TRUE
    ),
    sum(
      endpoint_data$mace == 1 &
        endpoint_data$mace5_revascularisation == 1,
      na.rm = TRUE
    )
  )
) %>%
  mutate(
    proportion_of_mace =
      n_positive_among_mace /
      sum(endpoint_data$mace == 1)
  )

# Mutually exclusive overlap of the two MACE-5 composites
composite_overlap <- endpoint_data %>%
  filter(
    mace == 1
  ) %>%
  count(
    composite_pattern,
    name = "n_participants"
  ) %>%
  mutate(
    proportion_of_mace =
      n_participants /
      sum(n_participants)
  ) %>%
  arrange(
    desc(n_participants)
  )

# Event variables used for detailed characterisation
event_vars <- c(
  "myocardial_infarction",
  "stroke",
  "decompensated_heart_failure",
  "unstable_angina",
  "coronary_revascularisation",
  "cardiovascular_death"
)

event_labels <- c(
  myocardial_infarction =
    "Myocardial infarction",
  stroke =
    "Cerebrovascular accident / stroke",
  decompensated_heart_failure =
    "Decompensated heart failure",
  unstable_angina =
    "Unstable angina",
  coronary_revascularisation =
    "Coronary revascularisation",
  cardiovascular_death =
    "Cardiovascular death"
)

# Event frequency among MACE cases
event_frequency_mace <- map_dfr(
  event_vars,
  function(event_name) {
    
    # Extract the complete event vector for all MACE-positive participants
    values <- endpoint_data[[event_name]][
      endpoint_data$mace == 1
    ]
    
    n_mace_total <- sum(
      endpoint_data$mace == 1
    )
    
    tibble(
      event = event_name,
      event_label = unname(
        event_labels[event_name]
      ),
      n_mace_participants = n_mace_total,
      n_event_non_missing = sum(
        !is.na(values)
      ),
      n_event_positive = sum(
        values == 1,
        na.rm = TRUE
      ),
      n_event_negative = sum(
        values == 0,
        na.rm = TRUE
      ),
      n_event_missing = sum(
        is.na(values)
      ),
      proportion_of_mace = sum(
        values == 1,
        na.rm = TRUE
      ) / n_mace_total,
      proportion_of_observed = if_else(
        sum(!is.na(values)) > 0,
        sum(
          values == 1,
          na.rm = TRUE
        ) / sum(!is.na(values)),
        NA_real_
      )
    )
  }
) %>%
  arrange(
    desc(n_event_positive)
  )

# Event frequency in both outcome groups
event_frequency_by_outcome <- map_dfr(
  event_vars,
  function(event_name) {
    endpoint_data %>%
      group_by(mace) %>%
      summarise(
        event = event_name,
        event_label = unname(
          event_labels[event_name]
        ),
        n_non_missing = sum(
          !is.na(
            .data[[event_name]]
          )
        ),
        n_event_positive = sum(
          .data[[event_name]] == 1,
          na.rm = TRUE
        ),
        event_prevalence = mean(
          .data[[event_name]] == 1,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
  }
)

# Patient-level event matrix for MACE cases
mace_event_matrix <- endpoint_data %>%
  filter(
    mace == 1
  ) %>%
  select(
    study_id,
    mace,
    mace5_hf,
    mace5_revascularisation,
    composite_pattern,
    all_of(event_vars)
  ) %>%
  mutate(
    n_recorded_qualifying_events = rowSums(
      across(
        all_of(event_vars),
        ~ replace_na(.x, 0)
      )
    ),
    has_any_recorded_qualifying_event =
      n_recorded_qualifying_events > 0
  ) %>%
  arrange(
    desc(n_recorded_qualifying_events),
    study_id
  )

# Coverage of the definitive MACE endpoint by available raw event fields
event_coverage_audit <- mace_event_matrix %>%
  summarise(
    n_mace = n(),
    n_with_at_least_one_recorded_event = sum(
      has_any_recorded_qualifying_event
    ),
    proportion_with_at_least_one_recorded_event =
      mean(has_any_recorded_qualifying_event),
    n_with_no_recorded_event = sum(
      !has_any_recorded_qualifying_event
    ),
    n_with_one_recorded_event = sum(
      n_recorded_qualifying_events == 1
    ),
    n_with_multiple_recorded_events = sum(
      n_recorded_qualifying_events >= 2
    ),
    median_recorded_events = median(
      n_recorded_qualifying_events
    ),
    maximum_recorded_events = max(
      n_recorded_qualifying_events
    )
  )

# Event-count distribution among MACE cases
event_count_distribution <- mace_event_matrix %>%
  count(
    n_recorded_qualifying_events,
    name = "n_mace_participants"
  ) %>%
  mutate(
    proportion_of_mace =
      n_mace_participants /
      sum(n_mace_participants)
  )

# Event co-occurrence matrix among MACE cases
#
# Missing values are treated as zero only for this descriptive co-occurrence
# table. The missingness table should be reviewed alongside the heatmap.
mace_event_binary <- mace_event_matrix %>%
  select(
    all_of(event_vars)
  ) %>%
  mutate(
    across(
      everything(),
      ~ replace_na(.x, 0)
    )
  ) %>%
  as.matrix()

event_cooccurrence_matrix <-
  t(mace_event_binary) %*%
  mace_event_binary

event_cooccurrence_table <-
  as.data.frame(
    event_cooccurrence_matrix
  ) %>%
  rownames_to_column(
    "event_1"
  ) %>%
  pivot_longer(
    cols = -event_1,
    names_to = "event_2",
    values_to = "n_cooccurring"
  ) %>%
  mutate(
    event_1_label = unname(
      event_labels[event_1]
    ),
    event_2_label = unname(
      event_labels[event_2]
    )
  )

# Distinct combinations of recorded events among MACE cases
event_combination_table <- mace_event_matrix %>%
  rowwise() %>%
  mutate(
    recorded_event_combination = {
      positive_events <- event_labels[
        c_across(
          all_of(event_vars)
        ) == 1
      ]

      if (length(positive_events) == 0) {
        "No qualifying raw event recorded"
      } else {
        paste(
          positive_events,
          collapse = " + "
        )
      }
    }
  ) %>%
  ungroup() %>%
  count(
    recorded_event_combination,
    name = "n_mace_participants"
  ) %>%
  mutate(
    proportion_of_mace =
      n_mace_participants /
      sum(n_mace_participants)
  ) %>%
  arrange(
    desc(n_mace_participants)
  )

# Clinical phenotype groups based on the two MACE-5 composites
#
# These are descriptive endpoint phenotypes rather than new outcomes.
mace_phenotype <- endpoint_data %>%
  filter(
    mace == 1
  ) %>%
  transmute(
    study_id,
    phenotype = factor(
      composite_pattern,
      levels = c(
        "MACE-5 with HF only",
        "MACE-5 with revascularisation only",
        "Both MACE-5 composites",
        "Incomplete composite data",
        "Neither MACE-5 composite"
      )
    )
  )

# Define clinical variable groups for MACE versus non-MACE description
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

pathology_vars <- names(
  model_source
) %>%
  str_subset(
    "^path_"
  )

clinical_review_vars <- unique(
  c(
    demographic_vars,
    lifestyle_vars,
    cardiometabolic_vars,
    liver_severity_vars,
    medication_vars,
    pathology_vars
  )
)

present_review_vars <- intersect(
  clinical_review_vars,
  names(model_source)
)

missing_review_vars <- setdiff(
  clinical_review_vars,
  names(model_source)
)

variable_group_lookup <- bind_rows(
  tibble(
    variable = demographic_vars,
    variable_group = "Demographic"
  ),
  tibble(
    variable = lifestyle_vars,
    variable_group = "Lifestyle"
  ),
  tibble(
    variable = cardiometabolic_vars,
    variable_group = "Cardiometabolic"
  ),
  tibble(
    variable = liver_severity_vars,
    variable_group = "Liver disease"
  ),
  tibble(
    variable = medication_vars,
    variable_group = "Medication"
  ),
  tibble(
    variable = pathology_vars,
    variable_group = "Pathology"
  )
) %>%
  distinct(
    variable,
    .keep_all = TRUE
  ) %>%
  filter(
    variable %in%
      present_review_vars
  )

# Identify categorical and continuous variables
continuous_override <- intersect(
  c(
    "age",
    "bmi",
    "alcohol_drinks",
    pathology_vars
  ),
  present_review_vars
)

categorical_vars <- present_review_vars[
  map_lgl(
    model_source[
      present_review_vars
    ],
    function(x) {
      !is.numeric(x) ||
        n_distinct(
          x,
          na.rm = TRUE
        ) <= 5
    }
  )
]

categorical_vars <- setdiff(
  categorical_vars,
  continuous_override
)

continuous_vars <- setdiff(
  present_review_vars,
  categorical_vars
)

# Continuous-variable descriptive table
continuous_table <- map_dfr(
  continuous_vars,
  function(variable) {
    values <- suppressWarnings(
      as.numeric(
        model_source[[variable]]
      )
    )

    tibble(
      mace = model_source$mace,
      value = values
    ) %>%
      group_by(mace) %>%
      summarise(
        variable = variable,
        n_non_missing = sum(
          !is.na(value)
        ),
        n_missing = sum(
          is.na(value)
        ),
        mean = mean(
          value,
          na.rm = TRUE
        ),
        sd = sd(
          value,
          na.rm = TRUE
        ),
        median = median(
          value,
          na.rm = TRUE
        ),
        q1 = quantile(
          value,
          probs = 0.25,
          na.rm = TRUE,
          names = FALSE
        ),
        q3 = quantile(
          value,
          probs = 0.75,
          na.rm = TRUE,
          names = FALSE
        ),
        .groups = "drop"
      )
  }
) %>%
  left_join(
    variable_group_lookup,
    by = "variable"
  ) %>%
  relocate(
    variable_group,
    variable,
    mace
  )

# Categorical-variable level and event-rate table
categorical_event_rates <- map_dfr(
  categorical_vars,
  function(variable) {
    model_source %>%
      transmute(
        mace,
        level = as.character(
          .data[[variable]]
        )
      ) %>%
      mutate(
        level = case_when(
          is.na(level) |
            str_trim(level) == "" ~
            "Missing",
          TRUE ~ level
        )
      ) %>%
      group_by(level) %>%
      summarise(
        variable = variable,
        level_total = n(),
        n_mace = sum(
          mace == 1
        ),
        n_non_mace = sum(
          mace == 0
        ),
        mace_rate = mean(
          mace == 1
        ),
        .groups = "drop"
      )
  }
) %>%
  left_join(
    variable_group_lookup,
    by = "variable"
  ) %>%
  relocate(
    variable_group,
    variable,
    level
  ) %>%
  arrange(
    variable_group,
    variable,
    desc(level_total)
  )

# Medication review
present_medication_vars <- intersect(
  medication_vars,
  names(model_source)
)

medication_by_outcome <- map_dfr(
  present_medication_vars,
  function(variable) {
    exposure <- safe_binary(
      model_source[[variable]]
    )

    tibble(
      medication = variable,
      mace = model_source$mace,
      exposure = exposure
    ) %>%
      group_by(
        medication,
        mace
      ) %>%
      summarise(
        n_non_missing = sum(
          !is.na(exposure)
        ),
        n_exposed = sum(
          exposure == 1,
          na.rm = TRUE
        ),
        exposure_prevalence = if_else(
          n_non_missing > 0,
          n_exposed / n_non_missing,
          NA_real_
        ),
        .groups = "drop"
      )
  }
)

medication_event_rates <- map_dfr(
  present_medication_vars,
  function(variable) {
    exposure <- safe_binary(
      model_source[[variable]]
    )

    tibble(
      medication = variable,
      exposure = exposure,
      mace = model_source$mace
    ) %>%
      filter(
        exposure %in% c(0, 1)
      ) %>%
      group_by(
        medication,
        exposure
      ) %>%
      summarise(
        n_participants = n(),
        n_mace = sum(
          mace == 1
        ),
        mace_rate = mean(
          mace == 1
        ),
        .groups = "drop"
      )
  }
)

# Missingness by MACE status
missingness_by_mace <- model_source %>%
  select(
    mace,
    all_of(
      present_review_vars
    )
  ) %>%
  mutate(
    across(
      -mace,
      ~ is.na(.x) |
        str_trim(
          as.character(.x)
        ) == ""
    )
  ) %>%
  pivot_longer(
    cols = -mace,
    names_to = "variable",
    values_to = "is_missing"
  ) %>%
  group_by(
    variable,
    mace
  ) %>%
  summarise(
    n = n(),
    n_missing = sum(
      is_missing
    ),
    missing_proportion = mean(
      is_missing
    ),
    .groups = "drop"
  ) %>%
  left_join(
    variable_group_lookup,
    by = "variable"
  ) %>%
  relocate(
    variable_group,
    variable,
    mace
  )

# Variable audit for organised pre-modelling review
variable_audit <- map_dfr(
  present_review_vars,
  function(variable) {
    x <- model_source[[variable]]

    missing_indicator <- is.na(x) |
      str_trim(
        as.character(x)
      ) == ""

    observed <- x[
      !missing_indicator
    ]

    level_counts <- table(
      as.character(observed)
    )

    n_unique <- n_distinct(
      observed
    )

    smallest_level_n <- if (
      length(level_counts) == 0
    ) {
      NA_integer_
    } else {
      min(level_counts)
    }

    largest_level_proportion <- if (
      length(level_counts) == 0
    ) {
      NA_real_
    } else {
      max(level_counts) /
        sum(level_counts)
    }

    tibble(
      variable = variable,
      inferred_role = if_else(
        variable %in% categorical_vars,
        "Categorical",
        "Continuous"
      ),
      n_non_missing = length(
        observed
      ),
      n_missing = sum(
        missing_indicator
      ),
      missing_proportion = mean(
        missing_indicator
      ),
      n_unique = n_unique,
      smallest_level_n =
        smallest_level_n,
      largest_level_proportion =
        largest_level_proportion,
      modelling_note = case_when(
        length(observed) == 0 ~
          "No observed data",

        mean(missing_indicator) >= 0.40 ~
          "High missingness: review before modelling",

        variable %in% categorical_vars &
          !is.na(smallest_level_n) &
          smallest_level_n < 5 ~
          paste0(
            "Sparse category: consider clinically ",
            "justified collapsing"
          ),

        variable %in% categorical_vars &
          !is.na(
            largest_level_proportion
          ) &
          largest_level_proportion >= 0.95 ~
          "Near-zero variation",

        TRUE ~
          "No immediate structural concern"
      )
    )
  }
) %>%
  left_join(
    variable_group_lookup,
    by = "variable"
  ) %>%
  relocate(
    variable_group,
    variable
  ) %>%
  arrange(
    variable_group,
    variable
  )

# Figure: overall MACE prevalence
mace_prevalence_plot <- model_source %>%
  mutate(
    outcome_group = factor(
      mace,
      levels = c(0, 1),
      labels = c(
        "Non-MACE",
        "MACE"
      )
    )
  ) %>%
  count(
    outcome_group,
    name = "n"
  ) %>%
  mutate(
    proportion = n /
      sum(n),
    label = paste0(
      n,
      " (",
      scales::percent(
        proportion,
        accuracy = 0.1
      ),
      ")"
    )
  ) %>%
  ggplot(
    aes(
      x = outcome_group,
      y = n
    )
  ) +
  geom_col() +
  geom_text(
    aes(label = label),
    vjust = -0.4
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.12)
    )
  ) +
  labs(
    title =
      "MACE prevalence in the complete modelling cohort",
    x = NULL,
    y = "Number of participants"
  ) +
  theme_minimal(
    base_size = 13
  )

ggsave(
  "outputs/figures/mace_cohort_prevalence.png",
  mace_prevalence_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# Figure: overlap of the two MACE-5 composites
composite_overlap_plot <- composite_overlap %>%
  ggplot(
    aes(
      x = reorder(
        composite_pattern,
        n_participants
      ),
      y = n_participants
    )
  ) +
  geom_col() +
  geom_text(
    aes(
      label = paste0(
        n_participants,
        " (",
        scales::percent(
          proportion_of_mace,
          accuracy = 0.1
        ),
        ")"
      )
    ),
    hjust = -0.1
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.25)
    )
  ) +
  labs(
    title = "Overlap of the two MACE-5 composites",
    subtitle =
      "Mutually exclusive groups among definitive MACE cases",
    x = NULL,
    y = "Number of participants"
  ) +
  theme_minimal(
    base_size = 13
  )

ggsave(
  "outputs/figures/mace_composite_overlap.png",
  composite_overlap_plot,
  width = 9,
  height = 5.5,
  dpi = 300
)

# Figure: recorded cardiovascular events among MACE cases
event_frequency_plot <- event_frequency_mace %>%
  ggplot(
    aes(
      x = reorder(
        event_label,
        n_event_positive
      ),
      y = n_event_positive
    )
  ) +
  geom_col() +
  geom_text(
    aes(
      label = paste0(
        n_event_positive,
        " (",
        scales::percent(
          proportion_of_mace,
          accuracy = 0.1
        ),
        ")"
      )
    ),
    hjust = -0.1
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.25)
    )
  ) +
  labs(
    title =
      "Recorded qualifying cardiovascular events among MACE cases",
    subtitle = paste0(
      "Raw event fields in the complete cohort; ",
      "participants may have multiple events"
    ),
    x = NULL,
    y = "MACE-positive participants"
  ) +
  theme_minimal(
    base_size = 13
  )

ggsave(
  "outputs/figures/mace_recorded_event_frequency.png",
  event_frequency_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Figure: event co-occurrence heatmap
event_cooccurrence_plot <- event_cooccurrence_table %>%
  ggplot(
    aes(
      x = event_1_label,
      y = event_2_label,
      fill = n_cooccurring
    )
  ) +
  geom_tile() +
  geom_text(
    aes(
      label = n_cooccurring
    )
  ) +
  scale_fill_viridis_c(
    option = "C",
    direction = -1
  ) +
  labs(
    title =
      "Co-occurrence of recorded cardiovascular events",
    subtitle =
      "Counts among definitive MACE cases in the complete cohort",
    x = NULL,
    y = NULL,
    fill = "Participants"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  "outputs/figures/mace_event_cooccurrence_heatmap.png",
  event_cooccurrence_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# Figure: number of recorded qualifying events per MACE participant
event_count_plot <- event_count_distribution %>%
  ggplot(
    aes(
      x = factor(
        n_recorded_qualifying_events
      ),
      y = n_mace_participants
    )
  ) +
  geom_col() +
  geom_text(
    aes(
      label = n_mace_participants
    ),
    vjust = -0.3
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.12)
    )
  ) +
  labs(
    title =
      "Number of recorded qualifying events per MACE participant",
    x =
      "Number of distinct recorded cardiovascular events",
    y =
      "MACE-positive participants"
  ) +
  theme_minimal(
    base_size = 13
  )

ggsave(
  "outputs/figures/mace_recorded_event_count_distribution.png",
  event_count_plot,
  width = 8,
  height = 5,
  dpi = 300
)

# Figure: medication prevalence by MACE status
if (nrow(medication_by_outcome) > 0) {
  medication_plot <- medication_by_outcome %>%
    mutate(
      outcome_group = factor(
        mace,
        levels = c(0, 1),
        labels = c(
          "Non-MACE",
          "MACE"
        )
      ),
      medication = medication %>%
        str_remove(
          "^med_"
        ) %>%
        str_replace_all(
          "_",
          " "
        ) %>%
        str_to_sentence()
    ) %>%
    ggplot(
      aes(
        x = medication,
        y = exposure_prevalence,
        fill = outcome_group
      )
    ) +
    geom_col(
      position = position_dodge(
        width = 0.8
      )
    ) +
    coord_flip() +
    scale_y_continuous(
      labels =
        scales::percent_format(
          accuracy = 1
        ),
      limits = c(0, 1)
    ) +
    labs(
      title =
        "Medication prevalence by MACE status",
      x = NULL,
      y =
        "Participants prescribed medication",
      fill = NULL
    ) +
    theme_minimal(
      base_size = 13
    )

  ggsave(
    "outputs/figures/mace_medication_prevalence.png",
    medication_plot,
    width = 9,
    height = 6,
    dpi = 300
  )
}

# Save output tables
write.csv(
  endpoint_definition,
  "outputs/tables/mace_endpoint_definition.csv",
  row.names = FALSE
)

write.csv(
  cohort_audit,
  "outputs/tables/mace_cohort_audit.csv",
  row.names = FALSE
)

write.csv(
  composite_frequency,
  "outputs/tables/mace_composite_frequency.csv",
  row.names = FALSE
)

write.csv(
  composite_overlap,
  "outputs/tables/mace_composite_overlap.csv",
  row.names = FALSE
)

write.csv(
  event_frequency_mace,
  "outputs/tables/mace_recorded_event_frequency.csv",
  row.names = FALSE
)

write.csv(
  event_frequency_by_outcome,
  "outputs/tables/mace_event_frequency_by_outcome.csv",
  row.names = FALSE
)

write.csv(
  mace_event_matrix,
  "outputs/tables/mace_patient_event_matrix.csv",
  row.names = FALSE
)

write.csv(
  event_coverage_audit,
  "outputs/tables/mace_raw_event_coverage_audit.csv",
  row.names = FALSE
)

write.csv(
  event_count_distribution,
  "outputs/tables/mace_recorded_event_count_distribution.csv",
  row.names = FALSE
)

write.csv(
  event_cooccurrence_table,
  "outputs/tables/mace_event_cooccurrence.csv",
  row.names = FALSE
)

write.csv(
  event_combination_table,
  "outputs/tables/mace_event_combinations.csv",
  row.names = FALSE
)

write.csv(
  endpoint_data,
  "outputs/tables/mace_linked_endpoint_and_event_data.csv",
  row.names = FALSE
)

write.csv(
  continuous_table,
  "outputs/tables/mace_table1_continuous.csv",
  row.names = FALSE
)

write.csv(
  categorical_event_rates,
  "outputs/tables/mace_categorical_event_rates.csv",
  row.names = FALSE
)

write.csv(
  medication_by_outcome,
  "outputs/tables/mace_medication_prevalence_by_outcome.csv",
  row.names = FALSE
)

write.csv(
  medication_event_rates,
  "outputs/tables/mace_medication_event_rates.csv",
  row.names = FALSE
)

write.csv(
  missingness_by_mace,
  "outputs/tables/mace_missingness_by_outcome.csv",
  row.names = FALSE
)

write.csv(
  variable_audit,
  "outputs/tables/mace_variable_audit.csv",
  row.names = FALSE
)

write.csv(
  tibble(
    variable = missing_review_vars
  ),
  "outputs/tables/mace_expected_variables_not_found.csv",
  row.names = FALSE
)

# Save all analysis objects
mace_characterisation_outputs <- list(
  endpoint_definition =
    endpoint_definition,
  event_column_audit =
    event_column_audit,
  cohort_audit =
    cohort_audit,
  composite_frequency =
    composite_frequency,
  composite_overlap =
    composite_overlap,
  event_frequency_mace =
    event_frequency_mace,
  event_frequency_by_outcome =
    event_frequency_by_outcome,
  mace_event_matrix =
    mace_event_matrix,
  event_coverage_audit =
    event_coverage_audit,
  event_count_distribution =
    event_count_distribution,
  event_cooccurrence_table =
    event_cooccurrence_table,
  event_combination_table =
    event_combination_table,
  endpoint_data =
    endpoint_data,
  mace_phenotype =
    mace_phenotype,
  continuous_table =
    continuous_table,
  categorical_event_rates =
    categorical_event_rates,
  medication_by_outcome =
    medication_by_outcome,
  medication_event_rates =
    medication_event_rates,
  missingness_by_mace =
    missingness_by_mace,
  variable_audit =
    variable_audit
)

saveRDS(
  mace_characterisation_outputs,
  "outputs/models/mace_clinical_event_characterisation_outputs.rds"
)

# Print the most important results
print(
  endpoint_definition
)

print(
  cohort_audit
)

print(
  composite_overlap
)

print(
  event_frequency_mace
)

print(
  event_coverage_audit
)

message(
  "MACE clinical-event characterisation complete. ",
  "Review the event-column audit and raw-event coverage audit before ",
  "interpreting the detailed event decomposition."
)
