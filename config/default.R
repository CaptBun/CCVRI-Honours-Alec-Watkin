# Define project paths and execution settings
# Paths are interpreted from the repository root. The source data remain local.
# The original seed, folds and analytical algorithms remain in the scripts.
ccrvi_config <- list(
  files = list(
    clinical = "data/raw/clindata.csv",
    outcome = "data/raw/outcomedata.csv",
    workbook = "data/raw/Alec Honours.xlsx",
    compiled_cohort = "data/processed/combined_multimodal_complete_cohort.csv",
    unfiltered_cohort = "data/processed/combined_multimodal_complete_cohort_unfiltered.csv",
    original_results = "outputs/models/aucfilter_multimodality_comparison_2layer_nested_cv_outputs.rds",
    no_aspirin_results = "outputs/models/aucfilter_aspirin_exclusion_2layer_nested_cv_outputs.rds"
  ),
  # Preserve the exhaustive comparison and tuning grids in the supplied scripts.
  model_run_mode = "exhaustive",
  alpha_grid = c(0.25, 0.5, 0.75, 1),
  top_k_grid = c(20, 30, 50, 75, 100, 150, 200),
  # Historical parallel settings are preserved. Reduce workers for limited RAM.
  # Switching between parallel and sequential execution can change RNG streams.
  parallel_models = TRUE,
  max_model_workers = 10L,
  max_rf_threads = 10L,
  # The main profile evaluates MASLD subgroups using the existing OOF predictions.
  # Retraining within MASLD is a separately requested exploratory analysis.
  run_exploratory_masld_model = FALSE,
  run_shap = FALSE
)
