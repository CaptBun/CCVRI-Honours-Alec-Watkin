# Copy this file to config/local.R and change values for your computer
# local.R is ignored by Git. Keep the default paths if using data/raw/.
# Values override config/default.R without changing the analysis scripts.
ccrvi_config$max_model_workers <- 2L
ccrvi_config$max_rf_threads <- 2L
# Use TRUE to reproduce the historical parallel branch when future/furrr exist.
ccrvi_config$parallel_models <- TRUE
# Example: keep a restricted input outside the repository.
# ccrvi_config$files$clinical <- "D:/restricted_data/clindata.csv"
# A screened run is for exploration and will not reproduce the final model panel.
# ccrvi_config$model_run_mode <- "screened"
