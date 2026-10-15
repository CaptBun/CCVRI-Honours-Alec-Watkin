# Load and validate settings without reading participant data
ccrvi_load_config <- function() {
  env <- new.env(parent = baseenv())
  sys.source("config/default.R", envir = env)
  if (file.exists("config/local.R")) sys.source("config/local.R", envir = env)
  cfg <- env$ccrvi_config
  if (!cfg$model_run_mode %in% c("screened", "pairs", "exhaustive")) {
    stop("model_run_mode must be screened, pairs or exhaustive.")
  }
  for (field in c("parallel_models", "run_exploratory_masld_model", "run_shap")) {
    if (!is.logical(cfg[[field]]) || length(cfg[[field]]) != 1L || is.na(cfg[[field]])) {
      stop(field, " must be a single TRUE/FALSE value.")
    }
  }
  for (field in c("max_model_workers", "max_rf_threads")) {
    value <- cfg[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != as.integer(value)) stop(field, " must be a positive integer.")
  }
  if (!is.numeric(cfg$alpha_grid) || !length(cfg$alpha_grid) ||
      anyNA(cfg$alpha_grid) || any(cfg$alpha_grid <= 0 | cfg$alpha_grid > 1)) {
    stop("alpha_grid must contain numeric values in (0, 1].")
  }
  if (!is.numeric(cfg$top_k_grid) || !length(cfg$top_k_grid) ||
      anyNA(cfg$top_k_grid) || any(!is.finite(cfg$top_k_grid)) ||
      any(cfg$top_k_grid < 1 | cfg$top_k_grid != floor(cfg$top_k_grid))) {
    stop("top_k_grid must contain positive integers.")
  }
  for (path in cfg$files) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
      stop("Every configured file path must be a non-empty character string.")
    }
  }
  if (length(unique(unlist(cfg$files))) != length(cfg$files)) {
    stop("Input and output paths must be distinct; check config/local.R.")
  }
  cfg
}
