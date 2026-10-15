# Execute one named stage and retain private run metadata
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("This helper is called by run_pipeline.R with stage and MASLD flag.")
if (!args[[2L]] %in% c("true", "false")) stop("Invalid MASLD flag.")
source("R/configuration.R")
source("R/pipeline_registry.R")
registry <- ccrvi_registry(ccrvi_load_config())
id <- args[[1L]]
if (!id %in% names(registry)) stop("Unknown stage: ", id)
Sys.setenv(CCRVI_STAGE = id, CCRVI_MASLD_EXPLORATORY = args[[2L]])
# Source in a dedicated environment: runner helper objects cannot be overwritten
# by similarly named variables in an analysis script.
analysis_env <- new.env(parent = globalenv())
source("tools/stage_logging.R")
warning_file <- file.path("outputs/logs", paste0("latest_", id, "_warnings.csv"))
ccrvi_with_warning_log(
  function() source(registry[[id]]$script, local = analysis_env, echo = FALSE),
  warning_file
)
metadata <- analysis_env$ccrvi_run_metadata(stage = id)
metadata$completed_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
metadata$session_info <- utils::sessionInfo()
metadata$exploratory_masld_retraining <- identical(args[[2L]], "true")
metadata$warning_count <- nrow(utils::read.csv(warning_file))
saveRDS(metadata, file.path("outputs/logs", paste0("latest_", id, "_metadata.rds")))
# Export the actual global validation-fold membership when the script defines it.
# These identifiers and mappings are private, and are excluded from Git.
if (exists("validation_folds_global", envir = analysis_env, inherits = FALSE) &&
    exists("model_source", envir = analysis_env, inherits = FALSE) &&
    "Study.ID" %in% names(analysis_env$model_source)) {
  folds <- analysis_env$validation_folds_global
  ids <- analysis_env$model_source$Study.ID
  mapping <- do.call(rbind, lapply(seq_along(folds), function(k) {
    index <- folds[[k]]
    data.frame(Study.ID = ids[index], cohort_row = index, validation_fold = k)
  }))
  mapping <- mapping[order(mapping$cohort_row), ]
  utils::write.csv(mapping, file.path("outputs/logs", paste0(id, "_validation_fold_map.csv")), row.names = FALSE)
}
cat("Stage completed successfully: ", id, "\n", sep = "")
cat("Captured warnings: ", metadata$warning_count, "; details: ", warning_file, "\n", sep = "")
