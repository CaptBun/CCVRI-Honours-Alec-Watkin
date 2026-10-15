# Initialise one script in a fresh R process
# Run from the repository root so all scripts share the same relative paths.
if (!file.exists("CCRVI.Rproj") || !file.exists("config/default.R")) {
  stop("Run from the repository root (the folder containing CCRVI.Rproj).")
}
source("R/configuration.R", local = TRUE)
ccrvi_config <- ccrvi_load_config()
# Profiles can explicitly select subgroup evaluation or exploratory retraining.
masld_override <- Sys.getenv("CCRVI_MASLD_EXPLORATORY", unset = "")
if (nzchar(masld_override)) {
  ccrvi_config$run_exploratory_masld_model <- identical(masld_override, "true")
}
for (path in c("data/processed", "outputs/logs", "outputs/models")) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
for (key in c("compiled_cohort", "unfiltered_cohort", "original_results", "no_aspirin_results")) {
  dir.create(dirname(ccrvi_config$files[[key]]), recursive = TRUE, showWarnings = FALSE)
}

# Record metadata for a run without storing raw participant records in the log
ccrvi_run_metadata <- function(stage = Sys.getenv("CCRVI_STAGE", unset = "manual")) {
  tracked <- sort(c(list.files("R", pattern = "\\.R$", full.names = TRUE),
                    list.files("scripts", pattern = "\\.R$", full.names = TRUE),
                    "config/default.R"))
  if (file.exists("config/local.R")) tracked <- c(tracked, "config/local.R")
  if (file.exists("renv.lock")) tracked <- c(tracked, "renv.lock")
  input_paths <- unlist(ccrvi_config$files[c("clinical", "outcome", "workbook", "compiled_cohort")])
  present <- file.exists(input_paths)
  pkgs <- c("tidyverse", "caret", "glmnet", "pROC", "PRROC", "future", "furrr", "ranger")
  versions <- vapply(pkgs, function(pkg) {
    if (nzchar(system.file(package = pkg))) as.character(utils::packageVersion(pkg)) else "not installed"
  }, character(1))
  list(stage = stage, config = ccrvi_config,
       input_md5 = tools::md5sum(input_paths[present]),
       code_md5 = tools::md5sum(tracked),
       R_version = R.version.string, platform = R.version$platform,
       detected_cores = parallel::detectCores(logical = TRUE), rng_kind = RNGkind(),
       package_versions = versions)
}

# Use a separate checkpoint namespace whenever relevant run provenance changes
# MD5 is used only to identify changed files, not as a security mechanism.
# This avoids silently accepting results created from other data or settings.
ccrvi_checkpoint_directory <- function(base_path) {
  metadata <- ccrvi_run_metadata()
  scratch_file <- tempfile(fileext = ".rds")
  on.exit(unlink(scratch_file), add = TRUE)
  saveRDS(metadata, scratch_file, version = 2, compress = FALSE)
  fingerprint <- unname(tools::md5sum(scratch_file))
  path <- file.path(base_path, fingerprint)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  provenance_file <- file.path(path, "run_provenance.rds")
  if (!file.exists(provenance_file)) saveRDS(metadata, provenance_file)
  message("Checkpoint namespace: ", path)
  path
}
