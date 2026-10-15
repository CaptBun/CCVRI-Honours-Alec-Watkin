# Plan or execute ordered analysis stages from the project root
# Each stage runs in a fresh R session so objects from earlier scripts cannot
# silently become dependencies. By default this prints a plan and does no work.
if (!file.exists("CCRVI.Rproj")) stop("Run from the directory containing CCRVI.Rproj.")
source("R/configuration.R")
source("R/pipeline_registry.R")
cfg <- ccrvi_load_config()
registry <- ccrvi_registry(cfg)
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript run_pipeline.R [profile] [--execute] [--from stage] [--stage stage] [--list]\n")
  cat("Profiles: ", paste(names(ccrvi_profiles), collapse = ", "), "\n", sep = "")
  cat("The default profile is main. Execution requires --execute. --stage runs only that stage.\n")
  quit(status = 0L)
}
if ("--list" %in% args) {
  for (name in names(registry)) cat(name, " : ", registry[[name]]$description, "\n", sep = "")
  cat("\nProfiles:\n")
  for (name in names(ccrvi_profiles)) cat(name, " : ", paste(ccrvi_profiles[[name]], collapse = " -> "), "\n", sep = "")
  quit(status = 0L)
}
execute <- FALSE
profile <- "main"
profile_seen <- FALSE
from <- NULL
single_stage <- NULL
i <- 1L
while (i <= length(args)) {
  token <- args[[i]]
  if (token == "--execute") {
    execute <- TRUE
  } else if (token %in% c("--from", "--stage")) {
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) stop("A stage name is required after ", token)
    i <- i + 1L
    if (token == "--from") from <- args[[i]] else single_stage <- args[[i]]
  } else if (startsWith(token, "--")) {
    stop("Unknown option: ", token)
  } else {
    if (profile_seen) stop("Only one profile may be specified.")
    profile <- token
    profile_seen <- TRUE
  }
  i <- i + 1L
}
if (!profile %in% names(ccrvi_profiles)) stop("Unknown profile: ", profile)
if (!is.null(from) && !is.null(single_stage)) stop("Use either --from or --stage.")
stages <- ccrvi_profiles[[profile]]
if (!is.null(single_stage)) stages <- single_stage
if (any(!stages %in% names(registry))) stop("Unknown stage: ", paste(setdiff(stages, names(registry)), collapse = ", "))
if (!is.null(from)) {
  start <- match(from, stages)
  if (is.na(start)) stop("--from stage is not in the selected profile.")
  stages <- stages[seq.int(start, length(stages))]
}
# The all and masld-exploratory profiles explicitly request subgroup retraining.
# Other profiles respect the user's configured setting (FALSE by default).
masld_exploratory <- if (profile %in% c("all", "masld-exploratory")) TRUE else cfg$run_exploratory_masld_model
cat("Profile: ", profile, "\nMode: ", if (execute) "EXECUTE" else "PLAN ONLY", "\n", sep = "")
cat("Comparison scope: ", cfg$model_run_mode, "; parallel requested: ", cfg$parallel_models,
    "; max workers: ", cfg$max_model_workers, "\n", sep = "")
cat("MASLD-only exploratory retraining: ", masld_exploratory, "\n\n", sep = "")
for (id in stages) {
  item <- registry[[id]]
  cat(id, " : ", item$script, "\n", sep = "")
  cat("  Requires: ", paste(item$inputs, collapse = ", "), "\n", sep = "")
}
# Determine which dependencies must already exist; later inputs may be produced
# by earlier stages in this same invocation.
available <- character()
missing <- character()
for (id in stages) {
  item <- registry[[id]]
  needed <- setdiff(item$inputs, available)
  missing <- union(missing, needed[!file.exists(needed)])
  available <- union(available, item$outputs)
}
if (length(missing)) {
  cat("\nInputs currently missing:\n", paste0("  ", missing, collapse = "\n"), "\n", sep = "")
}
if (!execute) {
  cat("\nPlan complete. Add --execute to run; no analysis or package installation occurred.\n")
  quit(status = 0L)
}
if (length(missing)) stop("Supply the missing inputs before executing this plan.")
packages <- unique(unlist(lapply(registry[stages], function(x) x$packages)))
if (cfg$parallel_models && any(stages %in% c("original", "no-aspirin", "sensitivity"))) {
  # Refuse a silent switch to the sequential RNG branch in the orchestrated run.
  packages <- union(packages, c("future", "furrr"))
}
if (cfg$run_shap && "candidate" %in% stages) packages <- union(packages, "fastshap")
missing_packages <- packages[!vapply(packages, function(pkg) nzchar(system.file(package = pkg)), logical(1))]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "),
                                 ". Restore the supplied package environment with renv::restore(); see README.md.")
# This code itself never installs or updates packages during analysis execution.
dir.create("outputs/logs", recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
for (id in stages) {
  item <- registry[[id]]
  absent <- item$inputs[!file.exists(item$inputs)]
  if (length(absent)) stop("Missing input before stage ", id, ": ", paste(absent, collapse = ", "))
  log_path <- file.path("outputs/logs", paste0(stamp, "_", id, "_", Sys.getpid(), ".log"))
  cat("\nRunning ", id, ". Log: ", log_path, "\n", sep = "")
  # --no-save/--no-restore preserve normal .Rprofile activation for renv.
  # All configured paths and profiles are reloaded independently in the child.
  status <- system2(rscript,
                    args = c("--no-save", "--no-restore", "tools/run_stage.R", shQuote(id),
                             if (masld_exploratory) "true" else "false"),
                    stdout = log_path, stderr = log_path)
  if (status != 0L) stop("Stage ", id, " failed. Inspect ", log_path)
  absent <- item$outputs[!file.exists(item$outputs)]
  if (length(absent)) stop("Stage ", id, " did not create its expected output: ", paste(absent, collapse = ", "))
  cat("Completed ", id, "\n", sep = "")
}
cat("\nRequested stages completed. Review outputs and compare with your reference results.\n")
