# Validate project wiring and R syntax without loading private study data
if (!file.exists("CCRVI.Rproj")) stop("Run from the repository root.")
files <- sort(c(list.files("R", pattern = "\\.R$", full.names = TRUE),
                list.files("scripts", pattern = "\\.R$", full.names = TRUE),
                list.files("tools", pattern = "\\.R$", full.names = TRUE),
                list.files("config", pattern = "\\.R$", full.names = TRUE),
                "run_pipeline.R"))
for (path in files) parse(file = path)
source("tools/stage_logging.R")
# Confirm that warnings are logged on both success and failure without changing
# the expression's return value or suppressing its error condition.
test_directory <- tempfile("ccrvi_validation_")
dir.create(test_directory)
success_log <- file.path(test_directory, "success.csv")
value <- suppressWarnings(ccrvi_with_warning_log(
  function() { warning("Synthetic validation warning"); 7L }, success_log
))
stopifnot(identical(value, 7L), nrow(utils::read.csv(success_log)) == 1L)
failure_log <- file.path(test_directory, "failure.csv")
failure <- tryCatch(suppressWarnings(ccrvi_with_warning_log(
  function() { warning("Synthetic warning before failure"); stop("Synthetic failure") },
  failure_log
)), error = function(e) e)
stopifnot(inherits(failure, "error"), nrow(utils::read.csv(failure_log)) == 1L)
unlink(test_directory, recursive = TRUE)
source("R/configuration.R")
source("R/pipeline_registry.R")
cfg <- ccrvi_load_config()
registry <- ccrvi_registry(cfg)
stopifnot(all(vapply(registry, function(x) file.exists(x$script), logical(1))))
stopifnot(all(unlist(ccrvi_profiles) %in% names(registry)))
# The main profile should supply each generated dependency before it is needed.
external <- unlist(cfg$files[c("clinical", "outcome", "workbook")])
produced <- external
for (id in ccrvi_profiles$main) {
  stage <- registry[[id]]
  unmet <- setdiff(stage$inputs, produced)
  if (length(unmet)) stop("Unmet dependency in main profile at ", id, ": ", paste(unmet, collapse = ", "))
  produced <- union(produced, stage$outputs)
}
# Test profile planning in child sessions, and ensure planning creates no outputs.
script <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
for (arguments in list(c("run_pipeline.R", "--list"),
                       c("run_pipeline.R", "main"),
                       c("run_pipeline.R", "main", "--from", "candidate"),
                       c("run_pipeline.R", "--stage", "rf"),
                       c("run_pipeline.R", "masld-exploratory"))) {
  result <- system2(script, args = arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) stop(paste(result, collapse = "\n"))
}
cat("PASS: parsed ", length(files), " R files; warning capture, main dependencies and runner plans are valid.\n", sep = "")
cat("This check does not establish numerical agreement with the study results.\n")
