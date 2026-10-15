# Inspect the files selected for a code submission without reading study data
# Use --all-files on an extracted publication archive. In an existing Git project
# the default inspects tracked/indexed files, allowing ignored data to stay local.
# --release additionally requires the actual licence and reproducible lockfile.
if (!file.exists("CCRVI.Rproj")) stop("Run from the project root.")
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--all-files", "--release"))) stop("Use --all-files and/or --release.")
all_files <- "--all-files" %in% args
paths <- NULL
if (!all_files && nzchar(Sys.which("git"))) {
  probe <- suppressWarnings(system2("git", c("rev-parse", "--is-inside-work-tree"), stdout = TRUE, stderr = TRUE))
  if (identical(probe, "true")) {
    paths <- system2("git", "ls-files", stdout = TRUE, stderr = TRUE)
    if (!length(paths)) stop("The Git index is empty. Stage only the intended code/documentation, then rerun; or use --all-files on the clean archive.")
  }
}
if (is.null(paths)) paths <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
paths <- gsub("\\\\", "/", paths)
paths <- paths[!grepl("(^|/)\\.git(/|$)", paths)]
allowed <- c("data/README.md", "data/raw/README.md", "data/processed/README.md", "outputs/README.md")
restricted_locations <- grepl("^(data/raw|data/processed|outputs|slide_exports)/", paths)
restricted_files <- grepl("\\.(csv|tsv|xlsx?|rds|rdata|zip|rda|parquet|feather)$", paths, ignore.case = TRUE)
local_state <- grepl("(^|/)(\\.Rhistory|\\.RData|\\.Ruserdata|\\.Renviron|\\.env(\\..*)?|local\\.R)$", paths)
unsafe <- paths[(restricted_locations | restricted_files | local_state) & !paths %in% allowed]
if (length(unsafe)) stop("Exclude these files before submission:\n", paste(unsafe, collapse = "\n"))
required <- c("README.md", "CCRVI.Rproj", "run_pipeline.R", ".gitignore",
              "docs/METHODS.md", "data/README.md", "renv.lock")
missing <- setdiff(required, paths)
if (length(missing)) stop("Required submission files are absent: ", paste(missing, collapse = ", "))
cat("PASS: ", length(paths), " selected files contain no restricted data/model files under these rules.\n", sep = "")
cat("This path/extension check does not replace review of document contents or Git history.\n")
blockers <- character()
if (!file.exists("LICENSE")) blockers <- c(blockers, "Add the agreed code LICENSE.")
if (!file.exists("renv.lock")) blockers <- c(blockers, "Capture the verified renv.lock from the working R environment.")
if (length(blockers)) {
  cat("Release items still pending:\n", paste0("  - ", blockers, collapse = "\n"), "\n", sep = "")
  if ("--release" %in% args) stop("Final release checks are incomplete.")
}
