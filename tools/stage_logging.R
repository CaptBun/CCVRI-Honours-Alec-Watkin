# Capture warning conditions without changing model fitting or random-number use
# The supplied analysis files retain their own warning behaviour. This wrapper
# records conditions that reach the main R process, including re-signalled
# future-worker warnings, so a successful exit does not hide their details.
ccrvi_with_warning_log <- function(fun, warning_file) {
  records <- list()
  dir.create(dirname(warning_file), recursive = TRUE, showWarnings = FALSE)
  on.exit({
    if (length(records)) {
      table <- do.call(rbind, records)
    } else {
      table <- data.frame(
        warning_number = integer(), message = character(), call = character(),
        stringsAsFactors = FALSE
      )
    }
    utils::write.csv(table, warning_file, row.names = FALSE)
  }, add = TRUE)
  withCallingHandlers(
    fun(),
    warning = function(condition) {
      expression <- conditionCall(condition)
      records[[length(records) + 1L]] <<- data.frame(
        warning_number = length(records) + 1L,
        message = conditionMessage(condition),
        call = if (is.null(expression)) "" else paste(deparse(expression), collapse = " "),
        stringsAsFactors = FALSE
      )
      # Do not muffle or convert warnings into errors: retain the original flow.
    }
  )
}
