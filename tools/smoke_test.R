# Exercise the actual elastic-net helper functions with wholly synthetic data
# This checks training-only preprocessing and function interoperability.
# It does not reproduce the study cohort or provide scientific validation.
if (!file.exists("CCRVI.Rproj")) stop("Run from the repository root.")
for (pkg in c("tidyverse", "glmnet", "pROC")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " before this test.")
}
library(tidyverse)
library(glmnet)
needed <- c("make_weights", "get_omics_features", "get_clinical_features", "select_auc_features", "prep_train_test", "fit_predict_enet")
helpers <- new.env(parent = globalenv())
# Evaluate only named function definitions, never the script's data reads or fits.
for (expression in parse("scripts/03_compare_modalities.R")) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% needed) {
    eval(expression, envir = helpers)
  }
}
stopifnot(all(vapply(needed, exists, logical(1), envir = helpers, inherits = FALSE)))
helpers$omics_prefixes <- c("lipo_", "lipid_", "sm_", "cyto_", "ba_", "path_")
helpers$use_glmnet_parallel <- FALSE
set.seed(123)
n <- 100L
y <- rep(c(0, 1), each = n / 2L)
synthetic <- data.frame(outcome = y, age = rnorm(n, 60, 8),
                        group = rep(c("A", "B"), length.out = n),
                        sm_signal = y + rnorm(n, 0, 0.3), sm_noise = rnorm(n),
                        lipid_noise = rnorm(n), constant = 1)
synthetic$age[c(3, 9)] <- NA_real_
# Use both classes in training and testing; no participant records are used.
training_index <- c(1:40, 51:90)
train <- synthetic[training_index, ]
test <- synthetic[-training_index, ]
test$group[[1L]] <- "Unseen category"
processed <- helpers$prep_train_test(train, test)
changed_test <- test
changed_test$age <- 1e8
changed_test$sm_signal <- -1e8
processed_changed <- helpers$prep_train_test(train, changed_test)
stopifnot(isTRUE(all.equal(processed$x_train, processed_changed$x_train)))
stopifnot(identical(colnames(processed$x_train), colnames(processed$x_test)))
stopifnot(all(is.finite(processed$x_train)), all(is.finite(processed$x_test)))
stopifnot(!"constant" %in% colnames(processed$x_train))
selected <- helpers$select_auc_features(processed$x_train, processed$y_train, top_k = 1)
stopifnot(identical(selected, "sm_signal"))
features <- unique(c(helpers$get_clinical_features(processed$x_train), selected))
fit <- helpers$fit_predict_enet(processed$x_train[, features, drop = FALSE], processed$y_train,
                                processed$x_test[, features, drop = FALSE], alpha_val = 0.5)
stopifnot(length(fit$pred) == nrow(test), all(is.finite(fit$pred)),
          all(fit$pred >= 0 & fit$pred <= 1), is.finite(fit$lambda_min))
cat("PASS: synthetic preprocessing, feature screening and weighted elastic-net fitting.\n")
