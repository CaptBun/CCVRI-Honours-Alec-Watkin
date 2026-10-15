# Repository stage: risk
# Source supplied as: 06_mace_risk_score_three_class_framework(5).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# 06 Exploratory MACE risk score and three-class stratification
#
# Project:
#   Cardiovascular Risk Stratification in Liver Cirrhosis
#
# Candidate model:
#   Expanded clinical (aspirin excluded) +
#   Lipoproteins + Lipids + Small molecules
#
# Purpose:
#   Translate the candidate model's pooled out-of-fold scores into:
#
#     1. A descriptive score from 0 to 100.
#     2. Three exploratory predicted-risk classes.
#     3. Class-level observed MACE rates and uncertainty.
#     4. Figures and tables suitable for a framework proposal.
#
# Important interpretation:
#   - The 0-100 score is a cohort-relative percentile score.
#   - It is not an absolute predicted probability.
#   - The Higher-risk boundary is defined by the pooled Youden J threshold.
#   - The Lower/Intermediate boundary is distributional only.
#   - None of these classes should be treated as clinically validated.
#
# Risk-class rules:
#   Lower predicted risk:
#     Bottom third of the cohort score distribution and below Youden J.
#
#   Intermediate predicted risk:
#     Above the lower distributional boundary but below Youden J.
#
#   Higher predicted risk:
#     At or above the pooled Youden J threshold.
#
# Required input:
#   outputs/models/lead_candidate_no_aspirin/
#     05_no_aspirin_candidate_outputs.rds
#
# Main outputs:
#   outputs/tables/risk_score/
#   outputs/figures/risk_score/
#   outputs/models/risk_score/

library(tidyverse)
library(pROC)

set.seed(123)


# Settings


candidate_output_file <- paste0(
  "outputs/models/lead_candidate_no_aspirin/",
  "05_no_aspirin_candidate_outputs.rds"
)

output_table_dir <- "outputs/tables/risk_score"
output_figure_dir <- "outputs/figures/risk_score"
output_model_dir <- "outputs/models/risk_score"

# Distributional boundary used only to distinguish Lower from Intermediate.
lower_risk_percentile_cutoff <- 33

# Number of bootstrap samples used to describe uncertainty in class-level
# observed event rates and in the Youden threshold.
n_bootstrap <- 2000
bootstrap_seed <- 123

dir.create(
  output_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  output_model_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Save helpers


save_table <- function(x, filename) {
  write.csv(
    x,
    file.path(output_table_dir, filename),
    row.names = FALSE
  )
}

save_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(output_figure_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}


# Statistical helpers


exact_binomial_ci <- function(events, total, conf_level = 0.95) {
  if (total == 0) {
    return(
      tibble(
        rate = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    )
  }

  fit <- binom.test(
    x = events,
    n = total,
    conf.level = conf_level
  )

  tibble(
    rate = events / total,
    ci_low = unname(fit$conf.int[1]),
    ci_high = unname(fit$conf.int[2])
  )
}

calculate_youden_threshold <- function(outcome, score) {
  roc_object <- pROC::roc(
    response = outcome,
    predictor = score,
    direction = "<",
    quiet = TRUE
  )

  best <- pROC::coords(
    roc_object,
    x = "best",
    best.method = "youden",
    ret = c(
      "threshold",
      "sensitivity",
      "specificity",
      "ppv",
      "npv"
    ),
    transpose = FALSE
  )

  tibble(
    threshold = as.numeric(best$threshold),
    sensitivity = as.numeric(best$sensitivity),
    specificity = as.numeric(best$specificity),
    ppv = as.numeric(best$ppv),
    npv = as.numeric(best$npv),
    youden_j = as.numeric(
      best$sensitivity + best$specificity - 1
    )
  )
}

# Convert a model score to a cohort-relative percentile score from 0 to 100.
#
# This preserves ranking but deliberately does not claim that the number is an
# absolute MACE probability.
make_percentile_score <- function(x) {
  if (length(x) <= 1) {
    return(rep(50, length(x)))
  }

  100 * (
    rank(
      x,
      ties.method = "average",
      na.last = "keep"
    ) - 1
  ) / (
    sum(!is.na(x)) - 1
  )
}


# Load Script 05 candidate outputs


if (!file.exists(candidate_output_file)) {
  stop(
    "Candidate output file not found: ",
    candidate_output_file
  )
}

candidate_outputs <- readRDS(candidate_output_file)

required_objects <- c(
  "predictions",
  "performance_summary",
  "settings"
)

missing_objects <- setdiff(
  required_objects,
  names(candidate_outputs)
)

if (length(missing_objects) > 0) {
  stop(
    "Script 05 output is missing required object(s): ",
    paste(missing_objects, collapse = ", ")
  )
}

candidate_predictions <- candidate_outputs$predictions %>%
  as_tibble() %>%
  mutate(
    outcome = as.numeric(outcome),
    pred = as.numeric(pred)
  )

if (!all(c("outcome", "pred") %in% names(candidate_predictions))) {
  stop(
    "The predictions object must contain outcome and pred columns."
  )
}

if (!"patient_index" %in% names(candidate_predictions)) {
  candidate_predictions <- candidate_predictions %>%
    mutate(patient_index = row_number())
}

candidate_label <- candidate_outputs$settings$candidate_label

# Use the exact threshold already reported in Script 05 where available.
candidate_threshold <- candidate_outputs$performance_summary$threshold[1]

if (
  length(candidate_threshold) == 0 ||
    is.na(candidate_threshold) ||
    !is.finite(candidate_threshold)
) {
  warning(
    "No valid threshold found in Script 05. ",
    "Recalculating pooled Youden J threshold."
  )

  candidate_threshold <- calculate_youden_threshold(
    candidate_predictions$outcome,
    candidate_predictions$pred
  )$threshold
}

message("Candidate model: ", candidate_label)
message(
  "Pooled Youden J threshold: ",
  round(candidate_threshold, 4)
)


# Construct the 0-100 score and risk classes

risk_score_data <- candidate_predictions %>%
  mutate(
    risk_score_100_raw = make_percentile_score(pred),
    risk_score_100 = round(risk_score_100_raw),
    above_youden = pred >= candidate_threshold,
    risk_class = case_when(
      above_youden ~ "Higher predicted risk",
      risk_score_100_raw <= lower_risk_percentile_cutoff ~
        "Lower predicted risk",
      TRUE ~ "Intermediate predicted risk"
    ),
    risk_class = factor(
      risk_class,
      levels = c(
        "Lower predicted risk",
        "Intermediate predicted risk",
        "Higher predicted risk"
      )
    ),
    observed_outcome = if_else(
      outcome == 1,
      "MACE",
      "No MACE"
    )
  )

# Translate the raw Youden threshold to its approximate percentile-score
# location for communication only. Classification continues to use raw pred.
youden_score_100 <- 100 * mean(
  risk_score_data$pred < candidate_threshold,
  na.rm = TRUE
)

if (youden_score_100 <= lower_risk_percentile_cutoff) {
  warning(
    "The Youden threshold lies at or below the Lower-risk percentile ",
    "boundary. The Intermediate group may be very small or absent."
  )
}


# Threshold summary


youden_metrics <- calculate_youden_threshold(
  risk_score_data$outcome,
  risk_score_data$pred
)

threshold_summary <- youden_metrics %>%
  mutate(
    candidate_model = candidate_outputs$settings$candidate_model,
    candidate_label = candidate_label,
    raw_model_threshold = candidate_threshold,
    approximate_score_100_threshold = youden_score_100,
    lower_risk_score_boundary =
      lower_risk_percentile_cutoff,
    lower_boundary_basis =
      "Bottom third of cohort score distribution",
    higher_boundary_basis =
      "Pooled Youden J threshold",
    score_interpretation =
      "Cohort-relative percentile score; not absolute probability",
    class_interpretation =
      "Exploratory framework; not clinically validated"
  ) %>%
  select(
    candidate_model,
    candidate_label,
    raw_model_threshold,
    approximate_score_100_threshold,
    lower_risk_score_boundary,
    sensitivity,
    specificity,
    ppv,
    npv,
    youden_j,
    lower_boundary_basis,
    higher_boundary_basis,
    score_interpretation,
    class_interpretation
  )

save_table(
  threshold_summary,
  "06A_risk_score_threshold_summary.csv"
)


# Patient-level score export


patient_score_export <- risk_score_data %>%
  select(
    patient_index,
    outcome,
    observed_outcome,
    pred,
    risk_score_100_raw,
    risk_score_100,
    risk_class,
    above_youden
  ) %>%
  arrange(desc(risk_score_100_raw))

save_table(
  patient_score_export,
  "06B_patient_level_risk_scores.csv"
)

# Save the empirical reference distribution required to score future patients
# on the same cohort-relative percentile scale.
score_reference_distribution <- risk_score_data %>%
  arrange(pred) %>%
  transmute(
    reference_rank = row_number(),
    reference_raw_model_score = pred,
    reference_percentile_score = risk_score_100_raw
  )

save_table(
  score_reference_distribution,
  "06B_score_reference_distribution.csv"
)


# Risk-class summary


risk_class_summary <- risk_score_data %>%
  group_by(risk_class) %>%
  summarise(
    n = n(),
    observed_mace_n = sum(outcome == 1),
    observed_no_mace_n = sum(outcome == 0),
    observed_mace_rate = mean(outcome == 1),
    score_min = min(risk_score_100_raw),
    score_q1 = quantile(risk_score_100_raw, 0.25),
    score_median = median(risk_score_100_raw),
    score_q3 = quantile(risk_score_100_raw, 0.75),
    score_max = max(risk_score_100_raw),
    raw_model_score_min = min(pred),
    raw_model_score_median = median(pred),
    raw_model_score_max = max(pred),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    ci = list(
      exact_binomial_ci(
        observed_mace_n,
        n
      )
    )
  ) %>%
  unnest(ci) %>%
  ungroup() %>%
  mutate(
    observed_mace_rate = rate
  ) %>%
  select(-rate)

save_table(
  risk_class_summary,
  "06C_risk_class_summary.csv"
)

# Compact presentation table
risk_class_presentation <- risk_class_summary %>%
  mutate(
    `Class definition` = case_when(
      risk_class == "Lower predicted risk" ~
        "Score <= 33",
      risk_class == "Intermediate predicted risk" ~
        paste0(
          "Score > 33 and raw model score < ",
          round(candidate_threshold, 3)
        ),
      risk_class == "Higher predicted risk" ~
        paste0(
          "Raw model score >= ",
          round(candidate_threshold, 3),
          " (~score ",
          round(youden_score_100, 1),
          ")"
        )
    )
  ) %>%
  transmute(
    `Risk class` = risk_class,
    `Class definition`,
    `Patients, n` = n,
    `Median score` = round(score_median),
    `MACE events, n` = observed_mace_n,
    `Observed MACE rate, %` =
      round(100 * observed_mace_rate, 1),
    `Exact 95% CI, %` = paste0(
      round(100 * ci_low, 1),
      "-",
      round(100 * ci_high, 1)
    )
  )

save_table(
  risk_class_presentation,
  "06C_risk_class_presentation.csv"
)


# Ordinal trend summaries


risk_trend_fit <- glm(
  outcome ~ risk_score_100_raw,
  data = risk_score_data,
  family = binomial
)

risk_trend_summary <- broom::tidy(
  risk_trend_fit,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  filter(term == "risk_score_100_raw") %>%
  transmute(
    score_increment = "Per 10-point increase",
    odds_ratio = estimate^10,
    ci_low = conf.low^10,
    ci_high = conf.high^10,
    p_value = p.value
  )

save_table(
  risk_trend_summary,
  "06D_risk_score_ordinal_trend.csv"
)


# Bootstrap uncertainty in the Youden threshold

set.seed(bootstrap_seed)

bootstrap_thresholds <- map_dfr(
  seq_len(n_bootstrap),
  function(bootstrap_i) {
    bootstrap_indices <- sample(
      seq_len(nrow(risk_score_data)),
      size = nrow(risk_score_data),
      replace = TRUE
    )

    bootstrap_data <- risk_score_data[
      bootstrap_indices,
      ,
      drop = FALSE
    ]

    if (n_distinct(bootstrap_data$outcome) < 2) {
      return(
        tibble(
          bootstrap = bootstrap_i,
          raw_threshold = NA_real_
        )
      )
    }

    threshold_i <- tryCatch(
      calculate_youden_threshold(
        bootstrap_data$outcome,
        bootstrap_data$pred
      )$threshold,
      error = function(e) NA_real_
    )

    tibble(
      bootstrap = bootstrap_i,
      raw_threshold = threshold_i
    )
  }
) %>%
  filter(
    !is.na(raw_threshold),
    is.finite(raw_threshold)
  ) %>%
  mutate(
    approximate_score_100 = map_dbl(
      raw_threshold,
      ~ 100 * mean(
        risk_score_data$pred < .x,
        na.rm = TRUE
      )
    )
  )

save_table(
  bootstrap_thresholds,
  "06E_bootstrap_youden_thresholds.csv"
)

bootstrap_threshold_summary <- bootstrap_thresholds %>%
  summarise(
    n_successful_bootstrap = n(),
    raw_threshold_median =
      median(raw_threshold),
    raw_threshold_ci_low =
      quantile(raw_threshold, 0.025),
    raw_threshold_ci_high =
      quantile(raw_threshold, 0.975),
    score_100_median =
      median(approximate_score_100),
    score_100_ci_low =
      quantile(approximate_score_100, 0.025),
    score_100_ci_high =
      quantile(approximate_score_100, 0.975)
  )

save_table(
  bootstrap_threshold_summary,
  "06E_bootstrap_youden_threshold_summary.csv"
)


# Figures


# 1. Patient scores by observed outcome and exploratory risk class
score_distribution_plot <- ggplot(
  risk_score_data,
  aes(
    x = observed_outcome,
    y = risk_score_100_raw,
    fill = risk_class
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.35
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    alpha = 0.75
  ) +
  geom_jitter(
    width = 0.08,
    alpha = 0.45,
    size = 1.5
  ) +
  geom_hline(
    yintercept = youden_score_100,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 1.5,
    y = youden_score_100,
    label = paste0(
      "Youden J boundary ≈ score ",
      round(youden_score_100)
    ),
    vjust = -0.7,
    size = 3.6
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20)
  ) +
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  labs(
    title = "Cohort-relative MACE risk score",
    subtitle = candidate_label,
    x = "Observed outcome",
    y = "Risk score (0-100)",
    fill = "Predicted-risk class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

save_plot(
  score_distribution_plot,
  "06A_patient_risk_score_distribution.png",
  width = 8,
  height = 6
)

# 2. Observed MACE rate by class with exact binomial confidence intervals
risk_class_rate_plot <- ggplot(
  risk_class_summary,
  aes(
    x = risk_class,
    y = 100 * observed_mace_rate,
    fill = risk_class
  )
) +
  geom_col(
    width = 0.68,
    alpha = 0.85
  ) +
  geom_errorbar(
    aes(
      ymin = 100 * ci_low,
      ymax = 100 * ci_high
    ),
    width = 0.16,
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = paste0(
        observed_mace_n,
        "/",
        n,
        "\n",
        round(100 * observed_mace_rate, 1),
        "%"
      )
    ),
    vjust = -0.35,
    size = 4
  ) +
  expand_limits(
    y = min(
      100,
      max(100 * risk_class_summary$ci_high) + 12
    )
  ) +
  labs(
    title = "Observed MACE rate by predicted-risk class",
    subtitle = paste0(
      "Higher-risk boundary uses pooled Youden J; ",
      "lower boundary is distributional"
    ),
    x = NULL,
    y = "Observed MACE rate (%)",
    fill = NULL,
    caption = paste0(
      "Exact binomial 95% confidence intervals. ",
      "Exploratory framework only; classes are not clinically validated."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 15,
      hjust = 1
    )
  )

save_plot(
  risk_class_rate_plot,
  "06B_observed_mace_rate_by_risk_class.png",
  width = 8,
  height = 6
)

# 3. Score-ranked patient plot
ranked_score_plot <- risk_score_data %>%
  arrange(risk_score_100_raw) %>%
  mutate(
    score_rank = row_number()
  ) %>%
  ggplot(
    aes(
      x = score_rank,
      y = risk_score_100_raw,
      colour = observed_outcome,
      shape = risk_class
    )
  ) +
  geom_hline(
    yintercept = lower_risk_percentile_cutoff,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = youden_score_100,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_point(
    alpha = 0.8,
    size = 2
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20)
  ) +
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  labs(
    title = "Ranked patient risk scores",
    subtitle = paste0(
      "Dotted line: lower distributional boundary; ",
      "dashed line: pooled Youden J boundary"
    ),
    x = "Patients ordered by score",
    y = "Risk score (0-100)",
    colour = "Observed outcome",
    shape = "Predicted-risk class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

save_plot(
  ranked_score_plot,
  "06C_ranked_patient_risk_scores.png",
  width = 10,
  height = 6
)

# 4. Bootstrap distribution of the Youden boundary on the 0-100 scale
bootstrap_threshold_plot <- ggplot(
  bootstrap_thresholds,
  aes(
    x = approximate_score_100
  )
) +
  geom_histogram(
    bins = 30,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = youden_score_100,
    linetype = "dashed",
    linewidth = 0.9
  ) +
  labs(
    title = "Bootstrap uncertainty in the Youden J boundary",
    subtitle = paste0(
      "Dashed line shows the pooled threshold at score ",
      round(youden_score_100)
    ),
    x = "Approximate threshold on 0-100 score scale",
    y = "Bootstrap samples"
  ) +
  theme_minimal(base_size = 12)

save_plot(
  bootstrap_threshold_plot,
  "06D_bootstrap_youden_score_boundary.png",
  width = 8,
  height = 5
)


# Save combined Script 06 output


risk_score_outputs <- list(
  settings = list(
    candidate_model =
      candidate_outputs$settings$candidate_model,
    candidate_label =
      candidate_label,
    source_file =
      candidate_output_file,
    score_method =
      "Empirical percentile rank of pooled out-of-fold model score",
    score_range =
      "0-100",
    lower_risk_percentile_cutoff =
      lower_risk_percentile_cutoff,
    higher_risk_threshold_method =
      "Pooled Youden J",
    interpretation =
      "Exploratory framework; not clinically validated"
  ),
  threshold_summary =
    threshold_summary,
  bootstrap_threshold_summary =
    bootstrap_threshold_summary,
  bootstrap_thresholds =
    bootstrap_thresholds,
  patient_scores =
    risk_score_data,
  score_reference_distribution =
    score_reference_distribution,
  risk_class_summary =
    risk_class_summary,
  risk_trend_summary =
    risk_trend_summary
)

saveRDS(
  risk_score_outputs,
  file.path(
    output_model_dir,
    "06_mace_risk_score_outputs.rds"
  )
)

message("Script 06 risk-score analysis complete.")
message(
  "Approximate Youden boundary on 0-100 scale: ",
  round(youden_score_100, 1)
)
message(
  "Lower-risk descriptive boundary: ",
  lower_risk_percentile_cutoff
)
message("Risk-class summary:")
print(risk_class_presentation)
