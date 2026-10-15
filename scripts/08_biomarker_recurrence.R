# Repository stage: biomarkers
# Source supplied as: 08_cross_model_biomarker_recurrence_analysis_final(3).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# 08 Cross-model biomarker recurrence, persistence and overlap analysis
#
# Project:
#   Cardiovascular Risk Stratification in Liver Cirrhosis
#
# Purpose:
#   Compare predictor selection, coefficient direction and stability across a
#   prespecified panel of clinically important models:
#
#   1. Clinical-only models.
#   2. Each individual omics-only model.
#   3. The locked candidate model.
#   4. A maximal aspirin-excluded model containing all omics modalities.
#
# Scientific interpretation:
#   This script does not train or tune new models. It reuses feature-selection
#   summaries generated from the outer folds of the nested cross-validation
#   pipeline. A predictor repeatedly selected across folds and across related
#   model specifications is treated as a more reproducible model signal than a
#   predictor selected in only one model or one fold.
#
#   Coefficients are estimated after fold-specific standardisation, but direct
#   magnitude comparisons across models remain approximate because each model
#   has a different predictor set, correlation structure and penalty. The main
#   cross-model evidence is therefore:
#
#     - selection frequency within each model
#     - presence across models
#     - coefficient-direction agreement
#     - rank consistency
#
# Required inputs:
#   outputs/models/aucfilter_multimodality_comparison_2layer_nested_cv_outputs.rds
#   outputs/models/aucfilter_aspirin_exclusion_2layer_nested_cv_outputs.rds
#
# Main outputs:
#   outputs/tables/biomarker_recurrence/
#   outputs/figures/biomarker_recurrence/
#   outputs/models/biomarker_recurrence/

library(tidyverse)

set.seed(123)

# Settings
original_output_candidates <- c(ccrvi_config$files$original_results)

no_aspirin_output_candidates <- c(ccrvi_config$files$no_aspirin_results)

output_table_dir <- "outputs/tables/biomarker_recurrence"
output_figure_dir <- "outputs/figures/biomarker_recurrence"
output_model_dir <- "outputs/models/biomarker_recurrence"

dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_model_dir, recursive = TRUE, showWarnings = FALSE)

# Models are deliberately prespecified rather than selected after inspecting
# biomarker results. This limits post-hoc cherry-picking of favourable features.
model_registry <- tribble(
  ~model, ~source, ~model_group, ~model_label, ~display_order,
  "clinical_baseline",
  "original",
  "Clinical only",
  "Baseline clinical",
  1,

  "clinical_expanded",
  "original",
  "Clinical only",
  "Expanded clinical",
  2,

  "clinical_expanded_no_aspirin",
  "no_aspirin",
  "Clinical only",
  "Expanded clinical\n(aspirin excluded)",
  3,

  "lipoproteins",
  "original",
  "Omics only",
  "Lipoproteins only",
  4,

  "lipids",
  "original",
  "Omics only",
  "Lipids only",
  5,

  "small_molecules",
  "original",
  "Omics only",
  "Small molecules only",
  6,

  "cytokines",
  "original",
  "Omics only",
  "Cytokines only",
  7,

  "bile_acids",
  "original",
  "Omics only",
  "Bile acids only",
  8,

  "pathology",
  "original",
  "Omics only",
  "Pathology only",
  9,

  paste0(
    "clinical_expanded_no_aspirin_",
    "lipoproteins_lipids_small_molecules"
  ),
  "no_aspirin",
  "Integrated",
  "Locked candidate",
  10,

  paste0(
    "clinical_expanded_no_aspirin_",
    "lipoproteins_lipids_small_molecules_",
    "cytokines_bile_acids_pathology"
  ),
  "no_aspirin",
  "Integrated",
  "All modalities",
  11
)

# Reuse the exact saved model identifiers throughout the script. Defining these
# objects once prevents literal-name errors in later overlap analyses.
candidate_model <- paste0(
  "clinical_expanded_no_aspirin_",
  "lipoproteins_lipids_small_molecules"
)

all_modalities_model <- paste0(
  "clinical_expanded_no_aspirin_",
  "lipoproteins_lipids_small_molecules_",
  "cytokines_bile_acids_pathology"
)

# Minimum within-model stability required for the main cross-model analyses.
# Three of five outer folds corresponds to a frequency of 0.60.
stable_selection_cutoff <- 0.60

# The heatmap is restricted to the most reproducible features so labels remain
# legible. Complete results are retained in the exported tables.
top_n_heatmap_features <- 35
top_n_recurrent_features <- 25
top_n_flow_features <- 18

# Save helpers
save_table <- function(x, filename) {
  write.csv(
    x,
    file.path(output_table_dir, filename),
    row.names = FALSE
  )
}

save_plot <- function(plot_object, filename, width = 10, height = 8) {
  ggsave(
    filename = file.path(output_figure_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

# General helpers
find_first_existing_file <- function(candidate_paths, description) {
  existing_paths <- candidate_paths[file.exists(candidate_paths)]

  if (length(existing_paths) == 0) {
    stop(
      description,
      " was not found. Expected one of:\n",
      paste(candidate_paths, collapse = "\n")
    )
  }

  existing_paths[1]
}

classify_feature_modality <- function(feature_name) {
  case_when(
    str_starts(feature_name, "lipo_") ~ "Lipoproteins",
    str_starts(feature_name, "lipid_") ~ "Lipids",
    str_starts(feature_name, "sm_") ~ "Small molecules",
    str_starts(feature_name, "cyto_") ~ "Cytokines",
    str_starts(feature_name, "ba_") ~ "Bile acids",
    str_starts(feature_name, "path_") ~ "Pathology",
    TRUE ~ "Clinical"
  )
}

clean_feature_label <- function(feature_name) {
  feature_name %>%
    str_remove("^lipo_") %>%
    str_remove("^lipid_") %>%
    str_remove("^sm_") %>%
    str_remove("^cyto_") %>%
    str_remove("^ba_") %>%
    str_remove("^path_") %>%
    str_replace_all("\\.", " ") %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

safe_divide <- function(numerator, denominator) {
  ifelse(
    is.na(denominator) | denominator == 0,
    NA_real_,
    numerator / denominator
  )
}

# Validate and standardise one Script 03 output object.
extract_source_tables <- function(output_object, source_label) {
  required_objects <- c(
    "modality_feature_votes",
    "pooled_modality_performance"
  )

  missing_objects <- setdiff(required_objects, names(output_object))

  if (length(missing_objects) > 0) {
    stop(
      "The ", source_label, " output object is missing: ",
      paste(missing_objects, collapse = ", ")
    )
  }

  feature_votes <- output_object$modality_feature_votes %>%
    as_tibble() %>%
    mutate(
      model = as.character(model),
      feature = as.character(feature),
      source = source_label
    )

  required_vote_columns <- c(
    "model",
    "feature",
    "selected_n",
    "selection_frequency",
    "mean_coefficient",
    "median_coefficient",
    "mean_abs_coefficient",
    "direction"
  )

  missing_vote_columns <- setdiff(
    required_vote_columns,
    names(feature_votes)
  )

  if (length(missing_vote_columns) > 0) {
    stop(
      "The ", source_label,
      " feature-vote table is missing column(s): ",
      paste(missing_vote_columns, collapse = ", ")
    )
  }

  performance <- output_object$pooled_modality_performance %>%
    as_tibble() %>%
    mutate(
      model = as.character(model),
      source = source_label
    )

  list(
    feature_votes = feature_votes,
    performance = performance
  )
}

# Load model outputs
original_output_file <- find_first_existing_file(
  original_output_candidates,
  "Original Script 03 output"
)

no_aspirin_output_file <- find_first_existing_file(
  no_aspirin_output_candidates,
  "Aspirin-excluded Script 03 output"
)

original_tables <- readRDS(original_output_file) %>%
  extract_source_tables("original")

no_aspirin_tables <- readRDS(no_aspirin_output_file) %>%
  extract_source_tables("no_aspirin")

all_feature_votes <- bind_rows(
  original_tables$feature_votes,
  no_aspirin_tables$feature_votes
)

all_performance <- bind_rows(
  original_tables$performance,
  no_aspirin_tables$performance
)

# Confirm the planned models are available before analysis.
available_models <- all_feature_votes %>%
  distinct(source, model)

model_availability <- model_registry %>%
  left_join(
    available_models %>% mutate(available = TRUE),
    by = c("source", "model")
  ) %>%
  mutate(available = replace_na(available, FALSE))

save_table(
  model_availability,
  "08A_prespecified_model_availability.csv"
)

missing_models <- model_availability %>%
  filter(!available)

if (nrow(missing_models) > 0) {
  stop(
    "The following prespecified models were not found in the saved Script 03 outputs:\n",
    paste0(
      missing_models$model_label,
      " [", missing_models$model, "]",
      collapse = "\n"
    ),
    "\nInspect 08A_prespecified_model_availability.csv or update model_registry."
  )
}

# Assemble one harmonised feature table for the prespecified model panel.
selected_model_features <- all_feature_votes %>%
  inner_join(
    model_registry,
    by = c("source", "model")
  ) %>%
  mutate(
    feature_modality = classify_feature_modality(feature),
    feature_label = clean_feature_label(feature),
    coefficient_sign = case_when(
      mean_coefficient > 0 ~ 1,
      mean_coefficient < 0 ~ -1,
      TRUE ~ 0
    ),
    stable_3_of_5 = selection_frequency >= stable_selection_cutoff,
    stability_weighted_magnitude =
      selection_frequency * mean_abs_coefficient
  ) %>%
  arrange(display_order, desc(selection_frequency), desc(mean_abs_coefficient))

save_table(
  selected_model_features,
  "08B_all_cross_model_feature_votes.csv"
)

# Attach validation performance so biomarker results can be interpreted alongside
# model discrimination rather than independently of it.
selected_model_performance <- all_performance %>%
  inner_join(
    model_registry,
    by = c("source", "model")
  ) %>%
  arrange(display_order)

save_table(
  selected_model_performance,
  "08C_selected_model_performance.csv"
)

# Build a complete feature-by-model matrix.
# Predictors absent from a model receive selection frequency zero. Their
# coefficient is left missing because zero would incorrectly imply that the
# predictor was offered to the model and penalised to zero in every fold.
feature_model_matrix <- expand_grid(
  feature = unique(selected_model_features$feature),
  model = model_registry$model
) %>%
  left_join(
    model_registry,
    by = "model"
  ) %>%
  left_join(
    selected_model_features %>%
      select(
        source,
        model,
        feature,
        feature_modality,
        feature_label,
        selected_n,
        selection_frequency,
        mean_coefficient,
        median_coefficient,
        mean_abs_coefficient,
        direction,
        coefficient_sign,
        stable_3_of_5,
        stability_weighted_magnitude
      ),
    by = c("source", "model", "feature")
  ) %>%
  group_by(feature) %>%
  mutate(
    feature_modality = first(na.omit(feature_modality)),
    feature_label = first(na.omit(feature_label))
  ) %>%
  ungroup() %>%
  mutate(
    selected_n = replace_na(selected_n, 0),
    selection_frequency = replace_na(selection_frequency, 0),
    stable_3_of_5 = replace_na(stable_3_of_5, FALSE),
    stability_weighted_magnitude =
      replace_na(stability_weighted_magnitude, 0)
  )

save_table(
  feature_model_matrix,
  "08D_complete_feature_by_model_matrix.csv"
)

# Recurrent biomarker summary
#
# No composite consensus score is calculated. The table reports the underlying
# evidence directly: the number of models containing each predictor, the number
# of models in which it was stable, mean fold-selection frequency, coefficient
# direction and direction agreement.
recurrent_summary <- selected_model_features %>%
  group_by(feature, feature_label, feature_modality) %>%
  summarise(
    models_selected_n = n_distinct(model),
    models_stable_n = n_distinct(model[stable_3_of_5]),
    total_outer_fold_selections = sum(selected_n, na.rm = TRUE),
    total_available_outer_folds = 5 * models_selected_n,
    pooled_selection_frequency = safe_divide(
      total_outer_fold_selections,
      total_available_outer_folds
    ),
    mean_model_selection_frequency = mean(selection_frequency, na.rm = TRUE),
    median_model_selection_frequency = median(selection_frequency, na.rm = TRUE),
    mean_coefficient = mean(mean_coefficient, na.rm = TRUE),
    mean_abs_coefficient = mean(mean_abs_coefficient, na.rm = TRUE),
    positive_models_n = sum(mean_coefficient > 0, na.rm = TRUE),
    negative_models_n = sum(mean_coefficient < 0, na.rm = TRUE),
    direction_agreement = safe_divide(
      max(positive_models_n, negative_models_n),
      positive_models_n + negative_models_n
    ),
    recurrent_direction = case_when(
      positive_models_n > negative_models_n ~ "Increase risk",
      negative_models_n > positive_models_n ~ "Decrease risk",
      TRUE ~ "Mixed or neutral"
    ),
    model_labels = paste(model_label[order(display_order)], collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    direction_consistent = direction_agreement == 1,
    stable_in_multiple_models = models_stable_n >= 2
  ) %>%
  arrange(
    desc(models_stable_n),
    desc(models_selected_n),
    desc(mean_model_selection_frequency),
    desc(mean_abs_coefficient)
  ) %>%
  mutate(recurrence_rank = row_number()) %>%
  select(recurrence_rank, everything())

save_table(
  recurrent_summary,
  "08E_cross_model_recurrent_biomarker_summary.csv"
)

recurrent_biomarkers <- recurrent_summary %>%
  filter(models_selected_n >= 2)

save_table(
  recurrent_biomarkers,
  "08F_recurrent_biomarkers_selected_in_multiple_models.csv"
)

direction_conflicts <- recurrent_summary %>%
  filter(models_selected_n >= 2, direction_agreement < 1) %>%
  arrange(direction_agreement, desc(models_selected_n))

save_table(
  direction_conflicts,
  "08G_predictors_with_cross_model_direction_conflict.csv"
)

# Consistent modality ordering across figures
modality_levels <- c(
  "Clinical",
  "Lipoproteins",
  "Lipids",
  "Small molecules",
  "Cytokines",
  "Bile acids",
  "Pathology"
)

# Figure 1: grouped cross-model heatmap
#
# Rows are grouped by modality. Point size represents outer-fold selection
# frequency and point colour represents the mean standardised coefficient.
# Predictors are ordered within modality by cross-model recurrence and stability.
heatmap_features <- recurrent_summary %>%
  filter(models_selected_n >= 2) %>%
  slice_head(n = top_n_heatmap_features) %>%
  pull(feature)

heatmap_order_table <- recurrent_summary %>%
  filter(feature %in% heatmap_features) %>%
  mutate(
    feature_modality = factor(feature_modality, levels = modality_levels)
  ) %>%
  arrange(
    feature_modality,
    models_stable_n,
    models_selected_n,
    mean_model_selection_frequency,
    mean_abs_coefficient
  ) %>%
  mutate(
    display_label = paste0(feature_label, " [", feature_modality, "]")
  )

heatmap_data <- feature_model_matrix %>%
  filter(feature %in% heatmap_features) %>%
  mutate(
    feature_modality = factor(feature_modality, levels = modality_levels),
    model_label = factor(
      model_label,
      levels = c(
        "Expanded clinical\n(aspirin excluded)",
        "Lipoproteins only",
        "Lipids only",
        "Small molecules only",
        "Cytokines only",
        "Bile acids only",
        "Pathology only",
        "Locked candidate",
        "All modalities",
        "Baseline clinical",
        "Expanded clinical"
      )
    ),
    display_label = paste0(feature_label, " [", feature_modality, "]"),
    display_label = factor(
      display_label,
      levels = heatmap_order_table$display_label
    )
  )

cross_model_heatmap <- ggplot(
  heatmap_data,
  aes(x = model_label, y = display_label)
) +
  geom_tile(fill = "grey97", colour = "white", linewidth = 0.4) +
  geom_point(
    data = heatmap_data %>% filter(selection_frequency > 0),
    aes(size = selection_frequency, colour = mean_coefficient),
    alpha = 0.95
  ) +
  facet_grid(
    rows = vars(feature_modality),
    scales = "free_y",
    space = "free_y",
    `switch` = "y"
  ) +
  scale_size_continuous(
    limits = c(0, 1),
    breaks = seq(0.2, 1, 0.2)
  ) +
  scale_colour_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-0.4, 0.4),
    oob = scales::squish,
    na.value = "grey92"
  ) +
  labs(
    title = "Cross-model stability and direction of recurrent predictors",
    subtitle = paste0(
      "Rows are grouped by modality; point size indicates outer-fold selection ",
      "frequency and colour indicates mean standardised coefficient"
    ),
    x = NULL,
    y = NULL,
    size = "Selection\nfrequency",
    colour = "Mean\ncoefficient"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

save_plot(
  cross_model_heatmap,
  "08H_cross_model_biomarker_heatmap_grouped.png",
  width = 12,
  height = max(8, 0.28 * length(heatmap_features) + 4)
)

# Figure 2: recurrent biomarker evidence plot
#
# This replaces the previous composite consensus-score ranking. Predictors are
# ranked directly by stable-model count, total model count and mean selection
# frequency. The x-axis reports mean fold-selection frequency, while point size
# reports the number of models in which the predictor was selected.
recurrent_plot_data <- recurrent_biomarkers %>%
  slice_head(n = top_n_recurrent_features) %>%
  mutate(
    display_label = paste0(feature_label, " [", feature_modality, "]"),
    display_label = factor(
      display_label,
      levels = rev(display_label)
    )
  )

recurrent_evidence_plot <- ggplot(
  recurrent_plot_data,
  aes(
    x = mean_model_selection_frequency,
    y = display_label,
    size = models_selected_n,
    colour = recurrent_direction
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = mean_model_selection_frequency,
      yend = display_label
    ),
    linewidth = 0.8,
    alpha = 0.65
  ) +
  geom_point(alpha = 0.95) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Recurrent biomarkers across model specifications",
    subtitle = paste0(
      "Ranked by stable-model count, model recurrence and mean fold-selection ",
      "frequency; no composite score is used"
    ),
    x = "Mean outer-fold selection frequency",
    y = NULL,
    colour = "Coefficient direction",
    size = "Models\nselected"
  ) +
  theme_minimal(base_size = 11)

save_plot(
  recurrent_evidence_plot,
  "08I_recurrent_biomarker_evidence.png",
  width = 10,
  height = max(7, 0.30 * nrow(recurrent_plot_data) + 2)
)

# Figure 3: directionality in biologically informative models
#
# Clinical-only benchmark panels are excluded here to keep the figure focused on
# the individual omics models, the locked candidate and the all-modality model.
directionality_models <- c(
  "lipoproteins",
  "lipids",
  "small_molecules",
  "cytokines",
  "bile_acids",
  "pathology",
  candidate_model,
  all_modalities_model
)

directionality_features <- recurrent_summary %>%
  filter(models_selected_n >= 2) %>%
  slice_head(n = top_n_recurrent_features) %>%
  pull(feature)

cross_model_directionality_data <- selected_model_features %>%
  filter(
    model %in% directionality_models,
    feature %in% directionality_features,
    selection_frequency >= 0.20
  ) %>%
  mutate(
    display_label = paste0(feature_label, " [", feature_modality, "]"),
    display_label = fct_reorder(display_label, mean_coefficient, .fun = mean),
    model_label = factor(
      model_label,
      levels = model_registry$model_label[
        match(directionality_models, model_registry$model)
      ]
    )
  )

cross_model_directionality_plot <- ggplot(
  cross_model_directionality_data,
  aes(
    x = mean_coefficient,
    y = display_label,
    colour = model_group,
    size = selection_frequency
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_segment(
    aes(x = 0, xend = mean_coefficient, yend = display_label),
    linewidth = 0.6,
    alpha = 0.35
  ) +
  geom_point(alpha = 0.9) +
  facet_wrap(vars(model_label), ncol = 2, scales = "free_y") +
  scale_size_continuous(
    limits = c(0, 1),
    breaks = seq(0.2, 1, 0.2)
  ) +
  labs(
    title = "Direction and stability of recurrent predictors",
    subtitle = paste0(
      "Individual omics models, the locked candidate and the all-modality model; ",
      "positive coefficients increase predicted MACE risk"
    ),
    x = "Mean standardised Elastic Net coefficient",
    y = NULL,
    colour = "Model group",
    size = "Selection\nfrequency"
  ) +
  theme_minimal(base_size = 10)

save_plot(
  cross_model_directionality_plot,
  "08J_recurrent_predictor_directionality_selected_models.png",
  width = 12,
  height = 14
)

# Figure 4: biomarker persistence across model contexts
#
# Models are arranged as a progression from clinical information, through
# individual omics models, to the locked integrated candidate and the maximal
# all-modality model. Gaps between stages make this progression explicit.
persistence_models <- c(
  "clinical_expanded_no_aspirin",
  "lipoproteins",
  "lipids",
  "small_molecules",
  "pathology",
  candidate_model,
  all_modalities_model
)

persistence_positions <- tibble(
  model = persistence_models,
  model_position = c(1, 3, 4, 5, 6, 8, 9),
  model_stage = c(
    "Clinical",
    rep("Individual omics", 4),
    "Integrated candidate",
    "All modalities"
  )
)

persistence_features <- recurrent_summary %>%
  filter(models_selected_n >= 2) %>%
  slice_head(n = min(20, top_n_recurrent_features)) %>%
  pull(feature)

persistence_data <- feature_model_matrix %>%
  filter(
    model %in% persistence_models,
    feature %in% persistence_features
  ) %>%
  left_join(persistence_positions, by = "model") %>%
  left_join(
    recurrent_summary %>% select(feature, recurrence_rank),
    by = "feature"
  ) %>%
  mutate(
    display_label = paste0(feature_label, " [", feature_modality, "]"),
    display_label = fct_reorder(display_label, -recurrence_rank),
    selected = selection_frequency > 0
  )

stage_labels <- persistence_positions %>%
  group_by(model_stage) %>%
  summarise(
    stage_x = mean(range(model_position)),
    .groups = "drop"
  )

biomarker_persistence_plot <- ggplot(
  persistence_data,
  aes(
    x = model_position,
    y = display_label,
    group = feature,
    colour = feature_modality
  )
) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = "grey97", colour = NA) +
  annotate("rect", xmin = 2.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
           fill = "grey98", colour = NA) +
  annotate("rect", xmin = 7.5, xmax = 8.5, ymin = -Inf, ymax = Inf,
           fill = "grey97", colour = NA) +
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = -Inf, ymax = Inf,
           fill = "grey98", colour = NA) +
  geom_line(
    data = persistence_data %>% filter(selected),
    linewidth = 0.6,
    alpha = 0.40
  ) +
  geom_point(
    data = persistence_data %>% filter(selected),
    aes(size = selection_frequency),
    alpha = 0.95
  ) +
  geom_text(
    data = stage_labels,
    aes(x = stage_x, y = Inf, label = model_stage),
    inherit.aes = FALSE,
    vjust = 1.4,
    fontface = "bold",
    size = 3.2
  ) +
  scale_x_continuous(
    breaks = persistence_positions$model_position,
    labels = model_registry$model_label[
      match(persistence_positions$model, model_registry$model)
    ],
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_size_continuous(
    limits = c(0, 1),
    breaks = seq(0.2, 1, 0.2)
  ) +
  labs(
    title = "Persistence of recurrent predictors across model contexts",
    subtitle = paste0(
      "Connected points show where the same predictor remains selected; point ",
      "size indicates outer-fold selection frequency"
    ),
    x = NULL,
    y = NULL,
    colour = "Predictor modality",
    size = "Selection\nfrequency"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 18, r = 10, b = 10, l = 10)
  )

save_plot(
  biomarker_persistence_plot,
  "08K_biomarker_persistence_across_models.png",
  width = 12,
  height = max(8, 0.34 * length(persistence_features) + 3)
)

# Figure 5: model-biomarker flow diagram
#
# This deterministic ggplot diagram connects selected model specifications to
# predictors that were stable in at least 3 of 5 outer folds. Predictors are
# grouped by source modality, and connection opacity represents outer-fold
# selection frequency.
#
# The figure describes model-selection recurrence. It does not imply causal or
# biological interactions between connected variables.

# Colours used consistently for predictor modalities
modality_colours <- c(
  "Clinical" = "#F8766D",
  "Lipoproteins" = "#00BA38",
  "Lipids" = "#7CAE00",
  "Small molecules" = "#C77CFF",
  "Pathology" = "#00B0F6"
)


flow_model_keys <- c(
  "clinical_expanded_no_aspirin",
  "lipoproteins",
  "lipids",
  "small_molecules",
  "pathology",
  candidate_model,
  all_modalities_model
)

# Display order follows the model-development narrative.
flow_model_order <- c(
  "Expanded clinical\n(aspirin excluded)",
  "Lipoproteins only",
  "Lipids only",
  "Small molecules only",
  "Pathology only",
  "Locked candidate",
  "All modalities"
)

# Use recurrent predictors only so the figure remains interpretable.
# The synthesis table retains the complete stable-predictor results.
flow_features <- recurrent_summary %>%
  filter(
    models_stable_n >= 2
  ) %>%
  pull(feature)

flow_membership <- selected_model_features %>%
  filter(
    model %in% flow_model_keys,
    stable_3_of_5,
    feature %in% flow_features
  ) %>%
  distinct(
    model,
    model_label,
    feature,
    feature_label,
    feature_modality,
    selection_frequency
  )

if (nrow(flow_membership) == 0) {
  stop(
    "No stable recurrent predictors were available for the model-biomarker ",
    "flow diagram."
  )
}

# Define modality order so clinically interpretable predictors appear first,
# followed by the molecular modalities.
flow_modality_order <- c(
  "Clinical",
  "Lipoproteins",
  "Lipids",
  "Small molecules",
  "Pathology"
)

# Create model positions on the left.
flow_model_nodes <- flow_membership %>%
  distinct(model, model_label) %>%
  mutate(
    model_label = factor(
      model_label,
      levels = flow_model_order
    )
  ) %>%
  arrange(model_label) %>%
  mutate(
    model_y = rev(seq_len(n())),
    model_x = 0
  ) %>%
  mutate(
    model_label = as.character(model_label)
  )

# Create predictor positions on the right, grouped by modality and then ordered
# by cross-model recurrence and mean selection frequency.
flow_predictor_nodes <- recurrent_summary %>%
  filter(feature %in% unique(flow_membership$feature)) %>%
  mutate(
    feature_modality = factor(
      feature_modality,
      levels = flow_modality_order
    )
  ) %>%
  arrange(
    feature_modality,
    desc(models_stable_n),
    desc(mean_model_selection_frequency),
    feature_label
  ) %>%
  mutate(
    predictor_y = rev(seq_len(n())),
    predictor_x = 1
  ) %>%
  mutate(
    feature_modality = as.character(feature_modality)
  )

# Add visible vertical spacing between modalities.
modality_offsets <- flow_predictor_nodes %>%
  distinct(feature_modality) %>%
  mutate(
    modality_offset = (row_number() - 1) * 1.2
  )

flow_predictor_nodes <- flow_predictor_nodes %>%
  left_join(
    modality_offsets,
    by = "feature_modality"
  ) %>%
  mutate(
    predictor_y = predictor_y - modality_offset
  )

# Join node coordinates onto each model-predictor connection.
flow_edges <- flow_membership %>%
  left_join(
    flow_model_nodes %>%
      select(
        model,
        model_x,
        model_y
      ),
    by = "model"
  ) %>%
  left_join(
    flow_predictor_nodes %>%
      select(
        feature,
        predictor_x,
        predictor_y,
        feature_modality
      ),
    by = c(
      "feature",
      "feature_modality"
    )
  )

# Modality labels are positioned beside each predictor group.
flow_modality_labels <- flow_predictor_nodes %>%
  group_by(feature_modality) %>%
  summarise(
    label_y = max(predictor_y) + 0.55,
    .groups = "drop"
  )

save_table(
  flow_edges,
  "08M_model_biomarker_flow_edges.csv"
)

save_table(
  flow_model_nodes,
  "08N_model_biomarker_flow_model_nodes.csv"
)

save_table(
  flow_predictor_nodes,
  "08O_model_biomarker_flow_predictor_nodes.csv"
)

flow_plot <- ggplot() +
  # Curved connections are less visually severe than diagonal straight lines
  # when several models retain the same predictor.
  geom_curve(
    data = flow_edges,
    aes(
      x = model_x,
      y = model_y,
      xend = predictor_x,
      yend = predictor_y,
      colour = feature_modality,
      alpha = selection_frequency,
      group = interaction(model, feature)
    ),
    curvature = 0.08,
    linewidth = 0.65,
    lineend = "round"
  ) +
  # Model nodes
  geom_point(
    data = flow_model_nodes,
    aes(
      x = model_x,
      y = model_y
    ),
    shape = 22,
    size = 5.5,
    fill = "grey20",
    colour = "grey20"
  ) +
  geom_text(
    data = flow_model_nodes,
    aes(
      x = model_x - 0.035,
      y = model_y,
      label = model_label
    ),
    hjust = 1,
    size = 3.5,
    lineheight = 0.9
  ) +
  # Predictor nodes
  geom_point(
    data = flow_predictor_nodes,
    aes(
      x = predictor_x,
      y = predictor_y,
      colour = feature_modality,
      size = models_stable_n
    ),
    alpha = 0.95
  ) +
  geom_text(
    data = flow_predictor_nodes,
    aes(
      x = predictor_x + 0.035,
      y = predictor_y,
      label = feature_label
    ),
    hjust = 0,
    size = 3.3
  ) +
  # Modality section labels
  geom_text(
    data = flow_modality_labels,
    aes(
      x = 1,
      y = label_y,
      label = feature_modality,
      colour = feature_modality
    ),
    hjust = 0,
    vjust = -0.3,
    fontface = "bold",
    size = 3.7,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = modality_colours,
    breaks = flow_modality_order,
    drop = FALSE
  ) +
  scale_alpha_continuous(
    range = c(0.25, 0.90),
    limits = c(0, 1),
    breaks = c(0.6, 0.8, 1.0)
  ) +
  scale_size_continuous(
    range = c(3, 7),
    breaks = sort(
      unique(flow_predictor_nodes$models_stable_n)
    )
  ) +
  coord_cartesian(
    xlim = c(-0.48, 1.58),
    clip = "off"
  ) +
  labs(
    title = "Stable predictors retained across model specifications",
    subtitle = paste0(
      "Connections identify predictors selected in at least 3 of 5 outer ",
      "folds; colour indicates source modality and opacity indicates ",
      "Outer-fold selection frequency"
    ),
    colour = "Predictor modality",
    alpha = "Selection frequency",
    size = "Stable models"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      size = 11,
      margin = margin(b = 14)
    ),
    legend.position = "right",
    plot.margin = margin(
      t = 15,
      r = 130,
      b = 15,
      l = 170
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    legend.background = element_rect(
      fill = "white",
      colour = NA
    )
  )

save_plot(
  flow_plot,
  "08P_model_biomarker_flow_diagram_final.png",
  width = 14,
  height = 10
)


# Final synthesis table
#
# This table condenses the central result into four model contexts: clinical,
# the relevant individual omics model, the locked candidate and the maximal
# all-modality model. A predictor is marked present only when it was stable in
# at least 3 of 5 outer folds in that context.
clinical_context_model <- "clinical_expanded_no_aspirin"
omics_context_models <- c(
  "lipoproteins",
  "lipids",
  "small_molecules",
  "cytokines",
  "bile_acids",
  "pathology"
)

context_membership <- selected_model_features %>%
  filter(stable_3_of_5) %>%
  mutate(
    context = case_when(
      model == clinical_context_model ~ "Clinical",
      model %in% omics_context_models ~ "Omics-only",
      model == candidate_model ~ "Locked candidate",
      model == all_modalities_model ~ "All modalities",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(context)) %>%
  distinct(feature, context, .keep_all = TRUE)

synthesis_features <- recurrent_summary %>%
  filter(
    models_stable_n >= 2 |
      feature %in% context_membership$feature[
        context_membership$context %in% c("Locked candidate", "All modalities")
      ]
  ) %>%
  pull(feature)

biomarker_synthesis_table <- recurrent_summary %>%
  filter(feature %in% synthesis_features) %>%
  select(
    feature,
    Biomarker = feature_label,
    Modality = feature_modality,
    models_selected_n,
    models_stable_n,
    mean_model_selection_frequency,
    mean_coefficient,
    direction_agreement,
    Direction = recurrent_direction
  ) %>%
  left_join(
    context_membership %>%
      select(feature, context) %>%
      mutate(present = "Yes") %>%
      pivot_wider(
        names_from = context,
        values_from = present,
        values_fill = "No"
      ),
    by = "feature"
  ) %>%
  mutate(
    across(
      any_of(c("Clinical", "Omics-only", "Locked candidate", "All modalities")),
      ~ replace_na(.x, "No")
    ),
    `Models selected, n` = models_selected_n,
    `Stable models, n` = models_stable_n,
    `Mean fold-selection frequency` = round(mean_model_selection_frequency, 2),
    `Mean standardised coefficient` = round(mean_coefficient, 3),
    `Direction agreement, %` = round(100 * direction_agreement, 0)
  ) %>%
  arrange(
    desc(`Stable models, n`),
    desc(`Models selected, n`),
    desc(`Mean fold-selection frequency`)
  ) %>%
  select(
    Biomarker,
    Modality,
    Clinical,
    `Omics-only`,
    `Locked candidate`,
    `All modalities`,
    Direction,
    `Models selected, n`,
    `Stable models, n`,
    `Mean fold-selection frequency`,
    `Mean standardised coefficient`,
    `Direction agreement, %`
  )

save_table(
  biomarker_synthesis_table,
  "08R_key_biomarker_synthesis_table.csv"
)

# Save a compact Markdown copy for direct review and thesis drafting.
markdown_table_lines <- c(
  paste0("| ", paste(names(biomarker_synthesis_table), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(biomarker_synthesis_table)), collapse = " | "), " |"),
  apply(
    biomarker_synthesis_table,
    1,
    function(row_values) {
      paste0("| ", paste(row_values, collapse = " | "), " |")
    }
  )
)

writeLines(
  markdown_table_lines,
  file.path(output_table_dir, "08S_key_biomarker_synthesis_table.md")
)

# Save compact outputs for later thesis figures and supplementary analyses.
recurrence_outputs <- list(
  settings = list(
    original_output_file = original_output_file,
    no_aspirin_output_file = no_aspirin_output_file,
    stable_selection_cutoff = stable_selection_cutoff,
    top_n_heatmap_features = top_n_heatmap_features,
    top_n_recurrent_features = top_n_recurrent_features,
    top_n_flow_features = top_n_flow_features,
    candidate_model = candidate_model,
    all_modalities_model = all_modalities_model
  ),
  model_registry = model_registry,
  model_availability = model_availability,
  selected_model_performance = selected_model_performance,
  selected_model_features = selected_model_features,
  feature_model_matrix = feature_model_matrix,
  recurrent_summary = recurrent_summary,
  recurrent_biomarkers = recurrent_biomarkers,
  direction_conflicts = direction_conflicts,
  persistence_data = persistence_data,
  flow_edges = flow_edges,
  flow_model_nodes = flow_model_nodes,
  flow_feature_nodes = flow_feature_nodes,
  upset_data = upset_data,
  upset_plot_created = upset_plot_created,
  biomarker_synthesis_table = biomarker_synthesis_table
)

saveRDS(
  recurrence_outputs,
  file.path(
    output_model_dir,
    "08_cross_model_biomarker_recurrence_outputs.rds"
  )
)

message("Script 08 cross-model biomarker recurrence analysis complete.")
message("Models compared: ", nrow(model_registry))
message("Unique selected predictors: ", n_distinct(selected_model_features$feature))
message("Predictors selected in at least two models: ", nrow(recurrent_biomarkers))
message("Direction conflicts identified: ", nrow(direction_conflicts))
message("Model-biomarker flow figure created: TRUE")
message("Stable-predictor UpSet plot created: ", upset_plot_created)
message("Synthesis-table rows: ", nrow(biomarker_synthesis_table))
message("Leading recurrent biomarkers:")

print(
  recurrent_summary %>%
    select(
      recurrence_rank,
      feature_label,
      feature_modality,
      models_selected_n,
      models_stable_n,
      mean_model_selection_frequency,
      direction_agreement,
      recurrent_direction
    ) %>%
    slice_head(n = 15)
)
