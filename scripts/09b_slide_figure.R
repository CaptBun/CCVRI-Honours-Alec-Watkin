# Repository stage: slides
# Source supplied as: 09B_top5_auc_slide_graph(2).R
# Run from the project root, preferably through run_pipeline.R.
source("R/bootstrap.R", local = TRUE)

# Presentation-specific forest plot for the five highest-performing models
#
# Purpose:
#   Recreate Figure 09B specifically for a widescreen presentation slide.
#   This script does not fit, tune or compare any models. It reads the final
#   participant-level results summary already produced by:
#
#     09_final_model_and_biomarker_synthesis_v3.R
#
# Required input:
#   outputs/tables/final_thesis_synthesis/
#     09B_top5_auc_comparison_vs_baseline.csv
#
# Outputs:
#   slide_exports/
#     09B_top5_auc_comparison_slide.png
#     09B_top5_auc_comparison_slide.pdf
#     09B_top5_auc_comparison_slide.svg  (when svglite is installed)

library(tidyverse)

# Define file locations
input_file <- file.path(
  "outputs",
  "tables",
  "final_thesis_synthesis",
  "09B_top5_auc_comparison_vs_baseline.csv"
)

# Presentation exports are written to a dedicated folder at the project root.
# Keeping these files separate avoids conflicts with an earlier figure that may
# be open in PowerPoint, Windows Photos or a PDF viewer.
output_directory <- file.path(getwd(), "slide_exports")

if (!file.exists(input_file)) {
  stop(
    "The required results table was not found:\n",
    input_file,
    "\n\nRun 09_final_model_and_biomarker_synthesis_v3.R first."
  )
}

if (!dir.exists(output_directory)) {
  directory_created <- dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = TRUE
  )

  if (!directory_created) {
    stop(
      "R could not create the slide export folder:\n",
      output_directory,
      "\nCheck that the project folder is writable."
    )
  }
}

# Test write access before attempting to open a graphics device. This produces
# a clear error if the project is read-only or controlled by synchronisation or
# security software.
write_test_file <- tempfile(
  pattern = "write_test_",
  tmpdir = output_directory,
  fileext = ".tmp"
)

if (!file.create(write_test_file)) {
  stop(
    "The slide export folder exists but is not writable:\n",
    output_directory,
    "\nChoose a writable project location or change output_directory."
  )
}

unlink(write_test_file)

# Avoid overwriting files that may currently be open in another application.
# The first run uses the clean base filename; later runs receive a timestamp.
make_available_output_path <- function(filename) {
  proposed_path <- file.path(output_directory, filename)

  if (!file.exists(proposed_path)) {
    return(proposed_path)
  }

  extension <- tools::file_ext(filename)
  filename_without_extension <- tools::file_path_sans_ext(filename)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  file.path(
    output_directory,
    paste0(filename_without_extension, "_", timestamp, ".", extension)
  )
}

# Define model names used in the completed analysis
baseline_model <- "clinical_baseline"
expanded_clinical_model <- "clinical_expanded"
candidate_model <- paste0(
  "clinical_expanded_no_aspirin_",
  "lipoproteins_lipids_small_molecules"
)

# Define presentation colours
unsw_yellow <- "#FFE600"
integrated_blue <- "#2F75B5"
expanded_orange <- "#E69F00"
baseline_grey <- "#767676"
confidence_interval_grey <- "#555555"
candidate_band <- "#FFF8C9"
text_colour <- "#202020"

# Convert the internal model identifiers into compact presentation labels
#
# Long model names make the original thesis figure difficult to read once it is
# placed on a slide. The abbreviations are defined in the figure caption so the
# model composition remains explicit without using multi-line paragraphs.
make_short_model_label <- function(model_name, rank) {
  if (model_name == baseline_model) {
    return("Baseline clinical")
  }

  if (model_name == expanded_clinical_model) {
    return("Expanded clinical")
  }

  modality_text <- model_name %>%
    str_remove("^clinical_expanded_no_aspirin_") %>%
    str_remove("^clinical_expanded_no_aspirin$") %>%
    str_replace_all("lipoproteins", "LP") %>%
    str_replace_all("small_molecules", "SM") %>%
    str_replace_all("bile_acids", "BA") %>%
    str_replace_all("lipids", "Lipids") %>%
    str_replace_all("cytokines", "Cytokines") %>%
    str_replace_all("pathology", "Pathology") %>%
    str_replace_all("_", " + ") %>%
    str_squish()

  prefix <- ifelse(is.na(rank), "", paste0(rank, ". "))

  if (model_name == candidate_model) {
    paste0(prefix, "Selected candidate: Expanded + ", modality_text)
  } else {
    paste0(prefix, "Expanded + ", modality_text)
  }
}

# Format p-values consistently while retaining values close to the conventional
# significance threshold. No significance symbols are used because these top-5
# comparisons were exploratory and the raw p-values are more informative.
format_p_value <- function(p_value) {
  case_when(
    is.na(p_value) ~ "Reference",
    p_value < 0.001 ~ "<0.001",
    TRUE ~ formatC(p_value, format = "f", digits = 3)
  )
}

# Read and validate the completed results
plot_data <- read.csv(input_file, check.names = FALSE) %>%
  as_tibble()

required_columns <- c(
  "rank",
  "model",
  "pooled_auc",
  "auc_ci_low",
  "auc_ci_high",
  "delong_p_value"
)

missing_columns <- setdiff(required_columns, names(plot_data))

if (length(missing_columns) > 0) {
  stop(
    "The results table is missing the following required columns:\n",
    paste(missing_columns, collapse = ", ")
  )
}

# Arrange the five integrated models first, followed by the two clinical
# benchmarks. The selected candidate is retained in its empirical AUROC rank.
plot_data <- plot_data %>%
  mutate(
    plot_order = case_when(
      !is.na(rank) ~ rank,
      model == expanded_clinical_model ~ 6,
      model == baseline_model ~ 7,
      TRUE ~ 8
    )
  ) %>%
  arrange(plot_order) %>%
  mutate(
    display_label = map2_chr(model, rank, make_short_model_label),
    display_label = factor(display_label, levels = rev(display_label)),
    is_candidate = model == candidate_model,
    model_class = case_when(
      is_candidate ~ "Selected candidate",
      model == expanded_clinical_model ~ "Expanded clinical",
      model == baseline_model ~ "Baseline clinical",
      TRUE ~ "Other integrated model"
    ),
    auc_ci_text = sprintf(
      "%.3f (%.3f\u2013%.3f)",
      pooled_auc,
      auc_ci_low,
      auc_ci_high
    ),
    p_value_text = format_p_value(delong_p_value)
  )

baseline_auc <- plot_data %>%
  filter(model == baseline_model) %>%
  pull(pooled_auc)

if (length(baseline_auc) != 1) {
  stop("Exactly one baseline clinical row is required in the results table.")
}

# Set graph and results-column positions dynamically from the observed CIs
#
# Only the AUROC region receives tick marks. The two unshaded regions to its
# right function as aligned text columns for the numerical estimate and paired
# DeLong p-value.
axis_min <- floor((min(plot_data$auc_ci_low, na.rm = TRUE) - 0.015) * 20) / 20
axis_min <- max(0.50, axis_min)

axis_max <- 0.95
axis_max <- min(0.95, axis_max)

auc_text_x <- axis_max + 0.035
p_value_x <- axis_max + 0.165
plot_x_max <- axis_max + 0.225

axis_breaks <- seq(axis_min, axis_max, by = 0.05)
candidate_y <- which(levels(plot_data$display_label) ==
  plot_data$display_label[plot_data$is_candidate])

# Create the slide-specific forest plot
slide_auc_plot <- ggplot(
  plot_data,
  aes(x = pooled_auc, y = display_label)
) +
  # A restrained background band makes the selected model immediately visible.
  annotate(
    "rect",
    xmin = axis_min,
    xmax = plot_x_max,
    ymin = candidate_y - 0.43,
    ymax = candidate_y + 0.43,
    fill = candidate_band,
    colour = NA
  ) +
  # The dashed reference line shows the pooled baseline clinical AUROC.
  geom_vline(
    xintercept = baseline_auc,
    linetype = "dashed",
    linewidth = 0.8,
    colour = baseline_grey
  ) +
  # Confidence intervals and point estimates are based only on held-out outer-
  # fold predictions from the completed nested cross-validation analysis.
  geom_segment(
    aes(
      x = auc_ci_low,
      xend = auc_ci_high,
      yend = display_label
    ),
    linewidth = 1.15,
    colour = confidence_interval_grey,
    lineend = "round"
  ) +
  geom_point(
    aes(fill = model_class, size = is_candidate),
    shape = 21,
    stroke = 1.0,
    colour = "white"
  ) +
  # Direct numerical columns remove the need to estimate values from the axis.
  geom_text(
    aes(x = auc_text_x, label = auc_ci_text, fontface = ifelse(is_candidate, "bold", "plain")),
    hjust = 0,
    size = 4.2,
    colour = text_colour,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = p_value_x, label = p_value_text, fontface = ifelse(is_candidate, "bold", "plain")),
    hjust = 0,
    size = 4.2,
    colour = text_colour,
    show.legend = FALSE
  ) +
  # Column headings are placed inside the plot so they remain aligned when the
  # image is resized in PowerPoint.
  annotate(
    "text",
    x = auc_text_x,
    y = nrow(plot_data) + 0.67,
    label = "AUROC (95% CI)",
    hjust = 0,
    fontface = "bold",
    size = 4.2,
    colour = text_colour
  ) +
  annotate(
    "text",
    x = p_value_x,
    y = nrow(plot_data) + 0.67,
    label = "Delong's p\n vs baseline",
    hjust = 0,
    fontface = "bold",
    size = 4.2,
    colour = text_colour
  ) +
  scale_fill_manual(
    values = c(
      "Other integrated model" = integrated_blue,
      "Selected candidate" = unsw_yellow,
      "Expanded clinical" = expanded_orange,
      "Baseline clinical" = baseline_grey
    )
  ) +
  scale_size_manual(values = c(`FALSE` = 4.8, `TRUE` = 6.4)) +
  scale_x_continuous(
    breaks = axis_breaks,
    labels = scales::number_format(accuracy = 0.01),
    limits = c(axis_min, plot_x_max),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.45, 1.05))) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Pooled outer-fold validation AUROC",
    y = NULL,
    caption = paste0(
      "Points show pooled outer-fold AUROC; bars show 95% DeLong confidence intervals.\n",
      "LP, lipoproteins; SM, small molecules; BA, bile acids."
    )
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.margin = margin(t = 16, r = 24, b = 10, l = 12),
    plot.caption = element_text(
      size = 10.5,
      colour = "#555555",
      hjust = 0,
      margin = margin(t = 12)
    ),
    axis.title.x = element_text(
      size = 13,
      colour = "#333333",
      margin = margin(t = 8)
    ),
    axis.text.x = element_text(size = 12, colour = "#444444"),
    axis.text.y = element_text(
      size = 12.5,
      colour = text_colour,
      lineheight = 0.95,
      hjust = 1
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#E5E5E5", linewidth = 0.55),
    legend.position = "none"
  )

# Save a high-resolution PNG for direct placement in PowerPoint and vector
# versions for resizing without loss of clarity.
#
# The base R devices are specified explicitly. This avoids Windows installations
# automatically selecting the ragg or Cairo devices, which can fail when their
# external graphics libraries are unavailable or cannot open the output stream.
png_output_file <- make_available_output_path(
  "09B_top5_auc_comparison_slide.png"
)

pdf_output_file <- make_available_output_path(
  "09B_top5_auc_comparison_slide.pdf"
)

ggsave(
  filename = png_output_file,
  plot = slide_auc_plot,
  width = 13.2,
  height = 6.8,
  units = "in",
  dpi = 320,
  device = grDevices::png,
  bg = "white"
)

ggsave(
  filename = pdf_output_file,
  plot = slide_auc_plot,
  width = 13.2,
  height = 6.8,
  units = "in",
  device = grDevices::pdf,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  svg_output_file <- make_available_output_path(
    "09B_top5_auc_comparison_slide.svg"
  )

  ggsave(
    filename = svg_output_file,
    plot = slide_auc_plot,
    width = 13.2,
    height = 6.8,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
} else {
  message(
    "The PNG and PDF were saved successfully. Install the svglite package ",
    "if an editable SVG is also required."
  )
}

message(
  "Presentation forest plot saved to: ",
  normalizePath(output_directory, winslash = "/", mustWork = FALSE)
)
