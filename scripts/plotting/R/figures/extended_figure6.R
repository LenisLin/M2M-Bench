suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source(file.path("plotting", "R", "figures", "figure3.R"))

extended_figure6_panel_ids <- c("EF6A", "EF6B")

task2_matrix_column_levels <- c(
  "Cosine | Gene",
  "Cosine | Pathway",
  "Pearson | Gene",
  "Pearson | Pathway",
  "MRR | Gene",
  "MRR | Pathway",
  "Hit@1 | Gene",
  "Hit@1 | Pathway",
  "Hit@5 | Gene",
  "Hit@5 | Pathway",
  "Hit@10 | Gene",
  "Hit@10 | Pathway"
)

load_cell_line_pattern_matrix <- function() {
  data <- read_analysis_csv(
    "figure3_task2_cell_line_pattern_summary.csv",
    required_cols = c(
      "dataset",
      "cell_line",
      "analysis_family",
      "direction",
      "representation",
      "metric_name",
      "n_targets",
      "n_queries_total",
      "mean_metric_value"
    )
  )
  data <- data[
    data$dataset == "LINCS" &
      data$representation %in% c("Gene", "Pathway") &
      data$analysis_family %in% c("group_concordance", "retrieval") &
      (data$analysis_family != "retrieval" | is.na(data$direction) | data$direction == "C2G") &
      !grepl("edist", data$metric_name, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  data$metric_label <- metric_display_name(data$metric_name)
  data$column_label <- factor(
    paste(data$metric_label, data$representation, sep = " | "),
    levels = task2_matrix_column_levels
  )
  data <- percent_rank_within(
    data,
    group_cols = c("column_label"),
    value_col = "mean_metric_value",
    out_col = "fill_value"
  )
  order_frame <- aggregate(mean_metric_value ~ cell_line, data = data, FUN = mean)
  order_frame <- order_frame[order(-order_frame$mean_metric_value, order_frame$cell_line), , drop = FALSE]
  data$cell_line <- factor(data$cell_line, levels = rev(order_frame$cell_line))
  data
}

load_target_pattern_matrix <- function(top_n_per_dataset = 24) {
  data <- read_analysis_csv(
    "figure3_task2_target_pattern_summary.csv",
    required_cols = c(
      "dataset",
      "target",
      "analysis_family",
      "direction",
      "representation",
      "metric_name",
      "n_cell_lines",
      "n_queries_total",
      "mean_metric_value"
    )
  )
  data <- data[
    data$representation %in% c("Gene", "Pathway") &
      data$analysis_family %in% c("group_concordance", "retrieval") &
      (data$analysis_family != "retrieval" | data$direction == "C2G") &
      !grepl("edist", data$metric_name, ignore.case = TRUE),
    ,
    drop = FALSE
  ]
  score_frame <- aggregate(mean_metric_value ~ dataset + target, data = data, FUN = mean)
  score_frame <- score_frame[order(score_frame$dataset, -score_frame$mean_metric_value, score_frame$target), , drop = FALSE]
  keep <- do.call(
    rbind,
    lapply(split(score_frame, score_frame$dataset), function(frame) head(frame, top_n_per_dataset))
  )
  data <- merge(data, keep[, c("dataset", "target"), drop = FALSE], by = c("dataset", "target"), all = FALSE)
  data$metric_label <- metric_display_name(data$metric_name)
  data$column_label <- factor(
    paste(data$metric_label, data$representation, sep = " | "),
    levels = task2_matrix_column_levels
  )
  data <- percent_rank_within(
    data,
    group_cols = c("dataset", "column_label"),
    value_col = "mean_metric_value",
    out_col = "fill_value"
  )
  order_frame <- aggregate(mean_metric_value ~ dataset + target, data = data, FUN = mean)
  order_frame <- order_frame[order(order_frame$dataset, -order_frame$mean_metric_value, order_frame$target), , drop = FALSE]
  order_frame$row_factor <- paste(order_frame$dataset, order_frame$target, sep = "::")
  data$row_factor <- factor(
    paste(data$dataset, data$target, sep = "::"),
    levels = rev(unique(order_frame$row_factor))
  )
  data
}

build_task2_matrix_heatmap <- function(data, panel_id, row_column, facet_formula = NULL, subtitle = NULL) {
  plot <- ggplot(data, aes(x = column_label, y = .data[[row_column]], fill = fill_value)) +
    geom_tile(colour = "white", linewidth = 0.26) +
    scale_x_discrete(labels = function(values) gsub(" \\| ", "\n", values)) +
    scale_fill_gradientn(
      colours = c("#F7F8FA", "#DCE8F2", "#9DB9D8", "#4E79A7"),
      values = c(0, 0.35, 0.7, 1),
      guide = compact_colorbar_guide(),
      name = "Within-column percentile"
    ) +
    labs(
      x = NULL,
      y = NULL,
      subtitle = subtitle
    ) +
    theme_m2m() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 5.7, lineheight = 0.92),
      axis.text.y = element_text(size = 6.3),
      plot.subtitle = element_text(size = 8.2, face = "bold", colour = "#526273"),
      legend.position = "bottom"
    )
  if (!is.null(facet_formula)) {
    plot <- plot + facet_grid(facet_formula, scales = "free_y", space = "free_y")
  }
  finish_panel_plot(plot, panel_id)
}

build_panel_EF6A <- function() {
  build_task2_matrix_heatmap(
    load_cell_line_pattern_matrix(),
    panel_id = "EF6A",
    row_column = "cell_line",
    subtitle = "LINCS full cell-line surface with one row per entity."
  )
}

build_panel_EF6B <- function() {
  data <- load_target_pattern_matrix()
  plot <- ggplot(data, aes(x = column_label, y = row_factor, fill = fill_value)) +
    geom_tile(colour = "white", linewidth = 0.26) +
    scale_x_discrete(labels = function(values) gsub(" \\| ", "\n", values)) +
    scale_fill_gradientn(
      colours = c("#F7F8FA", "#DCE8F2", "#9DB9D8", "#4E79A7"),
      values = c(0, 0.35, 0.7, 1),
      guide = compact_colorbar_guide(),
      name = "Within-column percentile"
    ) +
    scale_y_discrete(
      labels = function(values) sub("^.*::", "", values),
      expand = expansion(add = c(0.25, 1.10))
    ) +
    facet_grid(dataset ~ ., scales = "free_y", space = "free_y") +
    labs(
      x = NULL,
      y = NULL,
      subtitle = "Target-detail matrix with repeated y-axis rows collapsed to one row per target."
    ) +
    theme_m2m() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 5.7, lineheight = 0.92),
      axis.text.y = element_text(size = 6.0),
      strip.text.y = element_text(angle = 0, size = 8.0, face = "bold"),
      plot.subtitle = element_text(size = 8.2, face = "bold", colour = "#526273"),
      legend.position = "bottom"
    )
  finish_panel_plot(plot, "EF6B") +
    theme(
      plot.tag.position = c(-0.03, 1.02),
      plot.margin = margin(14, 10, 7, 16)
    )
}

compose_extended_figure6 <- function() {
  panels <- list(build_panel_EF6A(), build_panel_EF6B())
  compose_manuscript_figure("extended_figure6", panels, c("EF6A", "EF6B"), ncol = 1)
}

extended_figure6_panel_builders <- function() {
  list(
    `EF6A` = build_panel_EF6A,
    `EF6B` = build_panel_EF6B
  )
}

extended_figure6_panel_dimensions <- list(
  `EF6A` = list(width_mm = 178, height_mm = 150),
  `EF6B` = list(width_mm = 178, height_mm = 176)
)

render_extended_figure6 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_extended_figure6()
  render_composed_figure(
    "extended_figure6",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

render_extended_figure6_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  render_panel_set(
    figure_id = "extended_figure6",
    panel_ids = extended_figure6_panel_ids,
    panel_builders = extended_figure6_panel_builders(),
    panel_dimensions = extended_figure6_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

if (sys.nframe() == 0) {
  render_extended_figure6()
}
