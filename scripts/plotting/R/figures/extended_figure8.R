suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source(file.path("plotting", "R", "figures", "figure3.R"))

extended_figure8_panel_ids <- c("EF8A", "EF8B", "EF8C")

load_contextual_support_full <- function() {
  data <- read_plot_ready_csv("extended/extended_figure8_contextual_support_full.csv")
  data$representation_label <- factor(
    representation_detail_name(data$representation),
    levels = task2_representation_levels[task2_representation_levels %in% unique(representation_detail_name(data$representation))]
  )
  data$metric_label <- factor(
    metric_display_name(data$metric_name),
    levels = c("Cosine", "Pearson", "MRR", "Hit@1", "Hit@5", "Hit@10")
  )
  data$panel_label <- factor(
    paste(data$dataset, metric_display_name(data$metric_name), sep = " | "),
    levels = c(
      paste("LINCS", c("Cosine", "Pearson", "MRR", "Hit@1", "Hit@5", "Hit@10"), sep = " | "),
      paste("scPerturb", c("Cosine", "Pearson", "MRR", "Hit@1", "Hit@5", "Hit@10"), sep = " | ")
    )
  )
  data
}

build_panel_EF8A <- function() {
  data <- load_task2_representation_comparisons()
  data$row_label <- factor(as.character(data$row_label), levels = rev(unique(as.character(data$row_label))))
  data$label_x <- ifelse(data$effect_value >= 0, data$effect_value + 0.015, data$effect_value - 0.015)
  marker_data <- data[nzchar(data$marker), , drop = FALSE]
  marker_data$hjust_value <- ifelse(marker_data$effect_value >= 0, 0, 1)
  plot <- ggplot(data, aes(x = effect_value, y = row_label)) +
    geom_vline(xintercept = 0, linewidth = 0.34, linetype = "dashed", colour = "#D5DEE8") +
    geom_segment(
      aes(x = 0, xend = effect_value, yend = row_label, colour = analysis_family_label),
      linewidth = 1.05,
      lineend = "round",
      alpha = 0.86
    ) +
    geom_point(
      aes(fill = effect_value),
      shape = 21,
      size = 2.7,
      stroke = 0.34,
      colour = "#1F2D3A"
    ) +
    geom_text(
      data = marker_data,
      aes(x = label_x, y = row_label, label = marker),
      inherit.aes = FALSE,
      hjust = marker_data$hjust_value,
      size = 1.95,
      colour = "#6B1830"
    ) +
    scale_color_manual(values = c("Group concordance" = "#5E88C1", "Retrieval" = "#C76898"), drop = FALSE) +
    scale_fill_reference_diverging(name = "Signed median delta", midpoint = 0) +
    coord_cartesian(clip = "off") +
    labs(x = "Signed median delta (Gene vs Pathway)", y = NULL, colour = "Family") +
    theme_m2m() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(size = 6.8, face = "bold"),
      axis.text.y = element_text(size = 6.8),
      legend.position = "bottom"
    )
  finish_panel_plot(plot, "EF8A")
}

build_panel_EF8B <- function() {
  data <- load_task2_fm_local()
  data <- data[representation_detail_name(data$representation) %in% c("Gene", "Pathway", "Geneformer", "scBERT", "scFoundation", "scGPT", "STATE", "Tahoe-X1", "UCE"), , drop = FALSE]
  baseline <- data[data$is_gene_baseline, c("metric_panel", "plot_value", "plot_low", "plot_high"), drop = FALSE]
  names(baseline) <- c("metric_panel", "baseline_value", "baseline_low", "baseline_high")
  compare <- data[!data$is_gene_baseline, , drop = FALSE]
  compare <- merge(compare, baseline, by = "metric_panel", all.x = TRUE)
  compare$row_factor <- factor(
    representation_detail_name(compare$representation),
    levels = rev(c("Pathway", "Geneformer", "scBERT", "scFoundation", "scGPT", "STATE", "Tahoe-X1", "UCE"))
  )
  compare$baseline_y <- as.numeric(compare$row_factor) + 0.14
  compare$model_y <- as.numeric(compare$row_factor) - 0.14

  plot <- ggplot() +
    geom_segment(
      data = compare,
      aes(x = baseline_value, xend = plot_value, y = baseline_y, yend = model_y),
      inherit.aes = FALSE,
      linewidth = 0.72,
      colour = "#D5DEE8",
      lineend = "round"
    ) +
    geom_segment(
      data = compare,
      aes(x = baseline_low, xend = baseline_high, y = baseline_y, yend = baseline_y),
      inherit.aes = FALSE,
      linewidth = 1.0,
      colour = get_m2m_palette()$representation[["Gene"]]
    ) +
    geom_point(
      data = compare,
      aes(x = baseline_value, y = baseline_y),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.4,
      stroke = 0.42,
      fill = get_m2m_palette()$representation[["Gene"]],
      colour = "#1F2D3A"
    ) +
    geom_segment(
      data = compare,
      aes(x = plot_low, xend = plot_high, y = model_y, yend = model_y, colour = representation_label),
      linewidth = 1.0
    ) +
    geom_point(
      data = compare,
      aes(x = plot_value, y = model_y, fill = representation_label),
      shape = 21,
      size = 2.4,
      stroke = 0.42,
      colour = "#1F2D3A"
    ) +
    geom_text(
      data = compare[nzchar(compare$marker), , drop = FALSE],
      aes(x = Inf, y = model_y, label = marker),
      inherit.aes = FALSE,
      hjust = 1.02,
      size = 2.0,
      fontface = "bold",
      colour = "#7A0019"
    ) +
    scale_color_representation(drop = FALSE) +
    scale_fill_representation(drop = FALSE) +
    facet_wrap(~metric_panel, scales = "free_x", nrow = 1) +
    scale_y_continuous(
      breaks = seq_along(levels(compare$row_factor)),
      labels = levels(compare$row_factor),
      expand = expansion(mult = c(0.08, 0.14))
    ) +
    coord_cartesian(clip = "off") +
    labs(x = "Local performance with 95% CI", y = NULL, color = "Model", fill = "Model") +
    theme_m2m() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(size = 7, face = "bold"),
      axis.text.y = element_text(size = 7.1),
      strip.text.x = element_text(size = 8.2, face = "bold"),
      legend.position = "bottom"
    )
  finish_panel_plot(plot, "EF8B")
}

build_panel_EF8C <- function() {
  data <- load_contextual_support_full()
  data <- data[!is.na(data$representation_label), , drop = FALSE]
  plot <- ggplot(
    data,
    aes(x = representation_label, y = metric_value, fill = representation_label, colour = representation_label)
  ) +
    geom_boxplot(
      width = 0.62,
      alpha = 0.18,
      linewidth = 0.36,
      outlier.shape = NA
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.9,
      stroke = 0.32,
      colour = "#1F2D3A"
    ) +
    scale_fill_representation(drop = FALSE) +
    scale_color_representation(drop = FALSE) +
    facet_wrap(~panel_label, scales = "free_x", ncol = 3) +
    labs(x = NULL, y = "Contextual support reference metric", fill = "Representation", colour = "Representation") +
    theme_m2m() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 34, hjust = 1, size = 6.3),
      axis.text.y = element_text(size = 6.7),
      strip.text.x = element_text(size = 8.0, face = "bold"),
      legend.position = "bottom"
    )
  finish_panel_plot(plot, "EF8C")
}

compose_extended_figure8 <- function() {
  panels <- list(build_panel_EF8A(), build_panel_EF8B(), build_panel_EF8C())
  compose_manuscript_figure("extended_figure8", panels, c("EF8A", "EF8B", "EF8C"), ncol = 2)
}

extended_figure8_panel_builders <- function() {
  list(
    `EF8A` = build_panel_EF8A,
    `EF8B` = build_panel_EF8B,
    `EF8C` = build_panel_EF8C
  )
}

extended_figure8_panel_dimensions <- list(
  `EF8A` = list(width_mm = 196, height_mm = 124),
  `EF8B` = list(width_mm = 230, height_mm = 132),
  `EF8C` = list(width_mm = 236, height_mm = 182)
)

render_extended_figure8 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_extended_figure8()
  render_composed_figure(
    "extended_figure8",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

render_extended_figure8_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  render_panel_set(
    figure_id = "extended_figure8",
    panel_ids = extended_figure8_panel_ids,
    panel_builders = extended_figure8_panel_builders(),
    panel_dimensions = extended_figure8_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

if (sys.nframe() == 0) {
  render_extended_figure8()
}
