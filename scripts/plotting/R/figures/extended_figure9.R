suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source(file.path("plotting", "R", "shared", "io.R"))
source(file.path("plotting", "R", "shared", "palette.R"))
source(file.path("plotting", "R", "shared", "labels.R"))
source(file.path("plotting", "R", "shared", "theme.R"))
source(file.path("plotting", "R", "shared", "orderings.R"))
source(file.path("plotting", "R", "shared", "panels.R"))
source(file.path("plotting", "R", "shared", "assembly.R"))
source(file.path("plotting", "R", "shared", "export.R"))

extended_figure9_panel_ids <- c("EF9A", "EF9B")

load_support_vs_suitability <- function() {
  data <- read_plot_ready_csv(
    "extended/extended_figure9_support_vs_suitability.csv",
    required_cols = c(
      "surface_type",
      "dataset",
      "entity_display_label",
      "suitability_score",
      "support_log10",
      "shared_label",
      "rank_within_dataset_surface"
    )
  )
  data$surface_label <- factor(
    ifelse(data$surface_type == "cell_line", "Cell-line surface", "Target surface"),
    levels = c("Cell-line surface", "Target surface")
  )
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data$shared_label <- factor(
    ifelse(data$shared_label == "shared", "Shared", "Dataset-specific"),
    levels = c("Shared", "Dataset-specific")
  )
  data$label_bool <- data$rank_within_dataset_surface <= 3
  data$label_text <- ifelse(data$label_bool, truncate_label(data$entity_display_label, width = 18), "")
  data
}

load_target_exemplars <- function() {
  data <- read_plot_ready_csv(
    "extended/extended_figure9_target_exemplars.csv",
    required_cols = c(
      "selection_rank",
      "selection_label",
      "target",
      "dataset",
      "row_display_label",
      "representation",
      "enrichment_score",
      "support_n"
    )
  )
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data$representation <- factor(data$representation, levels = c("Gene", "Pathway"))
  order_frame <- unique(data[, c("selection_rank", "row_display_label"), drop = FALSE])
  order_frame <- order_frame[order(order_frame$selection_rank), , drop = FALSE]
  data$row_display_label <- factor(data$row_display_label, levels = rev(order_frame$row_display_label))
  data$support_label <- ifelse(
    is.na(data$support_n),
    "",
    paste0("n=", format(data$support_n, big.mark = ","))
  )
  data
}

build_panel_EF9A <- function() {
  data <- load_support_vs_suitability()
  plot <- ggplot(
    data,
    aes(x = support_log10, y = suitability_score, color = dataset, shape = shared_label)
  ) +
    geom_point(size = 2.35, alpha = 0.82) +
    geom_text(
      data = data[data$label_bool, , drop = FALSE],
      aes(label = label_text),
      size = 2.2,
      vjust = -0.55,
      show.legend = FALSE
    ) +
    facet_wrap(~surface_label, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c(LINCS = "#2E5EAA", scPerturb = "#B55D1F"), drop = FALSE) +
    scale_shape_manual(values = c(Shared = 16, `Dataset-specific` = 1), drop = FALSE) +
    labs(
      x = "log10(support size)",
      y = "Suitability summary",
      color = "Dataset",
      shape = "Entity class"
    ) +
    ggplot2::theme(
      strip.text = element_text(size = 8.8, face = "bold"),
      axis.text = element_text(size = 7.2),
      legend.position = "bottom"
    )
  finish_panel_plot(
    plot,
    "EF9A",
    caption_lines = "Diagnostic support plot to check whether high suitability is dominated by support size or persists away from the highest-support strata."
  )
}

build_panel_EF9B <- function() {
  data <- load_target_exemplars()
  gene <- data[data$representation == "Gene", c("dataset", "row_display_label", "enrichment_score"), drop = FALSE]
  pathway <- data[data$representation == "Pathway", c("dataset", "row_display_label", "enrichment_score"), drop = FALSE]
  names(gene)[names(gene) == "enrichment_score"] <- "gene_score"
  names(pathway)[names(pathway) == "enrichment_score"] <- "pathway_score"
  paired <- merge(gene, pathway, by = c("dataset", "row_display_label"), all = FALSE)
  support_ann <- unique(data[data$representation == "Gene", c("dataset", "row_display_label", "support_label"), drop = FALSE])

  plot <- ggplot(
    data,
    aes(x = enrichment_score, y = row_display_label, color = representation)
  ) +
    geom_segment(
      data = paired,
      aes(x = gene_score, xend = pathway_score, y = row_display_label, yend = row_display_label),
      inherit.aes = FALSE,
      linewidth = 0.8,
      color = "#D1D9E2",
      lineend = "round"
    ) +
    geom_point(size = 2.5, alpha = 0.92) +
    geom_text(
      data = support_ann,
      aes(x = Inf, y = row_display_label, label = support_label),
      inherit.aes = FALSE,
      hjust = 1.03,
      size = 2.3,
      color = "#444444"
    ) +
    facet_wrap(~dataset, ncol = 2, scales = "free_y") +
    coord_cartesian(clip = "off") +
    scale_color_representation(drop = FALSE) +
    labs(
      x = "Target suitability enrichment score",
      y = NULL,
      color = "Representation"
    ) +
    ggplot2::theme(
      strip.text = element_text(size = 8.8, face = "bold"),
      axis.text.x = element_text(size = 7.2),
      axis.text.y = element_text(size = 6.4),
      legend.position = "bottom"
    )
  finish_panel_plot(
    plot,
    "EF9B",
    caption_lines = "Target archetypes chosen to illustrate shared anchors, cross-dataset gaps, and dataset-specific high-suitability examples without turning the supplement into a causal mechanism claim."
  )
}

compose_extended_figure9 <- function() {
  panels <- list(build_panel_EF9A(), build_panel_EF9B())
  compose_manuscript_figure("extended_figure9", panels, c("EF9A", "EF9B"), ncol = 1)
}

extended_figure9_panel_builders <- function() {
  list(
    `EF9A` = build_panel_EF9A,
    `EF9B` = build_panel_EF9B
  )
}

extended_figure9_panel_dimensions <- list(
  `EF9A` = list(width_mm = 178, height_mm = 110),
  `EF9B` = list(width_mm = 178, height_mm = 135)
)

render_extended_figure9 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_extended_figure9()
  render_composed_figure(
    "extended_figure9",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

render_extended_figure9_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  render_panel_set(
    figure_id = "extended_figure9",
    panel_ids = extended_figure9_panel_ids,
    panel_builders = extended_figure9_panel_builders(),
    panel_dimensions = extended_figure9_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

if (sys.nframe() == 0) {
  render_extended_figure9()
}
