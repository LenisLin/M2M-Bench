suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source(file.path("plotting", "R", "figures", "figure2.R"))

extended_figure4_panel_ids <- c("EF4A", "EF4B")

load_shared_enrichment_heatmap <- function(path, entity_column, top_n_per_block) {
  data <- read_plot_ready_csv(
    path,
    required_cols = c(
      "dataset",
      entity_column,
      "entity_display_label",
      "perturbation_type_label",
      "shared_flag",
      "gene_enrichment_score",
      "pathway_enrichment_score",
      "gene_significance_label",
      "pathway_significance_label",
      "pair_mean_enrichment"
    )
  )
  data$shared_flag <- tolower(as.character(data$shared_flag)) == "true"
  data <- data[data$shared_flag, , drop = FALSE]
  if (!nrow(data)) {
    stop(sprintf("No shared rows available for %s.", path), call. = FALSE)
  }

  order_frame <- aggregate(
    abs(pair_mean_enrichment) ~ perturbation_type_label + entity_display_label,
    data = unique(data[, c("perturbation_type_label", "entity_display_label", "pair_mean_enrichment"), drop = FALSE]),
    FUN = max
  )
  names(order_frame)[names(order_frame) == "abs(pair_mean_enrichment)"] <- "rank_score"
  order_frame <- order_frame[order(order_frame$perturbation_type_label, -order_frame$rank_score, order_frame$entity_display_label), , drop = FALSE]
  selected <- do.call(
    rbind,
    lapply(split(order_frame, order_frame$perturbation_type_label), function(frame) head(frame, top_n_per_block))
  )
  data <- merge(
    data,
    selected[, c("perturbation_type_label", "entity_display_label"), drop = FALSE],
    by = c("perturbation_type_label", "entity_display_label"),
    all = FALSE
  )

  long <- rbind(
    data.frame(
      perturbation_type_label = data$perturbation_type_label,
      entity_display_label = data$entity_display_label,
      dataset = data$dataset,
      representation = "Gene",
      enrichment_score = data$gene_enrichment_score,
      significance_label = data$gene_significance_label,
      stringsAsFactors = FALSE
    ),
    data.frame(
      perturbation_type_label = data$perturbation_type_label,
      entity_display_label = data$entity_display_label,
      dataset = data$dataset,
      representation = "Pathway",
      enrichment_score = data$pathway_enrichment_score,
      significance_label = data$pathway_significance_label,
      stringsAsFactors = FALSE
    )
  )
  long$column_label <- factor(
    paste(long$dataset, long$representation, sep = " | "),
    levels = c("LINCS | Gene", "LINCS | Pathway", "scPerturb | Gene", "scPerturb | Pathway")
  )
  row_order <- aggregate(
    abs(enrichment_score) ~ perturbation_type_label + entity_display_label,
    data = long,
    FUN = max
  )
  names(row_order)[names(row_order) == "abs(enrichment_score)"] <- "rank_score"
  row_order <- row_order[order(row_order$perturbation_type_label, row_order$rank_score, row_order$entity_display_label), , drop = FALSE]
  row_order$row_factor <- paste(row_order$perturbation_type_label, row_order$entity_display_label, sep = "::")
  long$row_factor <- factor(
    paste(long$perturbation_type_label, long$entity_display_label, sep = "::"),
    levels = rev(unique(row_order$row_factor))
  )
  long
}

build_shared_enrichment_heatmap <- function(path, entity_column, panel_id, top_n_per_block) {
  data <- load_shared_enrichment_heatmap(path, entity_column, top_n_per_block)
  plot <- ggplot(data, aes(x = column_label, y = row_factor, fill = enrichment_score)) +
    geom_tile(colour = "white", linewidth = 0.28) +
    geom_text(
      data = data[nzchar(data$significance_label), , drop = FALSE],
      aes(label = significance_label),
      size = 1.9,
      colour = "#1F2D3A"
    ) +
    scale_fill_reference_diverging(name = "Log2 enrichment", midpoint = 0) +
    scale_y_discrete(labels = function(values) sub("^.*::", "", values)) +
    facet_grid(perturbation_type_label ~ ., scales = "free_y", space = "free_y") +
    labs(
      x = NULL,
      y = NULL,
      subtitle = "Shared-only heatmap; complete full-surface results remain in the support tables."
    ) +
    theme_m2m() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 20, hjust = 1, size = 6.9, face = "bold"),
      axis.text.y = element_text(size = 6.4),
      strip.text.y = element_text(angle = 0, size = 8.1, face = "bold"),
      plot.subtitle = element_text(size = 8.2, face = "bold", colour = "#526273"),
      legend.position = "bottom"
    )
  finish_panel_plot(plot, panel_id)
}

build_panel_EF4A <- function() {
  build_shared_enrichment_heatmap(
    "extended/extended_figure4_cell_line_high_concordance_full.csv",
    "cell_line",
    "EF4A",
    top_n_per_block = 13
  )
}

build_panel_EF4B <- function() {
  build_shared_enrichment_heatmap(
    "extended/extended_figure4_target_high_concordance_full.csv",
    "target_token",
    "EF4B",
    top_n_per_block = 30
  )
}

compose_extended_figure4 <- function() {
  panels <- list(build_panel_EF4A(), build_panel_EF4B())
  compose_manuscript_figure("extended_figure4", panels, c("EF4A", "EF4B"), ncol = 1)
}

extended_figure4_panel_builders <- function() {
  list(
    `EF4A` = build_panel_EF4A,
    `EF4B` = build_panel_EF4B
  )
}

extended_figure4_panel_dimensions <- list(
  `EF4A` = list(width_mm = 178, height_mm = 154),
  `EF4B` = list(width_mm = 178, height_mm = 188)
)

render_extended_figure4 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_extended_figure4()
  render_composed_figure(
    "extended_figure4",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

render_extended_figure4_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  render_panel_set(
    figure_id = "extended_figure4",
    panel_ids = extended_figure4_panel_ids,
    panel_builders = extended_figure4_panel_builders(),
    panel_dimensions = extended_figure4_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

if (sys.nframe() == 0) {
  render_extended_figure4()
}
