suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source(file.path("plotting", "R", "figures", "figure2.R"))

extended_figure3_panel_ids <- c("EF3A", "EF3B", "EF3C")

task1_bridge_detail_path <- function() {
  primary <- file.path(dirname(plot_ready_root()), "group_bridge", "task1_internal_vs_cross_group_bridge.csv")
  fallback <- "/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/group_bridge/task1_internal_vs_cross_group_bridge.csv"
  resolve_existing_input_path(primary, fallback)
}

load_task1_bridge_detail <- function() {
  data <- read_m2m_csv(
    task1_bridge_detail_path(),
    required_cols = c(
      "dataset",
      "cell_line",
      "target",
      "representation_detail",
      "metric_family",
      "metric_name",
      "cross_minus_internal"
    )
  )
  data <- data[
    data$representation_detail %in% c("Gene", "Pathway") &
      nzchar(data$cross_minus_internal) &
      nzchar(data$target) &
      data$target != "NA",
    ,
    drop = FALSE
  ]
  data$delta_value <- metric_plot_value(data$metric_name, data$cross_minus_internal)
  data <- data[is.finite(data$delta_value), , drop = FALSE]
  data$family_block <- factor(
    ifelse(data$metric_family == "group_concordance", "Group", "Retrieval"),
    levels = c("Group", "Retrieval")
  )
  data$metric_label <- factor(metric_display_name(data$metric_name), levels = c("Cosine", "Pearson", "eDist", "MRR", "Hit@1", "Hit@5", "Hit@10"))
  data$representation <- factor(data$representation_detail, levels = c("Gene", "Pathway"))
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data
}

bridge_quantile_summary <- function(data, group_cols) {
  q10 <- stats::aggregate(delta_value ~ ., data = data[, c(group_cols, "delta_value"), drop = FALSE], FUN = function(values) stats::quantile(values, probs = 0.1, na.rm = TRUE))
  q25 <- stats::aggregate(delta_value ~ ., data = data[, c(group_cols, "delta_value"), drop = FALSE], FUN = function(values) stats::quantile(values, probs = 0.25, na.rm = TRUE))
  q75 <- stats::aggregate(delta_value ~ ., data = data[, c(group_cols, "delta_value"), drop = FALSE], FUN = function(values) stats::quantile(values, probs = 0.75, na.rm = TRUE))
  q90 <- stats::aggregate(delta_value ~ ., data = data[, c(group_cols, "delta_value"), drop = FALSE], FUN = function(values) stats::quantile(values, probs = 0.9, na.rm = TRUE))
  names(q10)[names(q10) == "delta_value"] <- "q10_delta"
  names(q25)[names(q25) == "delta_value"] <- "q25_delta"
  names(q75)[names(q75) == "delta_value"] <- "q75_delta"
  names(q90)[names(q90) == "delta_value"] <- "q90_delta"
  out <- merge(q10, q25, by = group_cols, all = TRUE)
  out <- merge(out, q75, by = group_cols, all = TRUE)
  merge(out, q90, by = group_cols, all = TRUE)
}

build_bridge_pair_plot <- function(summary, panel_id, x_label) {
  gene <- summary[
    summary$representation == "Gene",
    c("dataset", "family_block", "metric_label", "entity_id", "entity_label", "delta_median"),
    drop = FALSE
  ]
  pathway <- summary[
    summary$representation == "Pathway",
    c("dataset", "family_block", "metric_label", "entity_id", "entity_label", "delta_median"),
    drop = FALSE
  ]
  names(gene)[names(gene) == "delta_median"] <- "gene_delta"
  names(pathway)[names(pathway) == "delta_median"] <- "pathway_delta"
  paired <- merge(
    gene,
    pathway,
    by = c("dataset", "family_block", "metric_label", "entity_id", "entity_label"),
    all = FALSE
  )
  paired$y_pos <- paired$entity_id

  plot <- ggplot(summary, aes(x = delta_median, y = entity_id, colour = representation)) +
    geom_vline(xintercept = 0, linewidth = 0.32, linetype = "dashed", colour = "#D7DEE8") +
    geom_segment(
      data = paired,
      aes(x = gene_delta, xend = pathway_delta, y = y_pos, yend = y_pos),
      inherit.aes = FALSE,
      linewidth = 0.62,
      colour = "#D5DEE8",
      lineend = "round"
    ) +
    geom_segment(
      aes(x = q25_delta, xend = q75_delta, y = entity_id, yend = entity_id),
      linewidth = 1.12,
      alpha = 0.80,
      lineend = "round"
    ) +
    geom_segment(
      aes(x = q10_delta, xend = q90_delta, y = entity_id, yend = entity_id),
      linewidth = 0.28,
      alpha = 0.40
    ) +
    geom_point(aes(fill = representation), shape = 21, size = 2.28, stroke = 0.36, colour = "#1F2D3A") +
    scale_color_representation(drop = FALSE) +
    scale_fill_representation(drop = FALSE) +
    scale_y_discrete(labels = function(values) sub("^.*::", "", values)) +
    facet_grid(dataset ~ family_block + metric_label, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    labs(x = x_label, y = NULL, colour = "Representation", fill = "Representation") +
    theme_m2m() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(size = 6.7, face = "bold"),
      axis.text.y = element_text(size = 6.4, face = "bold"),
      strip.text.x = element_text(size = 7.9, face = "bold"),
      strip.text.y = element_text(angle = 0, size = 7.9),
      panel.spacing.x = grid::unit(5.5, "mm"),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, panel_id)
}

summarise_bridge_by_cell_line <- function() {
  data <- load_task1_bridge_detail()
  group_cols <- c("dataset", "family_block", "metric_label", "representation", "cell_line")
  summary <- stats::aggregate(
    delta_value ~ dataset + family_block + metric_label + representation + cell_line,
    data = data,
    FUN = median
  )
  names(summary)[names(summary) == "delta_value"] <- "delta_median"
  ranges <- bridge_quantile_summary(data, group_cols)

  summary <- merge(
    summary,
    ranges,
    by = group_cols,
    all.x = TRUE
  )

  rank_frame <- stats::aggregate(delta_median ~ dataset + cell_line, data = summary, FUN = mean)
  rank_frame <- rank_frame[order(rank_frame$dataset, -rank_frame$delta_median, rank_frame$cell_line), , drop = FALSE]
  rank_levels <- rev(unique(paste(rank_frame$dataset, rank_frame$cell_line, sep = "::")))
  summary$entity_id <- factor(
    paste(summary$dataset, summary$cell_line, sep = "::"),
    levels = rank_levels
  )
  summary$entity_label <- summary$cell_line
  summary
}

select_bridge_targets <- function(mode = c("conservative", "fragile"), top_n = 12) {
  mode <- match.arg(mode)
  data <- load_task1_bridge_detail()
  target_rank <- stats::aggregate(delta_value ~ dataset + target, data = data, FUN = mean)
  target_rank <- target_rank[order(target_rank$dataset, if (mode == "conservative") -target_rank$delta_value else target_rank$delta_value, target_rank$target), , drop = FALSE]
  selected <- do.call(
    rbind,
    lapply(
      split(target_rank, target_rank$dataset),
      function(frame) head(frame, top_n)
    )
  )

  summary <- stats::aggregate(
    delta_value ~ dataset + family_block + metric_label + representation + target,
    data = data[data$target %in% selected$target & data$dataset %in% selected$dataset, , drop = FALSE],
    FUN = median
  )
  names(summary)[names(summary) == "delta_value"] <- "delta_median"
  target_data <- data[data$target %in% selected$target & data$dataset %in% selected$dataset, , drop = FALSE]
  ranges <- bridge_quantile_summary(target_data, c("dataset", "family_block", "metric_label", "representation", "target"))
  summary <- merge(
    summary,
    ranges,
    by = c("dataset", "family_block", "metric_label", "representation", "target"),
    all.x = TRUE
  )

  target_levels <- rev(unique(paste(selected$dataset, selected$target, sep = "::")))
  summary$entity_id <- factor(
    paste(summary$dataset, summary$target, sep = "::"),
    levels = target_levels
  )
  summary$entity_label <- summary$target
  summary
}

build_panel_EF3A <- function() {
  build_bridge_pair_plot(
    summarise_bridge_by_cell_line(),
    "EF3A",
    "Cross - internal median (positive = more robust)"
  )
}

build_panel_EF3B <- function() {
  build_bridge_pair_plot(
    select_bridge_targets("conservative", top_n = 12),
    "EF3B",
    "Cross - internal median (conservative targets)"
  )
}

build_panel_EF3C <- function() {
  build_bridge_pair_plot(
    select_bridge_targets("fragile", top_n = 12),
    "EF3C",
    "Cross - internal median (fragile targets)"
  )
}

compose_extended_figure3 <- function() {
  panels <- list(build_panel_EF3A(), build_panel_EF3B(), build_panel_EF3C())
  compose_manuscript_figure("extended_figure3", panels, c("EF3A", "EF3B", "EF3C"), ncol = 1)
}

extended_figure3_panel_builders <- function() {
  list(
    `EF3A` = build_panel_EF3A,
    `EF3B` = build_panel_EF3B,
    `EF3C` = build_panel_EF3C
  )
}

extended_figure3_panel_dimensions <- list(
  `EF3A` = list(width_mm = 220, height_mm = 142),
  `EF3B` = list(width_mm = 220, height_mm = 186),
  `EF3C` = list(width_mm = 220, height_mm = 186)
)

render_extended_figure3 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_extended_figure3()
  render_composed_figure(
    "extended_figure3",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

render_extended_figure3_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  render_panel_set(
    figure_id = "extended_figure3",
    panel_ids = extended_figure3_panel_ids,
    panel_builders = extended_figure3_panel_builders(),
    panel_dimensions = extended_figure3_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

if (sys.nframe() == 0) {
  render_extended_figure3()
}
