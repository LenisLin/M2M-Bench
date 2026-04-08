suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source(file.path("plotting", "R", "shared", "io.R"))
source(file.path("plotting", "R", "shared", "palette.R"))
source(file.path("plotting", "R", "shared", "labels.R"))
source(file.path("plotting", "R", "shared", "theme.R"))
source(file.path("plotting", "R", "shared", "orderings.R"))
source(file.path("plotting", "R", "shared", "panels.R"))
source(file.path("plotting", "R", "shared", "assembly.R"))
source(file.path("plotting", "R", "shared", "export.R"))

figure3_panel_ids <- c("3B", "3C", "3D", "3E", "3F")

task2_metric_panel <- function(analysis_family, metric_name) {
  prefix <- ifelse(analysis_family == "group_concordance", "Group", "C2G retrieval")
  paste(prefix, metric_display_name(metric_name), sep = " | ")
}

task2_metric_panel_levels <- c(
  "Group | Cosine",
  "Group | Pearson",
  "C2G retrieval | MRR",
  "C2G retrieval | Hit@1",
  "C2G retrieval | Hit@5",
  "C2G retrieval | Hit@10"
)

task2_metric_scoreboard_levels <- c(
  "Group::Cosine",
  "Group::Pearson",
  "Retrieval::MRR",
  "Retrieval::Hit@1",
  "Retrieval::Hit@5",
  "Retrieval::Hit@10"
)

task2_retrieval_metric_levels <- c("MRR", "Hit@1", "Hit@5", "Hit@10")

task2_specificity_tier_levels <- c(
  "The Cleanest Hits",
  "The Family Hits",
  "The Promiscuous Hits"
)

load_task2_performance_overview <- function() {
  data <- read_plot_ready_csv(
    "figure3/figure3_panel_3b_c2g_performance_overview.csv",
    required_cols = c(
      "analysis_family",
      "dataset",
      "cell_line",
      "direction",
      "representation",
      "metric_name",
      "metric_value"
    )
  )
  if (!("specificity_tier" %in% names(data))) {
    data$specificity_tier <- NA_character_
  }
  if (!("n_queries" %in% names(data)) && "n_total" %in% names(data)) {
    data$n_queries <- suppressWarnings(as.numeric(data$n_total))
  }
  if (!("n_units" %in% names(data)) && "n_targets_metric_valid" %in% names(data)) {
    data$n_units <- suppressWarnings(as.numeric(data$n_targets_metric_valid))
  }
  data <- data[data$representation %in% c("Gene", "Pathway"), , drop = FALSE]
  data$representation <- factor(data$representation, levels = c("Gene", "Pathway"))
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data$family_block <- factor(
    ifelse(data$analysis_family == "group_concordance", "Group", "Retrieval"),
    levels = c("Group", "Retrieval")
  )
  data$metric_panel <- factor(task2_metric_panel(data$analysis_family, data$metric_name), levels = task2_metric_panel_levels)
  data$metric_label <- factor(metric_display_name(data$metric_name), levels = c("Cosine", "Pearson", "MRR", "Hit@1", "Hit@5", "Hit@10"))
  data$specificity_tier <- trimws(as.character(data$specificity_tier))
  order_parts <- lapply(c("LINCS", "scPerturb"), function(ds) {
    sub <- data[data$dataset == ds, , drop = FALSE]
    ds_order <- aggregate(metric_value ~ cell_line, data = sub, FUN = mean)
    ds_order <- ds_order[order(-ds_order$metric_value, ds_order$cell_line), , drop = FALSE]
    paste(ds, ds_order$cell_line, sep = "::")
  })
  row_levels <- rev(unlist(order_parts))
  data$row_factor <- factor(
    paste(data$dataset, data$cell_line, sep = "::"),
    levels = row_levels
  )
  data$cell_line_label <- data$cell_line
  data$display_label <- metric_value_label(data$metric_name, data$metric_value, digits = 2)
  data
}

load_task2_suitability_panel <- function(path, entity_column) {
  data <- read_plot_ready_csv(
    path,
    required_cols = c(
      "dataset",
      entity_column,
      "row_id",
      "global_row_order",
      "row_rank",
      "entity_display_label",
      "suitability_score",
      "gene_enrichment_score",
      "pathway_enrichment_score",
      "gene_significance_label",
      "pathway_significance_label",
      "support_value",
      "support_n",
      "support_label",
      "pair_mean_enrichment",
      "shared_across_modalities_bool"
    )
  )
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data$shared_across_modalities_bool <- tolower(as.character(data$shared_across_modalities_bool)) == "true"
  row_keys <- paste(data$dataset, data$row_id, sep = "::")
  ordered_keys <- unique(row_keys[order(data$dataset, data$global_row_order, data$row_id)])
  data$row_factor <- factor(
    row_keys,
    levels = rev(ordered_keys)
  )
  data$facet_value <- data$dataset
  data$entity_label <- data$entity_display_label
  data$left_annotation_primary <- ifelse(data$shared_across_modalities_bool, "Shared", "")
  data$left_annotation_secondary <- paste0("#", data$row_rank)
  data
}

load_task2_cell_line_patterns <- function() {
  data <- load_task2_suitability_panel(
    "figure3/figure3_panel_3c_cell_line_pattern.csv",
    "cell_line"
  )
  data
}

load_task2_target_patterns <- function() {
  data <- load_task2_suitability_panel(
    "figure3/figure3_panel_3d_target_pattern_summary.csv",
    "target"
  )
  data
}

load_task2_representation_comparisons <- function() {
  data <- read_plot_ready_csv(
    "figure3/figure3_panel_3e_gene_vs_pathway_paired.csv",
    required_cols = c(
      "row_kind",
      "dataset",
      "analysis_family",
      "metric_name",
      "bh_q",
      "test_status",
      "effect_direction",
      "median_delta",
      "n_test_units"
    )
  )
  data <- unique(
    data[
      data$row_kind == "summary",
      c("dataset", "analysis_family", "metric_name", "bh_q", "test_status", "effect_direction", "median_delta", "n_test_units"),
      drop = FALSE
    ]
  )
  data$analysis_family_label <- factor(
    analysis_family_display_name(data$analysis_family),
    levels = c("Group concordance", "Retrieval")
  )
  data$row_label <- factor(
    paste(data$dataset, metric_display_name(data$metric_name), sep = " | "),
    levels = rev(unique(paste(data$dataset, metric_display_name(data$metric_name), sep = " | ")))
  )
  data$effect_value <- signed_effect_value(data$effect_direction, data$median_delta)
  data$marker <- comparison_stats_label(data$bh_q, data$test_status, data$n_test_units, data$median_delta)
  data
}

load_task2_fm_local <- function() {
  data <- read_plot_ready_csv(
    "figure3/figure3_panel_3f_fm_local_tradeoff.csv",
    required_cols = c(
      "representation",
      "analysis_family",
      "metric_name",
      "metric_panel",
      "metric_value",
      "gene_reference_metric_value",
      "gene_reference_plot_value",
      "ci_low",
      "ci_high",
      "n_units",
      "bh_q_vs_gene",
      "test_status_vs_gene",
      "is_gene_baseline"
    )
  )
  data$is_gene_baseline <- tolower(as.character(data$is_gene_baseline)) == "true"
  data$representation_group <- factor(representation_class_name(data$representation), levels = c("Gene", "Pathway", "FM"))
  data$representation_label <- factor(
    representation_detail_name(data$representation),
    levels = task2_representation_levels[task2_representation_levels %in% representation_detail_name(data$representation)]
  )
  data$fm_model_label <- representation_detail_name(data$representation)
  data$baseline_representation <- factor(data$representation_group, levels = c("Gene", "Pathway", "FM"))
  data$analysis_family_label <- factor(
    analysis_family_display_name(data$analysis_family),
    levels = c("Group concordance", "Retrieval")
  )
  data$surface_anchor_label <- factor(
    paste(data$dataset, data$cell_line, sep = " | "),
    levels = unique(paste(data$dataset, data$cell_line, sep = " | "))
  )
  data$metric_label <- factor(metric_display_name(data$metric_name), levels = c("Cosine", "Pearson", "MRR", "Hit@1", "Hit@5", "Hit@10"))
  data$metric_panel <- factor(data$metric_panel, levels = task2_metric_panel_levels)
  data$plot_value <- metric_plot_value(data$metric_name, data$metric_value)
  data$plot_low <- metric_plot_value(data$metric_name, data$ci_low)
  data$plot_high <- metric_plot_value(data$metric_name, data$ci_high)
  data$gene_reference_plot_value <- metric_plot_value(data$metric_name, data$gene_reference_metric_value)
  data$marker <- ifelse(data$test_status_vs_gene == "tested", significance_marker(data$bh_q_vs_gene), "")
  data$point_size <- ifelse(data$is_gene_baseline, 3.4, 2.7)
  data
}

load_task2_lincs_common_board <- function() {
  data <- load_task2_performance_overview()
  data <- data[data$dataset == "LINCS", , drop = FALSE]
  cell_order <- aggregate(metric_value ~ cell_line, data = data, FUN = mean)
  cell_order <- cell_order[order(-cell_order$metric_value, cell_order$cell_line), , drop = FALSE]
  row_levels <- rev(paste("LINCS", cell_order$cell_line, sep = "::"))
  data$row_factor <- factor(paste(data$dataset, data$cell_line, sep = "::"), levels = row_levels)
  data
}

build_lincs_common_scope_scoreboard <- function() {
  data <- load_task2_lincs_common_board()
  data$metric_family_label <- factor(
    paste(as.character(data$family_block), as.character(data$metric_label), sep = "::"),
    levels = task2_metric_scoreboard_levels
  )
  data$column_index <- as.numeric(data$metric_family_label)
  data$row_index <- as.numeric(data$row_factor)

  layout <- unique(
    data[, c("dataset", "family_block", "metric_family_label", "column_index", "row_factor", "row_index"), drop = FALSE]
  )
  layout$xmin <- layout$column_index - 0.46
  layout$xmax <- layout$column_index + 0.46
  layout$ymin <- layout$row_index - 0.44
  layout$ymax <- layout$row_index + 0.44

  panel_max <- ave(
    suppressWarnings(as.numeric(data$metric_value)),
    interaction(data$dataset, data$metric_family_label, drop = TRUE),
    FUN = function(values) {
      values <- values[is.finite(values)]
      if (!length(values)) {
        return(1)
      }
      max(values, na.rm = TRUE)
    }
  )
  panel_max[!is.finite(panel_max) | panel_max <= 0] <- 1
  data$bar_fraction <- pmax(0, pmin(1, suppressWarnings(as.numeric(data$metric_value)) / panel_max))
  data$bar_xmin <- data$column_index - 0.38
  data$bar_xmax <- data$bar_xmin + 0.72 * data$bar_fraction
  data$bar_y <- data$row_index + ifelse(data$representation == "Gene", 0.15, -0.15)
  data$point_x <- data$bar_xmax

  x_breaks <- seq_along(task2_metric_scoreboard_levels)
  x_labels <- gsub("^.*::", "", task2_metric_scoreboard_levels)
  y_breaks <- seq_along(levels(data$row_factor))
  y_labels <- gsub("^.*::", "", levels(data$row_factor))

  plot <- ggplot() +
    geom_rect(
      data = layout,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "#F8FAFC",
      colour = "white",
      linewidth = 0.28
    ) +
    geom_segment(
      data = data,
      aes(x = bar_xmin, xend = bar_xmax, y = bar_y, yend = bar_y, colour = representation),
      linewidth = 2.65,
      alpha = 0.96,
      lineend = "butt"
    ) +
    geom_point(
      data = data,
      aes(x = point_x, y = bar_y, fill = representation),
      shape = 21,
      size = 1.72,
      stroke = 0.32,
      colour = "#1F2D3A"
    ) +
    scale_color_representation(drop = FALSE) +
    scale_fill_representation(drop = FALSE) +
    facet_grid(. ~ family_block, scales = "free", space = "free") +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = NULL,
      colour = "Representation",
      fill = "Representation",
      subtitle = "LINCS common benchmark | rows = cell lines | lanes = Gene and Pathway"
    ) +
    theme_m2m() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold", size = 6.8),
      axis.text.y = element_text(face = "bold", size = 6.6),
      strip.text.x = element_text(size = 8.1, face = "bold"),
      plot.subtitle = element_text(size = 8.3, face = "bold", colour = "#526273"),
      legend.position = "bottom"
    )

  plot
}

prepare_task2_catalog_fraction <- function(plot_data, group_cols) {
  plot_data$group_id <- interaction(plot_data[group_cols], drop = TRUE, lex.order = TRUE)
  min_map <- tapply(plot_data$metric_value, plot_data$group_id, min, na.rm = TRUE)
  max_map <- tapply(plot_data$metric_value, plot_data$group_id, max, na.rm = TRUE)
  plot_data$value_min <- unname(min_map[as.character(plot_data$group_id)])
  plot_data$value_max <- unname(max_map[as.character(plot_data$group_id)])
  missing_group <- !is.finite(plot_data$value_min) | !is.finite(plot_data$value_max)
  plot_data$value_span <- ifelse(
    missing_group,
    NA_real_,
    pmax(plot_data$value_max - plot_data$value_min, 1e-6)
  )
  plot_data$bar_fraction <- ifelse(
    missing_group | is.na(plot_data$metric_value),
    NA_real_,
    ifelse(
      plot_data$value_max <= plot_data$value_min,
      1,
      pmax(0, pmin(1, (plot_data$metric_value - plot_data$value_min) / plot_data$value_span))
    )
  )
  plot_data
}

present_task2_specificity_tier_levels <- function(values) {
  value_text <- as.character(values)
  present <- unique(value_text[!is.na(value_text) & nzchar(value_text)])
  unexpected <- setdiff(present, task2_specificity_tier_levels)
  if (length(unexpected) > 0) {
    stop(
      sprintf(
        "Unexpected specificity_tier values in Figure 3B/EF5 input: %s",
        paste(sort(unexpected), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  task2_specificity_tier_levels[task2_specificity_tier_levels %in% present]
}

load_task2_scperturb_tier_board <- function() {
  data <- load_task2_performance_overview()
  data <- data[
    data$dataset == "scPerturb" &
      data$analysis_family == "retrieval" &
      data$cell_line == "K562" &
      data$direction == "C2G" &
      nzchar(data$specificity_tier),
    ,
    drop = FALSE
  ]
  if (nrow(data) == 0) {
    stop(
      "Figure 3B/EF5 requires scPerturb K562 retrieval rows with specificity_tier in figure3_panel_3b_c2g_performance_overview.csv.",
      call. = FALSE
    )
  }
  tier_levels <- present_task2_specificity_tier_levels(data$specificity_tier)
  data$specificity_tier <- factor(data$specificity_tier, levels = tier_levels)
  data$representation_label <- factor(
    representation_detail_name(data$representation),
    levels = c("Gene", "Pathway")
  )
  data$metric_label <- factor(metric_display_name(data$metric_name), levels = task2_retrieval_metric_levels)
  data$display_label <- metric_value_label(data$metric_name, data$metric_value, digits = 2)
  data$n_queries_label <- ifelse(
    is.na(data$n_queries),
    "",
    paste0("n=", format(round(suppressWarnings(as.numeric(data$n_queries))), big.mark = ","))
  )
  data
}

build_scperturb_tier_scoreboard <- function() {
  data <- load_task2_scperturb_tier_board()
  tier_levels <- present_task2_specificity_tier_levels(data$specificity_tier)
  plot_data <- unique(
    data[
      ,
      c("specificity_tier", "representation", "metric_label", "metric_value", "display_label"),
      drop = FALSE
    ]
  )
  plot_data$row_factor <- factor(
    as.character(plot_data$specificity_tier),
    levels = rev(tier_levels)
  )
  plot_data$metric_id <- as.numeric(plot_data$metric_label)
  plot_data$row_id <- as.numeric(plot_data$row_factor)
  panel_max <- ave(
    suppressWarnings(as.numeric(plot_data$metric_value)),
    plot_data$metric_label,
    FUN = function(values) {
      values <- values[is.finite(values)]
      if (!length(values)) {
        return(1)
      }
      max(values, na.rm = TRUE)
    }
  )
  panel_max[!is.finite(panel_max) | panel_max <= 0] <- 1
  plot_data$bar_fraction <- pmax(0, pmin(1, suppressWarnings(as.numeric(plot_data$metric_value)) / panel_max))
  plot_data$xmin <- plot_data$metric_id - 0.42
  plot_data$xmax <- plot_data$metric_id + 0.42
  plot_data$ymin <- plot_data$row_id - 0.34
  plot_data$ymax <- plot_data$row_id + 0.34
  plot_data$bar_xmin <- plot_data$xmin + 0.08
  plot_data$bar_xmax <- plot_data$bar_xmin + 0.68 * plot_data$bar_fraction
  plot_data$bar_y <- plot_data$row_id + ifelse(as.character(plot_data$representation) == "Gene", 0.14, -0.14)
  plot_data$point_x <- plot_data$bar_xmax
  background <- unique(plot_data[, c("metric_id", "row_id", "xmin", "xmax", "ymin", "ymax"), drop = FALSE])

  ggplot() +
    geom_rect(
      data = background,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#F8FAFC",
      color = "#E3E8EE",
      linewidth = 0.18
    ) +
    geom_segment(
      data = plot_data,
      aes(x = bar_xmin, xend = bar_xmax, y = bar_y, yend = bar_y, colour = representation),
      linewidth = 2.55,
      lineend = "butt"
    ) +
    geom_point(
      data = plot_data,
      aes(x = point_x, y = bar_y, fill = representation),
      shape = 21,
      size = 1.8,
      stroke = 0.32,
      colour = "#1F2D3A"
    ) +
    geom_text(
      data = plot_data,
      aes(x = metric_id, y = row_id + 0.19, label = display_label),
      size = 1.55,
      color = "#1F2D3A",
      fontface = "bold"
    ) +
    scale_color_representation(drop = FALSE) +
    scale_fill_representation(drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(task2_retrieval_metric_levels),
      labels = task2_retrieval_metric_levels,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(plot_data$row_factor)),
      labels = levels(plot_data$row_factor),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    labs(
      x = NULL,
      y = NULL,
      colour = "Representation",
      fill = "Representation",
      subtitle = "scPerturb K562 tier-local retrieval | rows = specificity tiers | lanes = Gene and Pathway"
    ) +
    theme_m2m() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6.8, face = "bold"),
      axis.text.y = element_text(size = 7.0, face = "bold"),
      plot.subtitle = element_text(size = 8.3, face = "bold", colour = "#526273"),
      legend.position = "bottom"
    )
}

build_scperturb_tier_detail_scoreboard <- function(tier_value, panel_id) {
  data <- load_task2_scperturb_tier_board()
  plot_data <- data[data$specificity_tier == tier_value, , drop = FALSE]
  if (nrow(plot_data) == 0) {
    return(
      panel_text_card(
        panel_id,
        "tier scoreboard missing",
        body_lines = c(
          sprintf("No scPerturb/K562 retrieval rows were found for %s.", tier_value),
          "Check figure3_panel_3b_c2g_performance_overview.csv and the Drug_meta.csv tier join."
        )
      )
    )
  }
  plot_data <- unique(
    plot_data[, c("representation_label", "metric_label", "metric_value", "display_label"), drop = FALSE]
  )
  plot_data <- prepare_task2_catalog_fraction(plot_data, c("metric_label"))
  plot_data$metric_id <- as.numeric(plot_data$metric_label)
  plot_data$row_id <- as.numeric(plot_data$representation_label)
  plot_data$xmin <- plot_data$metric_id - 0.42
  plot_data$xmax <- plot_data$metric_id + 0.42
  plot_data$ymin <- plot_data$row_id - 0.34
  plot_data$ymax <- plot_data$row_id + 0.34
  plot_data$bar_xmin <- plot_data$xmin + 0.08
  plot_data$bar_xmax <- plot_data$bar_xmin + 0.68 * plot_data$bar_fraction
  background <- unique(plot_data[, c("metric_id", "row_id", "xmin", "xmax", "ymin", "ymax"), drop = FALSE])

  ggplot() +
    geom_rect(
      data = background,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#F8FAFC",
      color = "#E3E8EE",
      linewidth = 0.18
    ) +
    geom_segment(
      data = plot_data,
      aes(x = bar_xmin, xend = bar_xmax, y = row_id, yend = row_id, colour = representation_label),
      linewidth = 2.7,
      lineend = "butt"
    ) +
    geom_text(
      data = plot_data,
      aes(x = metric_id, y = row_id + 0.18, label = display_label),
      size = 1.95,
      color = "#1F2D3A",
      fontface = "bold"
    ) +
    scale_color_representation(drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(levels(plot_data$metric_label)),
      labels = levels(plot_data$metric_label),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(plot_data$representation_label)),
      labels = levels(plot_data$representation_label),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    labs(
      x = NULL,
      y = NULL,
      colour = "Representation",
      subtitle = paste("scPerturb K562 |", tier_value)
    ) +
    theme_m2m() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6.8, face = "bold"),
      axis.text.y = element_text(size = 7.0),
      plot.subtitle = element_text(size = 8.3, face = "bold", colour = "#526273"),
      legend.position = "bottom"
    )
}

build_panel_3B <- function() {
  upper <- build_lincs_common_scope_scoreboard()
  lower <- build_scperturb_tier_scoreboard()
  combined <- (upper / lower) +
    patchwork::plot_layout(heights = c(1.2, 1), guides = "collect") &
    theme(legend.position = "bottom")
  finish_panel_plot(combined, "3B")
}

build_panel_3C <- function() {
  data <- load_task2_cell_line_patterns()
  paired_enrichment_exemplar_plot(
    data,
    panel_id = "3C",
    x_label = "Suitability summary (paired enrichment score)",
    annotation_primary = "left_annotation_primary",
    annotation_secondary = "left_annotation_secondary"
  )
}

build_panel_3D <- function() {
  data <- load_task2_target_patterns()
  paired_enrichment_exemplar_plot(
    data,
    panel_id = "3D",
    x_label = "Suitability summary (paired enrichment score)",
    annotation_primary = "left_annotation_primary",
    annotation_secondary = "left_annotation_secondary"
  )
}

build_panel_3E <- function() {
  data <- read_plot_ready_csv(
    "figure3/figure3_panel_3e_gene_vs_pathway_paired.csv",
    required_cols = c(
      "row_kind",
      "dataset",
      "analysis_family",
      "direction",
      "metric_name",
      "representation",
      "metric_value",
      "mean_value",
      "ci_low",
      "ci_high",
      "unit_id",
      "bh_q",
      "test_status",
      "n_test_units"
    )
  )
  data$representation <- factor(data$representation, levels = c("Gene", "Pathway"))
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data$metric_panel <- factor(task2_metric_panel(data$analysis_family, data$metric_name), levels = task2_metric_panel_levels)

  raw <- data[data$row_kind == "raw", , drop = FALSE]
  raw$plot_value <- metric_plot_value(raw$metric_name, raw$metric_value)
  raw_lincs <- raw[raw$dataset == "LINCS", , drop = FALSE]
  raw_scperturb <- raw[raw$dataset == "scPerturb", , drop = FALSE]

  gene <- raw[raw$representation == "Gene", c("dataset", "analysis_family", "direction", "metric_name", "metric_panel", "unit_id", "plot_value"), drop = FALSE]
  pathway <- raw[raw$representation == "Pathway", c("dataset", "analysis_family", "direction", "metric_name", "metric_panel", "unit_id", "plot_value"), drop = FALSE]
  names(gene)[names(gene) == "plot_value"] <- "gene_value"
  names(pathway)[names(pathway) == "plot_value"] <- "pathway_value"
  paired <- merge(gene, pathway, by = c("dataset", "analysis_family", "direction", "metric_name", "metric_panel", "unit_id"), all = FALSE)
  paired$gene_representation <- factor("Gene", levels = c("Gene", "Pathway"))
  paired$pathway_representation <- factor("Pathway", levels = c("Gene", "Pathway"))
  paired_lincs <- paired[paired$dataset == "LINCS", , drop = FALSE]
  paired_scperturb <- paired[paired$dataset == "scPerturb", , drop = FALSE]

  summary <- unique(
    data[data$row_kind == "summary", c("dataset", "analysis_family", "direction", "metric_name", "representation", "metric_value", "mean_value", "ci_low", "ci_high", "n_test_units", "bh_q", "test_status", "median_delta"), drop = FALSE]
  )
  summary$plot_value <- metric_plot_value(summary$metric_name, summary$metric_value)
  summary$plot_low <- metric_plot_value(summary$metric_name, summary$ci_low)
  summary$plot_high <- metric_plot_value(summary$metric_name, summary$ci_high)

  stats <- unique(
    summary[, c("dataset", "analysis_family", "direction", "metric_name", "bh_q", "test_status", "n_test_units", "median_delta"), drop = FALSE]
  )
  stats$dataset <- factor(stats$dataset, levels = c("LINCS", "scPerturb"))
  stats$metric_panel <- factor(task2_metric_panel(stats$analysis_family, stats$metric_name), levels = task2_metric_panel_levels)
  stats$marker <- comparison_stats_label(stats$bh_q, stats$test_status, stats$n_test_units, stats$median_delta)
  stats <- stats[nzchar(stats$marker), , drop = FALSE]

  plot <- ggplot(summary, aes(x = representation, y = plot_value, color = representation)) +
    geom_segment(
      data = paired_lincs,
      aes(
        x = gene_representation,
        xend = pathway_representation,
        y = gene_value,
        yend = pathway_value,
        group = unit_id
      ),
      inherit.aes = FALSE,
      linewidth = 0.10,
      color = "#A8B4C0",
      alpha = 0.018
    ) +
    geom_segment(
      data = paired_scperturb,
      aes(
        x = gene_representation,
        xend = pathway_representation,
        y = gene_value,
        yend = pathway_value,
        group = unit_id
      ),
      inherit.aes = FALSE,
      linewidth = 0.24,
      color = "#A8B4C0",
      alpha = 0.12
    ) +
    geom_point(
      data = raw_lincs,
      position = position_jitter(width = 0.045, height = 0),
      size = 0.38,
      alpha = 0.022
    ) +
    geom_point(
      data = raw_scperturb,
      position = position_jitter(width = 0.08, height = 0),
      size = 0.92,
      alpha = 0.18
    ) +
    geom_linerange(
      aes(ymin = plot_low, ymax = plot_high),
      linewidth = 0.98,
      alpha = 0.84
    ) +
    geom_point(
      aes(fill = representation),
      shape = 21,
      size = 3.15,
      stroke = 0.45,
      color = "#1F2D3A"
    ) +
    geom_text(
      data = stats,
      aes(x = 1.5, y = Inf, label = marker),
      inherit.aes = FALSE,
      vjust = 1.40,
      size = 1.72,
      fontface = "bold",
      color = "#A61E2B"
    ) +
    scale_color_representation(drop = FALSE) +
    scale_fill_representation(drop = FALSE) +
    facet_grid(dataset ~ metric_panel, scales = "free_y") +
    labs(x = NULL, y = "Performance-only Gene vs Pathway comparison", color = "Representation", fill = "Representation") +
    theme_m2m() +
    theme(
      panel.grid.major.x = element_blank(),
      strip.text.y = element_text(angle = 0),
      panel.spacing.y = unit(8.5, "mm"),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, "3E")
}

build_panel_3F <- function() {
  data <- load_task2_fm_local()
  reference_lines <- unique(data[, c("metric_panel", "gene_reference_plot_value"), drop = FALSE])
  plot <- ggplot(data, aes(x = plot_value, y = representation_label, color = representation_group)) +
    geom_vline(
      data = reference_lines,
      aes(xintercept = gene_reference_plot_value),
      linetype = "dashed",
      linewidth = 0.4,
      color = "#7A0019",
      alpha = 0.55
    ) +
    geom_segment(
      aes(x = plot_low, xend = plot_high, yend = representation_label),
      linewidth = 0.6,
      alpha = 0.55
    ) +
    geom_point(aes(shape = representation_group), size = 2.7, stroke = 0.45, fill = "white") +
    geom_text(
      data = data[nzchar(data$marker), , drop = FALSE],
      aes(x = Inf, label = marker),
      hjust = 1.04, size = 2.0, fontface = "bold", color = "#7A0019",
      show.legend = FALSE
    ) +
    scale_color_representation(drop = FALSE) +
    scale_shape_manual(values = c(Gene = 16, Pathway = 17, FM = 21), drop = FALSE) +
    facet_wrap(~metric_panel, scales = "free_x", ncol = 3) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Absolute local performance",
      y = NULL,
      color = "Method class",
      shape = "Method class",
      subtitle = "scPerturb | K562 only"
    ) +
    theme_m2m() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 7.1),
      plot.subtitle = element_text(size = 9, face = "bold", color = "#526273"),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, "3F")
}

compose_figure3 <- function() {
  panels <- list(
    build_panel_3B(),
    build_panel_3C(),
    build_panel_3D(),
    build_panel_3E(),
    build_panel_3F()
  )
  compose_manuscript_figure("figure3", panels, figure3_panel_ids, ncol = 2)
}

figure3_panel_builders <- function() {
  list(
    `3B` = build_panel_3B,
    `3C` = build_panel_3C,
    `3D` = build_panel_3D,
    `3E` = build_panel_3E,
    `3F` = build_panel_3F
  )
}

figure3_panel_dimensions <- list(
  `3B` = list(width_mm = 220, height_mm = 190),
  `3C` = list(width_mm = 200, height_mm = 128),
  `3D` = list(width_mm = 208, height_mm = 144),
  `3E` = list(width_mm = 200, height_mm = 130),
  `3F` = list(width_mm = 200, height_mm = 126)
)

render_figure3 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_figure3()
  render_composed_figure(
    "figure3",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export,
    panel_ids = figure3_panel_ids,
    remaining_placeholder_flags = character(),
    all_real_panels_included = TRUE
  )
}

render_figure3_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE,
  separate_legends = FALSE
) {
  render_panel_set(
    figure_id = "figure3",
    panel_ids = figure3_panel_ids,
    panel_builders = figure3_panel_builders(),
    panel_dimensions = figure3_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export,
    separate_legends = separate_legends
  )
}

if (sys.nframe() == 0) {
  render_figure3()
}
