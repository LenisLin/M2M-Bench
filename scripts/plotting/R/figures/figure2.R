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

figure2_panel_ids <- c("2A", "2B", "2C", "2D", "2E", "2F")

task1_representation_levels <- c(
  "Gene",
  "Pathway",
  "Geneformer",
  "scBERT",
  "scFoundation",
  "scGPT",
  "STATE",
  "Tahoe-X1",
  "Tahoe-X1 3B",
  "UCE",
  "FM PCA50",
  "FM PCA100",
  "FM PCA200"
)

task1_metric_panel <- function(analysis_family, metric_name) {
  family_prefix <- ifelse(analysis_family == "group_concordance", "Group", "Retrieval")
  paste(family_prefix, metric_display_name(metric_name), sep = " | ")
}

task1_metric_panel_levels <- c(
  "Group | Cosine",
  "Group | Pearson",
  "Group | eDist",
  "Retrieval | MRR",
  "Retrieval | Hit@1",
  "Retrieval | Hit@5",
  "Retrieval | Hit@10"
)

task1_metric_label_levels <- c("Cosine", "Pearson", "eDist", "MRR", "Hit@1", "Hit@5", "Hit@10")

task1_context_palette <- function() {
  c(
    "LINCS | Chemical" = "#4F81BD",
    "LINCS | Genetic" = "#7AA6D1",
    "scPerturb | Chemical" = "#C97B63",
    "scPerturb | Genetic" = "#8E7DBE"
  )
}

annotation_strip_palette <- function() {
  c(
    "Internal" = "#DCE7F3",
    "Cross" = "#F3E2D7",
    "LINCS" = "#4F81BD",
    "scPerturb" = "#C97B63",
    "LINCS -> scPerturb" = "#6D8FB5",
    "scPerturb -> LINCS" = "#A7717A",
    "Chemical" = "#E2EEF7",
    "Genetic" = "#E8DEEF"
  )
}

order_task1_representations <- function(values) {
  labels <- representation_detail_name(values)
  factor(labels, levels = rev(task1_representation_levels[task1_representation_levels %in% labels]))
}

task1_representation_group <- function(values) {
  labels <- as.character(values)
  ifelse(labels == "Gene", "Gene", ifelse(labels == "Pathway", "Pathway", "FM"))
}

pair_summary_offsets <- function(representation) {
  ifelse(as.character(representation) == "Gene", -0.22, 0.22)
}

pair_context_offsets <- function(group_id) {
  ave(
    seq_along(group_id),
    group_id,
    FUN = function(index) {
      if (length(index) == 1) {
        return(0)
      }
      seq(-0.12, 0.12, length.out = length(index))
    }
  )
}

build_pair_segments <- function(raw_data, context_column) {
  gene <- raw_data[raw_data$representation == "Gene", c(context_column, "analysis_family", "metric_name", "metric_panel", "unit_id", "plot_value"), drop = FALSE]
  pathway <- raw_data[raw_data$representation == "Pathway", c(context_column, "analysis_family", "metric_name", "metric_panel", "unit_id", "plot_value"), drop = FALSE]
  names(gene)[names(gene) == "plot_value"] <- "gene_value"
  names(pathway)[names(pathway) == "plot_value"] <- "pathway_value"
  wide <- merge(
    gene,
    pathway,
    by = c(context_column, "analysis_family", "metric_name", "metric_panel", "unit_id"),
    all = FALSE
  )
  wide$context_factor <- factor(wide[[context_column]], levels = levels(raw_data[[context_column]]))
  wide$y_base <- as.numeric(wide$context_factor)
  wide$y_offset <- pair_context_offsets(interaction(wide$metric_panel, wide$context_factor, drop = TRUE))
  wide$y_pos <- wide$y_base + wide$y_offset
  wide
}

pair_points_long <- function(wide_data) {
  gene_points <- wide_data[, c("metric_panel", "unit_id", "context_factor", "y_pos", "gene_value"), drop = FALSE]
  names(gene_points)[names(gene_points) == "gene_value"] <- "plot_value"
  gene_points$representation <- "Gene"
  pathway_points <- wide_data[, c("metric_panel", "unit_id", "context_factor", "y_pos", "pathway_value"), drop = FALSE]
  names(pathway_points)[names(pathway_points) == "pathway_value"] <- "plot_value"
  pathway_points$representation <- "Pathway"
  rbind(gene_points, pathway_points)
}

pair_plot <- function(data, panel_id, context_column, x_label) {
  summary <- data[data$row_kind == "summary", , drop = FALSE]
  summary$context_factor <- factor(summary[[context_column]], levels = levels(data[[context_column]]))
  summary$y_pos <- as.numeric(summary$context_factor) + pair_summary_offsets(summary$representation) * 0.65

  ann <- unique(summary[, c(context_column, "metric_panel", "bh_q", "test_status", "n_test_units")])
  ann$marker <- comparison_annotation_label(ann$bh_q, ann$test_status, ann$n_test_units)
  ann$context_factor <- factor(ann[[context_column]], levels = levels(data[[context_column]]))
  ann$y_pos <- as.numeric(ann$context_factor) + 0.34
  ann$x_pos <- Inf
  ann <- ann[nzchar(ann$marker), , drop = FALSE]

  plot <- ggplot() +
    geom_segment(
      data = summary,
      aes(x = plot_q10, xend = plot_q90, y = y_pos, yend = y_pos, color = representation),
      linewidth = 0.45,
      alpha = 0.80
    ) +
    geom_segment(
      data = summary,
      aes(x = plot_q25, xend = plot_q75, y = y_pos, yend = y_pos, color = representation),
      linewidth = 2.0,
      alpha = 0.90,
      lineend = "round"
    ) +
    geom_point(
      data = summary,
      aes(x = plot_q50, y = y_pos, fill = representation),
      shape = 21,
      size = 2.8,
      stroke = 0.4,
      color = "#1F2D3A"
    ) +
    geom_text(
      data = ann,
      aes(x = x_pos, y = y_pos, label = marker),
      inherit.aes = FALSE,
      hjust = 1.04,
      size = 2.25,
      fontface = "bold",
      color = "#A61E2B"
    ) +
    scale_color_representation(drop = FALSE) +
    scale_fill_representation(drop = FALSE) +
    facet_wrap(~metric_panel, scales = "free_x", ncol = 3) +
    scale_y_continuous(
      breaks = seq_along(levels(data[[context_column]])),
      labels = levels(data[[context_column]]),
      expand = expansion(mult = c(0.08, 0.12))
    ) +
    coord_cartesian(clip = "off") +
    labs(x = x_label, y = NULL, color = "Representation", fill = "Representation") +
    theme_m2m() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(face = "bold"),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, panel_id)
}

load_task1_scope <- function() {
  data <- read_plot_ready_csv(
    "figure2/figure2_panel_2a_task1_scope.csv",
    required_cols = c(
      "scope",
      "analysis_family",
      "slice_label",
      "dataset_or_direction",
      "perturbation_type",
      "representation_class",
      "scope_status",
      "coverage_denominator",
      "coverage_fraction",
      "count_annotation",
      "count_annotation_type"
    )
  )
  data$scope_block <- factor(scope_display_name(data$scope), levels = c("Internal", "Cross"))
  data$representation_class <- factor(data$representation_class, levels = c("Gene", "Pathway", "FM", "ALL"))
  display_order <- c(
    "LINCS | Chemical",
    "LINCS | Genetic",
    "scPerturb | Chemical",
    "scPerturb | Genetic",
    "LINCS_to_scPerturb | Genetic",
    "scPerturb_to_LINCS | Genetic",
    "Cross eligibility | Chemical",
    "Cross eligibility | Genetic"
  )
  data$slice_label <- factor(data$slice_label, levels = rev(display_order[display_order %in% unique(as.character(data$slice_label))]))
  data$count_label <- ifelse(
    data$count_annotation_type == "n_matched_keys",
    paste0("matched=", format(round(data$count_annotation), big.mark = ",")),
    paste0("n=", format(round(data$count_annotation), big.mark = ","))
  )
  data
}

load_task1_internal_performance <- function() {
  data <- read_plot_ready_csv(
    "figure2/figure2_panel_2b_internal_performance_overview.csv",
    required_cols = c(
      "analysis_family",
      "metric_name",
      "dataset_or_direction",
      "perturbation_type",
      "representation",
      "value"
    )
  )
  data$representation_label <- order_task1_representations(data$representation)
  data$dataset_label <- factor(dataset_display_name(data$dataset_or_direction), levels = c("LINCS", "scPerturb"))
  data$context_label <- factor(data$perturbation_type, levels = c("Chemical", "Genetic"))
  data$metric_panel <- factor(task1_metric_panel(data$analysis_family, data$metric_name), levels = task1_metric_panel_levels)
  data$display_label <- metric_value_label(data$metric_name, data$value, digits = 2)
  data
}

load_task1_cross_performance <- function() {
  data <- load_task1_merged_pairwise()
  data <- data[
    data$row_kind == "summary" &
      data$comparison_context %in% c("LINCS -> scPerturb", "scPerturb -> LINCS"),
    ,
    drop = FALSE
  ]
  data$fill_value <- percent_rank_within(
    data,
    group_cols = c("analysis_family", "metric_name"),
    value_col = "plot_value",
    out_col = "fill_value"
  )$fill_value
  data$context_label <- factor(
    data$comparison_context,
    levels = c("LINCS -> scPerturb", "scPerturb -> LINCS")
  )
  data$representation_label <- order_task1_representations(data$representation)
  data$value <- data$plot_value
  data
}

load_task1_performance <- function(mode = "internal") {
  data <- switch(
    mode,
    internal = load_task1_internal_performance(),
    cross = load_task1_cross_performance(),
    stop(sprintf("Unsupported Task 1 performance mode: %s", mode), call. = FALSE)
  )
  data$metric_label <- factor(
    metric_display_name(data$metric_name),
    levels = task1_metric_label_levels
  )
  data$analysis_family_label <- factor(
    analysis_family_display_name(data$analysis_family),
    levels = c("Group concordance", "Retrieval")
  )
  data
}

load_task1_merged_pairwise <- function() {
  data <- read_plot_ready_csv(
    "figure2/figure2_panel_2c_gene_vs_pathway_matched_units.csv",
    required_cols = c(
      "row_kind",
      "comparison_context",
      "analysis_family",
      "metric_name",
      "representation",
      "metric_value",
      "mean_value",
      "summary_mean",
      "summary_median",
      "q10_value",
      "q25_value",
      "q50_value",
      "q75_value",
      "q90_value",
      "ci_low",
      "ci_high",
      "unit_id",
      "bh_q",
      "test_status"
    )
  )
  data$comparison_context <- factor(
    data$comparison_context,
    levels = c("Internal LINCS", "Internal scPerturb", "LINCS -> scPerturb", "scPerturb -> LINCS")
  )
  data$representation <- factor(data$representation, levels = c("Gene", "Pathway"))
  data$metric_panel <- factor(task1_metric_panel(data$analysis_family, data$metric_name), levels = task1_metric_panel_levels)
  data$plot_value <- ifelse(
    data$row_kind == "summary",
    metric_plot_value(data$metric_name, data$q50_value),
    metric_plot_value(data$metric_name, data$metric_value)
  )
  data$plot_q10 <- metric_plot_value(data$metric_name, data$q10_value)
  data$plot_q25 <- metric_plot_value(data$metric_name, data$q25_value)
  data$plot_q50 <- metric_plot_value(data$metric_name, data$q50_value)
  data$plot_q75 <- metric_plot_value(data$metric_name, data$q75_value)
  data$plot_q90 <- metric_plot_value(data$metric_name, data$q90_value)
  data$plot_low <- metric_plot_value(data$metric_name, data$ci_low)
  data$plot_high <- metric_plot_value(data$metric_name, data$ci_high)
  data
}

load_task1_degradation <- function() {
  data <- read_plot_ready_csv(
    "figure2/figure2_panel_2d_internal_to_cross_degradation.csv",
    required_cols = c(
      "dataset",
      "family_block",
      "analysis_family",
      "metric_name",
      "representation",
      "scope",
      "value_transform",
      "n_units",
      "summary_mean",
      "summary_median",
      "q10_value",
      "q25_value",
      "q50_value",
      "q75_value",
      "q90_value",
      "bh_q",
      "test_status",
      "effect_direction",
      "median_delta",
      "n_test_units"
    )
  )
  data$representation <- factor(data$representation, levels = c("Gene", "Pathway"))
  data$family_block <- factor(data$family_block, levels = c("Group", "Retrieval"))
  data$metric_panel <- factor(task1_metric_panel(data$analysis_family, data$metric_name), levels = task1_metric_panel_levels)
  data$metric_display <- factor(metric_display_name(data$metric_name), levels = task1_metric_label_levels)
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data$scope <- factor(data$scope, levels = c("Internal", "Cross"))
  data$plot_median <- metric_plot_value(data$metric_name, data$q50_value)
  data$plot_q50 <- data$plot_median
  data$plot_q10 <- metric_plot_value(data$metric_name, data$q10_value)
  data$plot_q25 <- metric_plot_value(data$metric_name, data$q25_value)
  data$plot_q75 <- metric_plot_value(data$metric_name, data$q75_value)
  data$plot_q90 <- metric_plot_value(data$metric_name, data$q90_value)
  data
}

build_task1_internal_to_cross_degradation_plot <- function(data, panel_id, y_label) {
  data$representation <- factor(data$representation, levels = c("Gene", "Pathway"))
  data$family_block <- factor(data$family_block, levels = c("Group", "Retrieval"))
  data$metric_panel <- factor(task1_metric_panel(data$analysis_family, data$metric_name), levels = task1_metric_panel_levels)
  data$metric_display <- factor(metric_display_name(data$metric_name), levels = task1_metric_label_levels)
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data$scope <- factor(data$scope, levels = c("Internal", "Cross"))
  data$scope_group <- factor(
    paste(data$representation, data$scope, sep = "::"),
    levels = c(
      "Gene::Internal",
      "Gene::Cross",
      "Pathway::Internal",
      "Pathway::Cross"
    )
  )
  if (!("plot_q50" %in% colnames(data)) && "plot_median" %in% colnames(data)) {
    data$plot_q50 <- data$plot_median
  }

  panel_max <- stats::aggregate(plot_q90 ~ dataset + family_block + metric_display, data = data, FUN = max)
  panel_min <- stats::aggregate(plot_q10 ~ dataset + family_block + metric_display, data = data, FUN = min)
  stats <- unique(
    data[
      ,
      c("dataset", "family_block", "metric_display", "bh_q", "test_status", "n_test_units", "median_delta"),
      drop = FALSE
    ]
  )
  stats$inline_label <- comparison_stats_inline_label(
    stats$bh_q,
    test_status = stats$test_status,
    n_test_units = stats$n_test_units,
    median_delta = stats$median_delta
  )
  stats <- stats[nzchar(stats$inline_label), , drop = FALSE]
  stats <- merge(stats, panel_max, by = c("dataset", "family_block", "metric_display"), all.x = TRUE)
  stats <- merge(stats, panel_min, by = c("dataset", "family_block", "metric_display"), all.x = TRUE)
  stats$y_pad <- pmax(0.05, 0.10 * (stats$plot_q90 - stats$plot_q10))
  stats$y_pos <- stats$plot_q90 + stats$y_pad

  dodge <- position_dodge(width = 0.72)

  plot <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.32, colour = "#DCE3EB") +
    geom_linerange(
      data = data,
      aes(
        x = metric_display,
        ymin = plot_q10,
        ymax = plot_q90,
        colour = representation,
        linetype = scope,
        group = scope_group
      ),
      position = dodge,
      linewidth = 0.30,
      alpha = 0.34
    ) +
    geom_linerange(
      data = data,
      aes(
        x = metric_display,
        ymin = plot_q25,
        ymax = plot_q75,
        colour = representation,
        linetype = scope,
        group = scope_group
      ),
      position = dodge,
      linewidth = 0.72,
      alpha = 0.74
    ) +
    geom_point(
      data = data,
      aes(
        x = metric_display,
        y = plot_q50,
        fill = representation,
        shape = scope,
        group = scope_group
      ),
      position = dodge,
      size = 2.45,
      stroke = 0.36,
      colour = "#1F2D3A"
    ) +
    geom_text(
      data = stats,
      aes(
        x = metric_display,
        y = y_pos,
        label = inline_label
      ),
      inherit.aes = FALSE,
      size = 1.56,
      lineheight = 0.9,
      colour = "#6B1830",
      show.legend = FALSE
    ) +
    scale_color_representation(drop = FALSE) +
    scale_fill_representation(drop = FALSE) +
    scale_linetype_manual(values = c(Internal = "solid", Cross = "22"), drop = FALSE) +
    scale_shape_manual(values = c(Internal = 21, Cross = 24), drop = FALSE) +
    facet_grid(dataset ~ family_block, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = y_label,
      colour = "Representation",
      fill = "Representation",
      linetype = "Scope",
      shape = "Scope"
    ) +
    theme_m2m() +
    theme(
      axis.text.x = element_text(face = "bold", size = 6.9),
      axis.text.y = element_text(size = 7.0),
      strip.text.x = element_text(size = 8.2, face = "bold"),
      strip.text.y = element_text(angle = 0, size = 8.2),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = grid::unit(6, "mm"),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, panel_id)
}

load_task1_internal_matched_support <- function() {
  data <- load_task1_merged_pairwise()
  data <- data[
    data$row_kind == "raw" & grepl("^Internal ", as.character(data$comparison_context)),
    ,
    drop = FALSE
  ]
  data$dataset <- factor(
    ifelse(grepl("LINCS", as.character(data$comparison_context)), "LINCS", "scPerturb"),
    levels = c("LINCS", "scPerturb")
  )
  data$family_block <- factor(
    ifelse(data$analysis_family == "group_concordance", "Group", "Retrieval"),
    levels = c("Group", "Retrieval")
  )
  data$metric_display <- factor(metric_display_name(data$metric_name), levels = task1_metric_label_levels)
  data$plot_value <- metric_plot_value(data$metric_name, data$metric_value)
  data$x_pos <- as.numeric(data$metric_display) + ifelse(data$representation == "Gene", -0.18, 0.18)
  data
}

build_task1_internal_matched_support_plot <- function(data, panel_id) {
  gene <- data[
    data$representation == "Gene",
    c("dataset", "family_block", "metric_display", "metric_name", "unit_id", "x_pos", "plot_value"),
    drop = FALSE
  ]
  pathway <- data[
    data$representation == "Pathway",
    c("dataset", "family_block", "metric_display", "metric_name", "unit_id", "x_pos", "plot_value"),
    drop = FALSE
  ]
  names(gene)[names(gene) == "x_pos"] <- "gene_x"
  names(gene)[names(gene) == "plot_value"] <- "gene_value"
  names(pathway)[names(pathway) == "x_pos"] <- "pathway_x"
  names(pathway)[names(pathway) == "plot_value"] <- "pathway_value"
  paired <- merge(
    gene,
    pathway,
    by = c("dataset", "family_block", "metric_display", "metric_name", "unit_id"),
    all = FALSE
  )

  plot <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.32, colour = "#DCE3EB") +
    geom_segment(
      data = paired,
      aes(
        x = gene_x,
        xend = pathway_x,
        y = gene_value,
        yend = pathway_value,
        group = unit_id
      ),
      inherit.aes = FALSE,
      linewidth = 0.28,
      colour = "#D5DEE8",
      lineend = "round"
    ) +
    geom_point(
      data = data,
      aes(x = x_pos, y = plot_value, fill = representation),
      shape = 21,
      size = 1.9,
      stroke = 0.32,
      colour = "#1F2D3A",
      alpha = 0.82
    ) +
    scale_color_representation(drop = FALSE) +
    scale_fill_representation(drop = FALSE) +
    facet_grid(dataset ~ family_block, scales = "free_x", space = "free_x") +
    scale_x_continuous(
      breaks = seq_along(task1_metric_label_levels),
      labels = task1_metric_label_levels,
      expand = expansion(mult = c(0.04, 0.05))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = "Raw matched-triplet value",
      colour = "Representation",
      fill = "Representation"
    ) +
    theme_m2m() +
    theme(
      axis.text.x = element_text(face = "bold", size = 6.9),
      axis.text.y = element_text(size = 7.0),
      strip.text.x = element_text(size = 8.2, face = "bold"),
      strip.text.y = element_text(angle = 0, size = 8.2),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = grid::unit(6, "mm"),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, panel_id)
}

load_task1_enrichment_pairs <- function(path, entity_column) {
  data <- read_plot_ready_csv(
    path,
    required_cols = c(
      "dataset",
      entity_column,
      "row_id",
      "global_row_order",
      "row_order",
      "entity_display_label",
      "row_block",
      "col_block",
      "perturbation_type_label",
      "selection_tail",
      "selection_note",
      "gene_enrichment_score",
      "pathway_enrichment_score",
      "gene_significance_label",
      "pathway_significance_label",
      "support_value",
      "support_n",
      "support_label",
      "pair_mean_enrichment",
      "shared_across_modalities_bool",
      "shared_flag"
    )
  )
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data$shared_across_modalities_bool <- tolower(as.character(data$shared_across_modalities_bool)) == "true"
  if ("shared_flag" %in% colnames(data)) {
    data$shared_flag <- tolower(as.character(data$shared_flag)) == "true"
  } else {
    data$shared_flag <- data$shared_across_modalities_bool
  }
  data$row_block <- factor(data$row_block, levels = c("Chemical", "Genetic"))
  data$col_block <- factor(data$col_block, levels = c("LINCS", "scPerturb"))
  row_keys <- paste(data$col_block, data$row_block, data$row_id, sep = "::")
  ordered_keys <- unique(row_keys[order(data$col_block, data$row_block, data$global_row_order, data$row_id)])
  data$row_factor <- factor(
    row_keys,
    levels = rev(ordered_keys)
  )
  data$facet_value <- data$col_block
  data$entity_label <- data$entity_display_label
  data$left_annotation_primary <- ifelse(data$shared_flag, "Shared", "")
  data$left_annotation_secondary <- ifelse(
    is.na(data$selection_tail) | !nzchar(data$selection_tail),
    "",
    tools::toTitleCase(as.character(data$selection_tail))
  )
  data
}

build_bar_scoreboard <- function(data, panel_id) {
  data <- data[!is.na(data$metric_panel) & !is.na(data$context_label) & !is.na(data$representation_label) & !is.na(data$dataset_label), , drop = FALSE]
  rep_levels <- rev(task1_representation_levels[task1_representation_levels %in% unique(task1_representation_levels)])
  layout <- expand.grid(
    dataset_label = levels(data$dataset_label),
    metric_panel = levels(data$metric_panel),
    context_label = levels(data$context_label),
    representation_label = rep_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  layout$dataset_label <- factor(layout$dataset_label, levels = levels(data$dataset_label))
  layout$metric_panel <- factor(layout$metric_panel, levels = levels(data$metric_panel))
  layout$context_label <- factor(layout$context_label, levels = levels(data$context_label))
  layout$representation_label <- factor(layout$representation_label, levels = rep_levels)

  board <- merge(
    layout,
    data,
    by = c("dataset_label", "metric_panel", "context_label", "representation_label"),
    all.x = TRUE,
    sort = FALSE
  )
  board$representation_group <- factor(
    task1_representation_group(board$representation_label),
    levels = c("Gene", "Pathway", "FM")
  )
  board$column_id <- as.numeric(board$context_label)
  board$row_id <- as.numeric(board$representation_label)
  board$xmin <- board$column_id - 0.40
  board$xmax <- board$column_id + 0.40
  board$ymin <- board$row_id - 0.31
  board$ymax <- board$row_id + 0.31
  board$cell_status <- ifelse(
    !is.na(board$value),
    "data",
    ifelse(
      as.character(board$dataset_label) == "LINCS" & !(as.character(board$representation_group) %in% c("Gene", "Pathway")),
      "not_applicable",
      "empty"
    )
  )
  board$fill_display <- suppressWarnings(as.numeric(board$fill_value))
  panel_max <- ave(
    suppressWarnings(as.numeric(board$value)),
    interaction(board$dataset_label, board$metric_panel, drop = TRUE),
    FUN = function(values) {
      values <- values[is.finite(values)]
      if (!length(values)) {
        return(1)
      }
      max(values, na.rm = TRUE)
    }
  )
  panel_max[!is.finite(panel_max) | panel_max <= 0] <- 1
  board$bar_fraction <- ifelse(
    board$cell_status == "data",
    pmax(0, pmin(1, suppressWarnings(as.numeric(board$value)) / panel_max)),
    0
  )
  board$bar_xmin <- board$xmin + 0.05
  board$bar_xmax <- board$bar_xmin + 0.70 * board$bar_fraction
  board$bar_y <- board$ymin + 0.12
  board$na_label <- ifelse(board$cell_status == "not_applicable", "N/A", "")
  board$value_label <- ifelse(board$cell_status == "data", board$display_label, board$na_label)

  plot <- ggplot() +
    geom_rect(
      data = board,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_display),
      color = "#E3E8EE",
      linewidth = 0.22
    ) +
    geom_segment(
      data = board[board$cell_status == "data", , drop = FALSE],
      aes(x = bar_xmin, xend = bar_xmax, y = bar_y, yend = bar_y, color = representation_group),
      linewidth = 2.15,
      alpha = 0.95,
      lineend = "round"
    ) +
    geom_text(
      data = board,
      aes(x = column_id, y = row_id + 0.03, label = value_label),
      size = 2.15,
      fontface = "bold",
      color = "#22313F",
      lineheight = 0.9
    ) +
    scale_fill_gradientn(
      colours = c("#F8FAFC", "#EAF0F6", "#CEDCEB", "#9DB9D8"),
      values = c(0, 0.45, 0.75, 1),
      na.value = "#F3F5F7",
      guide = "none"
    ) +
    scale_color_manual(values = c(Gene = "#5E88C1", Pathway = "#C76898", FM = "#9C7AC7"), drop = FALSE) +
    facet_grid(dataset_label ~ metric_panel) +
    scale_x_continuous(
      breaks = seq_along(levels(data$context_label)),
      labels = levels(data$context_label),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(data$representation_label)),
      labels = levels(data$representation_label),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(x = NULL, y = NULL, color = "Method class") +
    theme_m2m() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 6.9, lineheight = 0.88),
      axis.text.y = element_text(face = "bold", size = 7.0),
      panel.grid = element_blank(),
      strip.text.x = element_text(size = 8.3, face = "bold"),
      strip.text.y = element_text(angle = 0, size = 8.3),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, panel_id)
}

build_annotation_strip_plot <- function(data) {
  annotation_data <- unique(
    data[, c("row_factor", "left_annotation_scope", "left_annotation_dataset", "left_annotation_perturbation"), drop = FALSE]
  )
  scope_rows <- data.frame(
    row_factor = annotation_data$row_factor,
    annotation_type = "Scope",
    annotation_value = annotation_data$left_annotation_scope,
    label_value = annotation_data$left_annotation_scope,
    stringsAsFactors = FALSE
  )
  dataset_rows <- data.frame(
    row_factor = annotation_data$row_factor,
    annotation_type = "Dataset",
    annotation_value = annotation_data$left_annotation_dataset,
    label_value = annotation_data$left_annotation_dataset,
    stringsAsFactors = FALSE
  )
  perturbation_rows <- data.frame(
    row_factor = annotation_data$row_factor,
    annotation_type = "Type",
    annotation_value = annotation_data$left_annotation_perturbation,
    label_value = annotation_data$left_annotation_perturbation,
    stringsAsFactors = FALSE
  )
  long <- rbind(scope_rows, dataset_rows, perturbation_rows)
  long$annotation_type <- factor(long$annotation_type, levels = c("Scope", "Dataset", "Type"))

  ggplot(long, aes(x = annotation_type, y = row_factor, fill = annotation_value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = label_value), size = 2.15, lineheight = 0.88, color = "#1F2D3A") +
    scale_fill_manual(values = annotation_strip_palette(), drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      plot.margin = margin(10, 2, 8, 10)
    )
}

build_paired_enrichment_panel <- function(path, entity_column, panel_id, facet_row = NULL, facet_col = NULL) {
  data <- load_task1_enrichment_pairs(path, entity_column)
  paired_enrichment_exemplar_plot(
    data,
    panel_id = panel_id,
    annotation_primary = "left_annotation_primary",
    annotation_secondary = "left_annotation_secondary",
    facet_row = facet_row,
    facet_col = facet_col
  )
}

build_panel_2A <- function() {
  data <- load_task1_scope()
  label_data <- aggregate(coverage_fraction ~ scope_block + slice_label, data = data, FUN = sum)
  label_data <- merge(
    label_data,
    unique(data[, c("scope_block", "slice_label", "count_label")]),
    by = c("scope_block", "slice_label"),
    all.x = TRUE
  )

  plot <- ggplot(data, aes(x = coverage_fraction, y = slice_label, fill = representation_class, alpha = scope_status)) +
    geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", color = "#B7C1CB") +
    geom_col(position = "stack", width = 0.72, color = "white", linewidth = 0.25) +
    geom_text(
      data = label_data,
      aes(x = 1.12, y = slice_label, label = count_label),
      inherit.aes = FALSE,
      hjust = 1,
      size = 2.35,
      color = "#6B7D8E",
      fontface = "plain"
    ) +
    facet_grid(scope_block ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(
      values = c(Gene = "#4F81BD", Pathway = "#C97B63", FM = "#8E7DBE", ALL = "#9FA6AD"),
      drop = FALSE
    ) +
    scale_alpha_manual(
      values = c(materialized = 1, available = 0.96, excluded_by_support_gate = 0.45, not_applicable_scope = 0.22),
      guide = "none"
    ) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(
      limits = c(0, 1.15),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Fraction of lawful representation surfaces", y = NULL, fill = "Representation") +
    theme_m2m() +
    theme(
      strip.text.y = element_text(angle = 0),
      axis.text.y = element_text(face = "bold"),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, "2A")
}

build_panel_2B <- function() {
  build_bar_scoreboard(load_task1_internal_performance(), "2B")
}

build_task1_crossover_summary <- function() {
  data <- load_task1_merged_pairwise()
  summary <- data[data$row_kind == "summary", , drop = FALSE]
  summary$metric_label <- factor(
    metric_display_name(summary$metric_name),
    levels = task1_metric_label_levels
  )
  summary$metric_panel <- factor(
    task1_metric_panel(summary$analysis_family, summary$metric_name),
    levels = task1_metric_panel_levels
  )
  summary$analysis_family_label <- factor(
    analysis_family_display_name(summary$analysis_family),
    levels = c("Group concordance", "Retrieval")
  )
  summary$comparison_context <- factor(
    summary$comparison_context,
    levels = c("Internal LINCS", "Internal scPerturb", "LINCS -> scPerturb", "scPerturb -> LINCS")
  )
  summary$context_index <- as.numeric(summary$comparison_context)
  summary$y_pos <- summary$context_index + ifelse(summary$representation == "Gene", -0.13, 0.13)

  gene <- summary[
    summary$representation == "Gene",
    c("comparison_context", "analysis_family", "analysis_family_label", "metric_name", "metric_label", "metric_panel", "plot_q50"),
    drop = FALSE
  ]
  pathway <- summary[
    summary$representation == "Pathway",
    c("comparison_context", "analysis_family", "analysis_family_label", "metric_name", "metric_label", "metric_panel", "plot_q50"),
    drop = FALSE
  ]
  names(gene)[names(gene) == "plot_q50"] <- "gene_q50"
  names(pathway)[names(pathway) == "plot_q50"] <- "pathway_q50"
  paired <- merge(
    gene,
    pathway,
    by = c("comparison_context", "analysis_family", "analysis_family_label", "metric_name", "metric_label", "metric_panel"),
    all = FALSE
  )
  paired$context_index <- as.numeric(paired$comparison_context)
  paired$gene_y <- paired$context_index - 0.13
  paired$pathway_y <- paired$context_index + 0.13

  stats <- unique(
    summary[
      ,
      c("comparison_context", "analysis_family", "analysis_family_label", "metric_name", "metric_label", "metric_panel", "bh_q", "test_status", "n_test_units", "median_delta"),
      drop = FALSE
    ]
  )
  stats$marker <- comparison_stats_label(stats$bh_q, stats$test_status, stats$n_test_units, stats$median_delta)
  stats$y_pos <- as.numeric(stats$comparison_context) + 0.36

  list(summary = summary, paired = paired, stats = stats[nzchar(stats$marker), , drop = FALSE])
}

build_task1_internal_only_summary <- function() {
  data <- load_task1_merged_pairwise()
  summary <- data[
    data$row_kind == "summary" & grepl("^Internal ", as.character(data$comparison_context)),
    ,
    drop = FALSE
  ]
  summary$dataset <- factor(
    ifelse(grepl("LINCS", as.character(summary$comparison_context)), "LINCS", "scPerturb"),
    levels = c("LINCS", "scPerturb")
  )
  summary$family_block <- factor(
    ifelse(summary$analysis_family == "group_concordance", "Group", "Retrieval"),
    levels = c("Group", "Retrieval")
  )
  summary$metric_display <- factor(metric_display_name(summary$metric_name), levels = task1_metric_label_levels)

  stats <- unique(
    summary[
      ,
      c("dataset", "family_block", "metric_name", "metric_display", "bh_q", "test_status", "n_test_units", "median_delta"),
      drop = FALSE
    ]
  )
  stats$inline_label <- comparison_stats_inline_label(
    stats$bh_q,
    test_status = stats$test_status,
    n_test_units = stats$n_test_units,
    median_delta = stats$median_delta
  )
  stats <- stats[nzchar(stats$inline_label), , drop = FALSE]

  list(summary = summary, stats = stats)
}

build_task1_crossover_plot <- function(plot_data, panel_id, x_label) {
  plot_data$stats$inline_label <- comparison_stats_inline_label(
    plot_data$stats$bh_q,
    test_status = plot_data$stats$test_status,
    n_test_units = plot_data$stats$n_test_units,
    median_delta = plot_data$stats$median_delta
  )
  plot_data$stats <- plot_data$stats[nzchar(plot_data$stats$inline_label), , drop = FALSE]
  plot <- ggplot() +
    geom_segment(
      data = plot_data$paired,
      aes(x = gene_q50, xend = pathway_q50, y = gene_y, yend = pathway_y),
      inherit.aes = FALSE,
      linewidth = 0.70,
      colour = "#D5DEE8",
      lineend = "round"
    ) +
    geom_segment(
      data = plot_data$summary,
      aes(x = plot_q10, xend = plot_q90, y = y_pos, yend = y_pos, colour = representation),
      linewidth = 0.30,
      alpha = 0.46
    ) +
    geom_segment(
      data = plot_data$summary,
      aes(x = plot_q25, xend = plot_q75, y = y_pos, yend = y_pos, colour = representation),
      linewidth = 1.35,
      alpha = 0.94,
      lineend = "round"
    ) +
    geom_point(
      data = plot_data$summary,
      aes(x = plot_q50, y = y_pos, fill = representation),
      shape = 21,
      size = 2.45,
      stroke = 0.38,
      colour = "#1F2D3A"
    ) +
    geom_text(
      data = plot_data$stats,
      aes(x = Inf, y = y_pos, label = inline_label),
      inherit.aes = FALSE,
      hjust = 1.02,
      size = 1.56,
      colour = "#6B1830",
      lineheight = 0.92
    ) +
    scale_color_representation(drop = FALSE) +
    scale_fill_representation(drop = FALSE) +
    facet_wrap(~metric_panel, scales = "free_x", nrow = 1) +
    scale_y_continuous(
      breaks = seq_along(levels(plot_data$summary$comparison_context)),
      labels = levels(plot_data$summary$comparison_context),
      expand = expansion(mult = c(0.08, 0.18))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = x_label,
      y = NULL,
      colour = "Representation",
      fill = "Representation"
    ) +
    theme_m2m() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(face = "bold", size = 6.8),
      axis.text.x = element_text(size = 6.8, face = "bold"),
      strip.text.x = element_text(size = 7.7, face = "bold"),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, panel_id)
}

build_panel_2C <- function() {
  build_task1_internal_to_cross_degradation_plot(
    load_task1_degradation(),
    panel_id = "2C",
    y_label = "Internal-to-cross degradation summary"
  )
}

build_panel_2D <- function() {
  plot_data <- build_task1_internal_only_summary()
  data <- plot_data$summary
  data$plot_q10 <- metric_plot_value(data$metric_name, data$q10_value)
  data$plot_q25 <- metric_plot_value(data$metric_name, data$q25_value)
  data$plot_q50 <- metric_plot_value(data$metric_name, data$q50_value)
  data$plot_q75 <- metric_plot_value(data$metric_name, data$q75_value)
  data$plot_q90 <- metric_plot_value(data$metric_name, data$q90_value)
  stats <- plot_data$stats
  panel_max <- stats::aggregate(plot_q90 ~ dataset + family_block + metric_display, data = data, FUN = max)
  panel_min <- stats::aggregate(plot_q10 ~ dataset + family_block + metric_display, data = data, FUN = min)
  stats <- merge(stats, panel_max, by = c("dataset", "family_block", "metric_display"), all.x = TRUE)
  stats <- merge(stats, panel_min, by = c("dataset", "family_block", "metric_display"), all.x = TRUE)
  stats$y_pad <- pmax(0.05, 0.10 * (stats$plot_q90 - stats$plot_q10))
  stats$y_pos <- stats$plot_q90 + stats$y_pad

  dodge <- position_dodge(width = 0.52)

  plot <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.32, colour = "#DCE3EB") +
    geom_linerange(
      data = data,
      aes(
        x = metric_display,
        ymin = plot_q10,
        ymax = plot_q90,
        colour = representation,
        group = representation
      ),
      position = dodge,
      linewidth = 0.30,
      alpha = 0.34
    ) +
    geom_linerange(
      data = data,
      aes(
        x = metric_display,
        ymin = plot_q25,
        ymax = plot_q75,
        colour = representation,
        group = representation
      ),
      position = dodge,
      linewidth = 0.72,
      alpha = 0.74
    ) +
    geom_point(
      data = data,
      aes(
        x = metric_display,
        y = plot_q50,
        fill = representation,
        group = representation
      ),
      position = dodge,
      shape = 21,
      size = 2.45,
      stroke = 0.36,
      colour = "#1F2D3A"
    ) +
    geom_text(
      data = stats,
      aes(
        x = metric_display,
        y = y_pos,
        label = inline_label
      ),
      inherit.aes = FALSE,
      size = 1.56,
      lineheight = 0.9,
      colour = "#6B1830",
      show.legend = FALSE
    ) +
    scale_color_representation(drop = FALSE) +
    scale_fill_representation(drop = FALSE) +
    facet_grid(dataset ~ family_block, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL,
      y = "Internal Gene-vs-Pathway summary",
      colour = "Representation",
      fill = "Representation"
    ) +
    theme_m2m() +
    theme(
      axis.text.x = element_text(face = "bold", size = 6.9),
      axis.text.y = element_text(size = 7.0),
      strip.text.x = element_text(size = 8.2, face = "bold"),
      strip.text.y = element_text(angle = 0, size = 8.2),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = grid::unit(6, "mm"),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, "2D")
}

build_panel_2E <- function() {
  build_paired_enrichment_panel(
    "figure2/figure2_panel_2e_cell_line_high_concordance_summary.csv",
    "cell_line",
    "2E",
    facet_row = "row_block",
    facet_col = "col_block"
  )
}

build_panel_2F <- function() {
  build_paired_enrichment_panel(
    "figure2/figure2_panel_2f_target_high_concordance_summary.csv",
    "target_token",
    "2F",
    facet_row = "row_block",
    facet_col = "col_block"
  )
}

compose_figure2 <- function() {
  panels <- list(
    build_panel_2A(),
    build_panel_2B(),
    build_panel_2C(),
    build_panel_2D(),
    build_panel_2E(),
    build_panel_2F()
  )
  compose_manuscript_figure("figure2", panels, figure2_panel_ids, ncol = 2)
}

figure2_panel_builders <- function() {
  list(
    `2A` = build_panel_2A,
    `2B` = build_panel_2B,
    `2C` = build_panel_2C,
    `2D` = build_panel_2D,
    `2E` = build_panel_2E,
    `2F` = build_panel_2F
  )
}

figure2_panel_dimensions <- list(
  `2A` = list(width_mm = 200, height_mm = 120),
  `2B` = list(width_mm = 200, height_mm = 120),
  `2C` = list(width_mm = 240, height_mm = 126),
  `2D` = list(width_mm = 218, height_mm = 156),
  `2E` = list(width_mm = 210, height_mm = 160),
  `2F` = list(width_mm = 210, height_mm = 160)
)

render_figure2 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_figure2()
  render_composed_figure(
    "figure2",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export,
    panel_ids = figure2_panel_ids,
    remaining_placeholder_flags = character(),
    all_real_panels_included = TRUE
  )
}

render_figure2_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE,
  separate_legends = FALSE
) {
  render_panel_set(
    figure_id = "figure2",
    panel_ids = figure2_panel_ids,
    panel_builders = figure2_panel_builders(),
    panel_dimensions = figure2_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export,
    separate_legends = separate_legends
  )
}

if (sys.nframe() == 0) {
  render_figure2()
}
