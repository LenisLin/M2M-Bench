suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source(file.path("plotting", "R", "figures", "figure3.R"))

extended_figure7_panel_ids <- c("EF7A", "EF7B", "EF7C", "EF7D")

modifier_bool <- function(values) {
  tolower(as.character(values)) %in% c("true", "t", "1")
}

modifier_panel_subtitle <- function(representation_value, modifier_type_value) {
  paste(
    representation_value,
    "|",
    ifelse(modifier_type_value == "time", "Time effect", "Dose effect")
  )
}

load_task2_modifier_spearman_stats <- function() {
  data <- read_analysis_csv(
    "figure3_task2_c2g_dose_time_spearman_stats.csv",
    required_cols = c(
      "dataset",
      "cell_line",
      "representation",
      "target",
      "modifier_type",
      "metric",
      "spearman_rho",
      "p_value",
      "bh_q",
      "n_obs",
      "testable_bool",
      "untestable_reason"
    )
  )
  data$metric_name <- data$metric
  data <- data[
    data$representation %in% c("Gene", "Pathway") &
      data$modifier_type %in% c("time", "dose"),
    ,
    drop = FALSE
  ]
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data$representation <- factor(data$representation, levels = c("Gene", "Pathway"))
  data$modifier_type <- factor(data$modifier_type, levels = c("time", "dose"))
  data$spearman_rho <- suppressWarnings(as.numeric(data$spearman_rho))
  data$bh_q <- suppressWarnings(as.numeric(data$bh_q))
  data$n_obs <- suppressWarnings(as.numeric(data$n_obs))
  data$testable_bool <- modifier_bool(data$testable_bool)
  data$point_y <- ifelse(
    !is.na(data$bh_q),
    -log10(pmax(data$bh_q, 1e-12)),
    NA_real_
  )
  data$significance_class <- ifelse(
    data$testable_bool & !is.na(data$bh_q) & data$bh_q <= 0.05 & data$spearman_rho < 0,
    "Negative q<=0.05",
    ifelse(
      data$testable_bool & !is.na(data$bh_q) & data$bh_q <= 0.05 & data$spearman_rho > 0,
      "Positive q<=0.05",
      "NS"
    )
  )
  data$significance_class <- factor(
    data$significance_class,
    levels = c("Negative q<=0.05", "NS", "Positive q<=0.05")
  )
  data
}

build_modifier_dataset_annotations <- function(all_data, tested_data) {
  datasets <- unique(as.character(all_data$dataset))
  rows <- lapply(
    datasets,
    function(dataset_name) {
      all_subset <- all_data[as.character(all_data$dataset) == dataset_name, , drop = FALSE]
      tested_subset <- tested_data[as.character(tested_data$dataset) == dataset_name, , drop = FALSE]
      untestable_subset <- all_subset[!all_subset$testable_bool, , drop = FALSE]
      reason <- ""
      if (nrow(untestable_subset) > 0) {
        reason_counts <- sort(table(untestable_subset$untestable_reason), decreasing = TRUE)
        reason <- names(reason_counts)[[1]]
      }
      label <- paste0(
        "tested=", format(nrow(tested_subset), big.mark = ","),
        "\nuntestable=", format(nrow(untestable_subset), big.mark = ",")
      )
      if (nzchar(reason)) {
        label <- paste0(label, "\nmain skip=", truncate_label(reason, width = 26))
      }
      data.frame(
        dataset = factor(dataset_name, levels = levels(all_data$dataset)),
        x = Inf,
        y = Inf,
        label = label,
        stringsAsFactors = FALSE
      )
    }
  )
  do.call(rbind, rows)
}

build_modifier_highlights <- function(tested_data) {
  highlight <- tested_data[tested_data$bh_q <= 0.05, , drop = FALSE]
  if (nrow(highlight) == 0) {
    return(highlight)
  }
  highlight <- highlight[
    order(as.character(highlight$dataset), highlight$bh_q, -abs(highlight$spearman_rho), highlight$target),
    ,
    drop = FALSE
  ]
  highlight$label_rank <- ave(
    seq_len(nrow(highlight)),
    highlight$dataset,
    FUN = seq_along
  )
  highlight <- highlight[highlight$label_rank <= 3, , drop = FALSE]
  highlight$target_label <- truncate_label(highlight$target, width = 16)
  highlight
}

build_spearman_volcano_panel <- function(representation_value, modifier_type_value, panel_id) {
  data <- load_task2_modifier_spearman_stats()
  data <- data[
    as.character(data$representation) == representation_value &
      as.character(data$modifier_type) == modifier_type_value,
    ,
    drop = FALSE
  ]
  if (nrow(data) == 0) {
    return(
      panel_text_card(
        panel_id,
        "spearman volcano missing",
        body_lines = c(
          sprintf("No %s rows were found for %s.", modifier_type_value, representation_value),
          "Generate figure3_task2_c2g_dose_time_spearman_stats.csv before rendering EF7."
        )
      )
    )
  }

  tested <- data[
    data$testable_bool &
      is.finite(data$spearman_rho) &
      is.finite(data$point_y),
    ,
    drop = FALSE
  ]
  if (nrow(tested) == 0) {
    return(
      panel_text_card(
        panel_id,
        "spearman volcano empty",
        body_lines = c(
          sprintf("All %s rows for %s are currently untestable.", modifier_type_value, representation_value),
          "Review minimum-observation and modifier-variation gates."
        )
      )
    )
  }

  annotations <- build_modifier_dataset_annotations(data, tested)
  highlights <- build_modifier_highlights(tested)

  plot <- ggplot(tested, aes(x = spearman_rho, y = point_y)) +
    geom_hline(yintercept = -log10(0.05), linewidth = 0.28, linetype = "dashed", colour = "#D5DEE8") +
    geom_vline(xintercept = 0, linewidth = 0.28, linetype = "dashed", colour = "#D5DEE8") +
    geom_point(aes(colour = significance_class, size = n_obs), alpha = 0.86) +
    geom_text(
      data = highlights,
      aes(x = spearman_rho, y = point_y, label = target_label),
      inherit.aes = FALSE,
      check_overlap = TRUE,
      nudge_y = 0.12,
      size = 2.05,
      colour = "#44576B"
    ) +
    geom_text(
      data = annotations,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.04,
      vjust = 1.08,
      size = 2.1,
      lineheight = 0.95,
      colour = "#526273"
    ) +
    facet_wrap(~dataset, nrow = 1) +
    scale_color_manual(
      values = c(
        `Negative q<=0.05` = "#C76898",
        NS = "#9AA8B6",
        `Positive q<=0.05` = "#5E88C1"
      ),
      drop = FALSE
    ) +
    scale_size_continuous(range = c(1.2, 2.8), breaks = c(3, 5, 10)) +
    coord_cartesian(xlim = c(-1, 1), clip = "off") +
    labs(
      x = "Spearman rho",
      y = "-log10(BH q)",
      colour = "Direction / FDR",
      size = "n obs",
      subtitle = modifier_panel_subtitle(representation_value, modifier_type_value)
    ) +
    theme_m2m() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6.8, face = "bold"),
      axis.text.y = element_text(size = 6.7),
      strip.text = element_text(size = 8.0, face = "bold"),
      plot.subtitle = element_text(size = 8.4, face = "bold", colour = "#526273"),
      legend.position = "bottom"
    )

  finish_panel_plot(plot, panel_id)
}

build_panel_EF7A <- function() {
  build_spearman_volcano_panel("Gene", "time", "EF7A")
}

build_panel_EF7B <- function() {
  build_spearman_volcano_panel("Gene", "dose", "EF7B")
}

build_panel_EF7C <- function() {
  build_spearman_volcano_panel("Pathway", "time", "EF7C")
}

build_panel_EF7D <- function() {
  build_spearman_volcano_panel("Pathway", "dose", "EF7D")
}

compose_extended_figure7 <- function() {
  panels <- list(build_panel_EF7A(), build_panel_EF7B(), build_panel_EF7C(), build_panel_EF7D())
  compose_manuscript_figure("extended_figure7", panels, c("EF7A", "EF7B", "EF7C", "EF7D"), ncol = 2)
}

extended_figure7_panel_builders <- function() {
  list(
    `EF7A` = build_panel_EF7A,
    `EF7B` = build_panel_EF7B,
    `EF7C` = build_panel_EF7C,
    `EF7D` = build_panel_EF7D
  )
}

extended_figure7_panel_dimensions <- list(
  `EF7A` = list(width_mm = 206, height_mm = 150),
  `EF7B` = list(width_mm = 206, height_mm = 150),
  `EF7C` = list(width_mm = 206, height_mm = 150),
  `EF7D` = list(width_mm = 206, height_mm = 150)
)

render_extended_figure7 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_extended_figure7()
  render_composed_figure(
    "extended_figure7",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

render_extended_figure7_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  render_panel_set(
    figure_id = "extended_figure7",
    panel_ids = extended_figure7_panel_ids,
    panel_builders = extended_figure7_panel_builders(),
    panel_dimensions = extended_figure7_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

if (sys.nframe() == 0) {
  render_extended_figure7()
}
