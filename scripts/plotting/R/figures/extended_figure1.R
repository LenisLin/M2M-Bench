suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source(file.path("plotting", "R", "figures", "figure1.R"))

extended_figure1_panel_ids <- c("EF1A", "EF1B", "EF1C", "EF1D", "EF1E")

load_extended_figure1_registry <- function() {
  data <- read_plot_ready_csv("extended/extended_figure1_support_registry.csv")
  data$item_label_short <- truncate_label(data$item_label, width = 34)
  data$item_label_short <- factor(data$item_label_short, levels = rev(unique(data$item_label_short)))
  data
}

build_panel_EF1A <- function() {
  palette <- get_m2m_palette()
  plot <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0.72, ymax = 1, fill = palette$section_band[["task_design_fill"]], color = NA) +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 0.72, fill = "#FAFBFD", color = NA) +
    annotate("text", x = 0.04, y = 0.91, hjust = 0, vjust = 1, label = "Registry framing", size = 4.2, fontface = "bold", color = "#1F2D3A") +
    annotate(
      "text",
      x = 0.04,
      y = 0.64,
      hjust = 0,
      vjust = 1,
      label = paste(
        "Figure 1 remains doc-derived.",
        "The registry panels below show the governing support files, scope rules, and result objects that feed the downstream manuscript figures.",
        sep = "\n"
      ),
      size = 3.5,
      color = "#334155",
      lineheight = 1.05
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
  finish_panel_plot(plot, "EF1A", legend_position = "none")
}

build_registry_panel <- function(panel_id, registry_type, facet_field = "panel_id") {
  invisible(facet_field)
  data <- load_extended_figure1_registry()
  data <- data[data$registry_type == registry_type, , drop = FALSE]
  plot <- ggplot(
    data,
    aes(x = status, y = item_label_short, fill = status)
  ) +
    geom_tile(color = "#FFFFFF", linewidth = 0.28) +
    geom_text(aes(label = basename(source_file)), size = 2.15, color = "#111111", lineheight = 0.88) +
    scale_fill_status(drop = FALSE) +
    ggplot2::theme(
      axis.text.x = element_text(size = 7.1, face = "bold"),
      axis.text.y = element_text(size = 7.1),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  finish_panel_plot(plot, panel_id)
}

build_panel_EF1B <- function() {
  build_registry_panel("EF1B", "coverage")
}

build_panel_EF1C <- function() {
  build_registry_panel("EF1C", "lawful_scope")
}

build_panel_EF1D <- function() {
  build_registry_panel("EF1D", "representation_modifier")
}

build_panel_EF1E <- function() {
  build_registry_panel("EF1E", "result_object")
}

compose_extended_figure1 <- function() {
  panels <- list(
    build_panel_EF1A(),
    build_panel_EF1B(),
    build_panel_EF1C(),
    build_panel_EF1D(),
    build_panel_EF1E()
  )
  compose_manuscript_figure(
    "extended_figure1",
    panels,
    c("EF1A", "EF1B", "EF1C", "EF1D", "EF1E"),
    ncol = 2
  )
}

extended_figure1_panel_builders <- function() {
  list(
    `EF1A` = build_panel_EF1A,
    `EF1B` = build_panel_EF1B,
    `EF1C` = build_panel_EF1C,
    `EF1D` = build_panel_EF1D,
    `EF1E` = build_panel_EF1E
  )
}

extended_figure1_panel_dimensions <- list(
  `EF1A` = list(width_mm = 178, height_mm = 74),
  `EF1B` = list(width_mm = 178, height_mm = 104),
  `EF1C` = list(width_mm = 178, height_mm = 104),
  `EF1D` = list(width_mm = 178, height_mm = 104),
  `EF1E` = list(width_mm = 178, height_mm = 104)
)

render_extended_figure1 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_extended_figure1()
  render_composed_figure(
    "extended_figure1",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

render_extended_figure1_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  render_panel_set(
    figure_id = "extended_figure1",
    panel_ids = extended_figure1_panel_ids,
    panel_builders = extended_figure1_panel_builders(),
    panel_dimensions = extended_figure1_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

if (sys.nframe() == 0) {
  render_extended_figure1()
}
