suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source(file.path("plotting", "R", "figures", "figure2.R"))

extended_figure2_panel_ids <- c("EF2C")

build_panel_EF2C <- function() {
  build_task1_internal_matched_support_plot(
    load_task1_internal_matched_support(),
    panel_id = "EF2C"
  )
}

compose_extended_figure2 <- function() {
  panels <- list(build_panel_EF2C())
  compose_manuscript_figure("extended_figure2", panels, c("EF2C"), ncol = 1)
}

extended_figure2_panel_builders <- function() {
  list(
    `EF2C` = build_panel_EF2C
  )
}

extended_figure2_panel_dimensions <- list(
  `EF2C` = list(width_mm = 230, height_mm = 122)
)

render_extended_figure2 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_extended_figure2()
  render_composed_figure(
    "extended_figure2",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

render_extended_figure2_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  render_panel_set(
    figure_id = "extended_figure2",
    panel_ids = extended_figure2_panel_ids,
    panel_builders = extended_figure2_panel_builders(),
    panel_dimensions = extended_figure2_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

if (sys.nframe() == 0) {
  render_extended_figure2()
}
