suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source(file.path("plotting", "R", "figures", "figure3.R"))

extended_figure5_panel_ids <- c("EF5A", "EF5B", "EF5C")

build_panel_EF5A <- function() {
  build_scperturb_tier_detail_scoreboard("The Cleanest Hits", "EF5A")
}

build_panel_EF5B <- function() {
  build_scperturb_tier_detail_scoreboard("The Family Hits", "EF5B")
}

build_panel_EF5C <- function() {
  build_scperturb_tier_detail_scoreboard("The Promiscuous Hits", "EF5C")
}

compose_extended_figure5 <- function() {
  panels <- list(build_panel_EF5A(), build_panel_EF5B(), build_panel_EF5C())
  compose_manuscript_figure("extended_figure5", panels, c("EF5A", "EF5B", "EF5C"), ncol = 2)
}

extended_figure5_panel_builders <- function() {
  list(
    `EF5A` = build_panel_EF5A,
    `EF5B` = build_panel_EF5B,
    `EF5C` = build_panel_EF5C
  )
}

extended_figure5_panel_dimensions <- list(
  `EF5A` = list(width_mm = 206, height_mm = 140),
  `EF5B` = list(width_mm = 206, height_mm = 140),
  `EF5C` = list(width_mm = 206, height_mm = 140)
)

render_extended_figure5 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_extended_figure5()
  render_composed_figure(
    "extended_figure5",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

render_extended_figure5_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  render_panel_set(
    figure_id = "extended_figure5",
    panel_ids = extended_figure5_panel_ids,
    panel_builders = extended_figure5_panel_builders(),
    panel_dimensions = extended_figure5_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

if (sys.nframe() == 0) {
  render_extended_figure5()
}
