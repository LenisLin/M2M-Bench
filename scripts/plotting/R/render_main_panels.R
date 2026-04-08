suppressPackageStartupMessages({
  source(file.path("plotting", "R", "figures", "figure1.R"))
  source(file.path("plotting", "R", "figures", "figure2.R"))
  source(file.path("plotting", "R", "figures", "figure3.R"))
})

report_panel_exports <- function(results) {
  for (meta in results) {
    cat(
      paste(
        meta[["figure_id"]],
        meta[["panel_id"]],
        meta[["output_path"]],
        meta[["width_mm"]],
        meta[["height_mm"]],
        meta[["device_used"]],
        meta[["font_family"]],
        sep = "\t"
      ),
      "\n",
      sep = ""
    )
  }
}

report_panel_exports(render_figure1_panels(export = TRUE))
report_panel_exports(render_figure2_panels(export = TRUE))
report_panel_exports(render_figure3_panels(export = TRUE))
