panel_summary_lines <- function(source_reports) {
  unlist(
    lapply(
      source_reports,
      function(report) {
        if (!report$exists) {
          return(sprintf("%s missing; generate before final render", basename(report$path)))
        }
        sprintf(
          "%s | %s rows | %s cols",
          basename(report$path),
          format(nrow(report$data), big.mark = ","),
          format(ncol(report$data), big.mark = ",")
        )
      }
    ),
    use.names = FALSE
  )
}

panel_contract_card <- function(panel_id, plot_type, source_specs, body_lines = character()) {
  labels <- panel_labels(panel_id)
  reports <- lapply(source_specs, inspect_panel_source)
  summary_lines <- panel_summary_lines(reports)
  source_lines <- vapply(reports, function(report) source_label(report), character(1))
  caption_lines <- c(
    sprintf("Plot type: %s", plot_type),
    if (length(source_lines) > 0) sprintf("Sources: %s", paste(source_lines, collapse = ", ")) else NULL,
    summary_lines,
    body_lines,
    if (nzchar(labels$scope_constraint_text)) labels$scope_constraint_text else NULL,
    if (nzchar(labels$panel_footer)) labels$panel_footer else NULL
  )
  card_data <- data.frame(x = 0, y = 0)
  ggplot2::ggplot(card_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_blank() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      hjust = 0,
      vjust = 1,
      size = 3.4,
      family = "Arial",
      label = paste(c(labels$panel_title, labels$panel_subtitle, caption_lines), collapse = "\n")
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::labs(title = panel_id)
}

panel_text_card <- function(panel_id, plot_type, body_lines = character()) {
  labels <- panel_labels(panel_id)
  card_data <- data.frame(x = 0, y = 0)
  ggplot2::ggplot(card_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_blank() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      hjust = 0,
      vjust = 1,
      size = 3.4,
      family = "Arial",
      label = paste(
        c(
          labels$panel_title,
          labels$panel_subtitle,
          sprintf("Plot type: %s", plot_type),
          body_lines,
          if (nzchar(labels$scope_constraint_text)) labels$scope_constraint_text else NULL,
          if (nzchar(labels$panel_footer)) labels$panel_footer else NULL
        ),
        collapse = "\n"
      )
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::labs(title = panel_id)
}

compose_manuscript_figure <- function(figure_id, panels, panel_ids, ncol = NULL) {
  invisible(figure_id)
  invisible(panel_ids)
  patchwork::wrap_plots(panels, ncol = ncol, guides = "collect") &
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 12, 10, 12),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.justification = "center"
    )
}
