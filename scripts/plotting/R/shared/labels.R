panel_label_config_path <- function() {
  file.path("plotting", "R", "config", "panel_labels.csv")
}

panel_label_table <- NULL

get_panel_label_table <- function() {
  if (is.null(panel_label_table)) {
    panel_label_table <<- read_m2m_csv(
      panel_label_config_path(),
      required_cols = c(
        "figure_id",
        "panel_id",
        "panel_title",
        "panel_subtitle",
        "panel_footer",
        "scope_constraint_text"
      )
    )
  }
  panel_label_table
}

panel_labels <- function(panel_id) {
  table <- get_panel_label_table()
  row <- table[table$panel_id == panel_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop(sprintf("Expected exactly one label row for panel %s", panel_id), call. = FALSE)
  }
  as.list(row[1, ])
}

panel_labels_for_figure <- function(figure_id) {
  table <- get_panel_label_table()
  rows <- table[table$figure_id == figure_id, , drop = FALSE]
  if (nrow(rows) == 0) {
    stop(sprintf("Expected at least one label row for figure %s", figure_id), call. = FALSE)
  }
  rows
}

non_empty_unique <- function(values) {
  values <- values[!is.na(values) & nzchar(values)]
  unique(values)
}

figure_scope_footer <- function(panel_ids) {
  table <- get_panel_label_table()
  rows <- table[match(panel_ids, table$panel_id), , drop = FALSE]
  if (any(is.na(rows$panel_id))) {
    missing <- panel_ids[is.na(rows$panel_id)]
    stop(
      sprintf("Missing panel label rows for: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  lines <- c()
  for (index in seq_len(nrow(rows))) {
    row <- rows[index, , drop = FALSE]
    lines <- c(
      lines,
      non_empty_unique(row$scope_constraint_text),
      non_empty_unique(row$panel_footer)
    )
  }
  paste(non_empty_unique(lines), collapse = "\n")
}

figure_display_title <- function(figure_id) {
  if (grepl("^extended_figure", figure_id)) {
    return(sprintf("Extended Figure %s", sub("^extended_figure", "", figure_id)))
  }
  sprintf("Figure %s", sub("^figure", "", figure_id))
}
