analysis_root <- function() {
  Sys.getenv(
    "M2M_ANALYSIS_ROOT",
    unset = "/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis"
  )
}

default_analysis_root <- function() {
  "/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis"
}

plot_ready_root <- function() {
  Sys.getenv(
    "M2M_PLOT_READY_ROOT",
    unset = "/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/plot_ready"
  )
}

default_plot_ready_root <- function() {
  "/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/plot_ready"
}

figures_root <- function() {
  Sys.getenv(
    "M2M_FIGURES_ROOT",
    unset = "/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/figures"
  )
}

resolved_analysis_root <- function() {
  primary <- analysis_root()
  fallback <- default_analysis_root()
  candidate <- resolve_existing_input_path(primary, fallback)
  normalizePath(candidate, winslash = "/", mustWork = FALSE)
}

resolved_plot_ready_root <- function() {
  primary <- plot_ready_root()
  fallback <- default_plot_ready_root()
  candidate <- resolve_existing_input_path(primary, fallback)
  normalizePath(candidate, winslash = "/", mustWork = FALSE)
}

resolved_figures_root <- function() {
  normalizePath(figures_root(), winslash = "/", mustWork = FALSE)
}

plotting_config_root <- function() {
  file.path("plotting", "R", "config")
}

analysis_path <- function(name) {
  file.path(analysis_root(), name)
}

plot_ready_path <- function(name) {
  file.path(plot_ready_root(), name)
}

resolve_existing_input_path <- function(primary_path, fallback_path = NULL) {
  candidates <- c(primary_path)
  if (!is.null(fallback_path) && !identical(primary_path, fallback_path)) {
    candidates <- c(candidates, fallback_path)
  }
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }
  primary_path
}

resolved_analysis_path <- function(name) {
  resolve_existing_input_path(
    analysis_path(name),
    file.path(default_analysis_root(), name)
  )
}

resolved_plot_ready_path <- function(name) {
  resolve_existing_input_path(
    plot_ready_path(name),
    file.path(default_plot_ready_root(), name)
  )
}

figure_output_dir <- function(figure_name) {
  file.path(figures_root(), figure_name)
}

assert_required_cols <- function(data, cols) {
  missing <- setdiff(cols, colnames(data))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(data)
}

read_m2m_csv <- function(path, required_cols = NULL) {
  if (!file.exists(path)) {
    stop(sprintf("Required CSV does not exist: %s", path), call. = FALSE)
  }
  data <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(required_cols)) {
    assert_required_cols(data, required_cols)
  }
  data
}

read_analysis_csv <- function(name, required_cols = NULL) {
  read_m2m_csv(resolved_analysis_path(name), required_cols = required_cols)
}

read_plot_ready_csv <- function(name, required_cols = NULL) {
  read_m2m_csv(resolved_plot_ready_path(name), required_cols = required_cols)
}

resolve_panel_source <- function(source) {
  kind <- source$kind
  name <- source$name
  if (identical(kind, "analysis")) {
    path <- resolved_analysis_path(name)
    return(list(path = path, data = read_m2m_csv(path)))
  }
  if (identical(kind, "plot_ready")) {
    path <- resolved_plot_ready_path(name)
    return(list(path = path, data = read_m2m_csv(path)))
  }
  stop(sprintf("Unknown source kind: %s", kind), call. = FALSE)
}

inspect_panel_source <- function(source) {
  kind <- source$kind
  name <- source$name
  path <- if (identical(kind, "analysis")) {
    resolved_analysis_path(name)
  } else if (identical(kind, "plot_ready")) {
    resolved_plot_ready_path(name)
  } else {
    stop(sprintf("Unknown source kind: %s", kind), call. = FALSE)
  }
  exists <- file.exists(path)
  data <- NULL
  if (exists) {
    required_cols <- source$required_cols
    data <- read_m2m_csv(path, required_cols = required_cols)
  }
  list(
    kind = kind,
    name = name,
    path = path,
    exists = exists,
    data = data
  )
}

source_label <- function(source) {
  sprintf("%s/%s", source$kind, source$name)
}
