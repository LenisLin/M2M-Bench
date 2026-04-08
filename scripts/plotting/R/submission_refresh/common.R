suppressPackageStartupMessages({
  submission_refresh_repo_root <- function() {
    env_root <- Sys.getenv("M2M_REPO_ROOT", unset = "")
    if (nzchar(env_root)) {
      return(env_root)
    }
    if (file.exists(file.path(getwd(), "plotting", "R", "submission_refresh", "common.R"))) {
      return(getwd())
    }
    stop(
      "Could not resolve the M2M repo root. Run from the repo root or set M2M_REPO_ROOT.",
      call. = FALSE
    )
  }

  if (!identical(Sys.getenv("CONDA_DEFAULT_ENV", unset = ""), "Spatial")) {
    stop(
      "Submission-refresh R code must run inside the conda environment 'Spatial'.",
      call. = FALSE
    )
  }
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

source(file.path(submission_refresh_repo_root(), "plotting", "R", "shared", "theme.R"))

submission_refresh_data_root <- function() {
  Sys.getenv(
    "M2M_SUBMISSION_REFRESH_DATA_ROOT",
    unset = "/mnt/NAS_21T/ProjectData/M2M/R_Vis_Ready"
  )
}

submission_refresh_output_root <- function() {
  Sys.getenv(
    "M2M_SUBMISSION_REFRESH_OUTPUT_ROOT",
    unset = "/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/submission_refresh"
  )
}

submission_refresh_palette <- function() {
  list(
    standard = "#0072B2",
    systema = "#D55E00",
    gene = "#0072B2",
    pathway = "#009E73",
    robust_high = "#009E73",
    protocol_sensitive = "#CC79A7",
    robust_low = "#D55E00",
    neutral_fill = "#EEF2F6",
    neutral_line = "#B8C4CF",
    dark_text = "#1F2D3A",
    grid = "#E7ECF1",
    effect_low = "#0072B2",
    effect_mid = "#F7F8FA",
    effect_high = "#D55E00"
  )
}

submission_view_palette <- function() {
  pal <- submission_refresh_palette()
  c(Standard = pal$standard, Systema = pal$systema)
}

submission_track_palette <- function() {
  pal <- submission_refresh_palette()
  c(Gene = pal$gene, Pathway = pal$pathway)
}

submission_class_palette <- function() {
  pal <- submission_refresh_palette()
  c(
    Robust_High = pal$robust_high,
    Protocol_Sensitive = pal$protocol_sensitive,
    Robust_Low = pal$robust_low
  )
}

submission_signed_scale <- function(name) {
  pal <- submission_refresh_palette()
  scale_fill_gradient2(
    low = pal$effect_low,
    mid = pal$effect_mid,
    high = pal$effect_high,
    midpoint = 0,
    name = name
  )
}

submission_signed_colour_scale <- function(name) {
  pal <- submission_refresh_palette()
  scale_colour_gradient2(
    low = pal$effect_low,
    mid = pal$effect_mid,
    high = pal$effect_high,
    midpoint = 0,
    name = name
  )
}

theme_submission_refresh <- function(base_size = 9.5) {
  theme_m2m(base_size = base_size) +
    theme(
      panel.grid.major.x = element_line(color = submission_refresh_palette()$grid, linewidth = 0.2),
      panel.grid.major.y = element_line(color = submission_refresh_palette()$grid, linewidth = 0.2),
      plot.tag = element_text(
        face = "bold",
        size = base_size + 2,
        color = submission_refresh_palette()$dark_text
      ),
      plot.tag.position = c(0, 1),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
}

tag_submission_panels <- function(plot) {
  plot + plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 13,
        color = submission_refresh_palette()$dark_text
      ),
      plot.tag.position = c(0, 1)
    )
  )
}

submission_refresh_csv <- function(relative_path, required_cols = NULL) {
  path <- file.path(submission_refresh_data_root(), relative_path)
  if (!file.exists(path)) {
    stop(sprintf("Submission refresh source file is missing: %s", path), call. = FALSE)
  }
  data <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      stop(
        sprintf(
          "File %s is missing required columns: %s",
          basename(path),
          paste(missing_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  data
}

submission_refresh_dir <- function() {
  out_dir <- submission_refresh_output_root()
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

submission_refresh_pdf_path <- function(filename) {
  file.path(submission_refresh_dir(), filename)
}

save_submission_refresh_pdf <- function(plot, filename, width_mm, height_mm) {
  output_path <- submission_refresh_pdf_path(filename)
  device_fn <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(
    filename = output_path,
    plot = plot,
    device = device_fn,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    units = "in",
    bg = "white"
  )
  output_path
}

metric_label_submission <- function(values) {
  map <- c(
    Cosine = "Pairwise cosine",
    DEG_PCC = "DEG PCC",
    JaccAbs = "Jaccard abs",
    Mean_MRR = "Mean MRR",
    Mean_Success = "Mean success",
    NegLogEDist = "E-distance score",
    CentroidCosine = "Centroid cosine",
    ED_score = "E-distance score",
    Success_CRISPR2Drug = "CRISPR -> Drug success",
    Mean_CentroidCosine = "Mean centroid cosine",
    Mean_NegEDist = "Mean E-distance score",
    Median_Rank = "Median rank"
  )
  unknown_values <- setdiff(unique(values), names(map))
  if (length(unknown_values) > 0) {
    stop(
      sprintf(
        "Unknown metric identifiers in submission refresh renderer: %s",
        paste(sort(unknown_values), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  out <- unname(map[values])
  out
}

protocol_label_submission <- function(values) {
  recode <- c(Dose_Corr = "Dose", Time_Corr = "Time")
  out <- unname(recode[values])
  out[is.na(out)] <- values[is.na(out)]
  out
}

clean_label <- function(values) {
  gsub("_", " ", values, fixed = TRUE)
}
