suppressPackageStartupMessages({
  source(file.path("plotting", "R", "figures", "figure1.R"))
  source(file.path("plotting", "R", "figures", "figure2.R"))
  source(file.path("plotting", "R", "figures", "figure3.R"))
  source(file.path("plotting", "R", "figures", "extended_figure1.R"))
  source(file.path("plotting", "R", "figures", "extended_figure2.R"))
  source(file.path("plotting", "R", "figures", "extended_figure4.R"))
  source(file.path("plotting", "R", "figures", "extended_figure5.R"))
  source(file.path("plotting", "R", "figures", "extended_figure6.R"))
  source(file.path("plotting", "R", "figures", "extended_figure7.R"))
  source(file.path("plotting", "R", "figures", "extended_figure8.R"))
  source(file.path("plotting", "R", "figures", "extended_figure9.R"))
})

panel_render_specs <- function() {
  list(
    list(figure_id = "figure1", panel_ids = figure1_panel_ids, render_fn = render_figure1_panels),
    list(figure_id = "figure2", panel_ids = figure2_panel_ids, render_fn = render_figure2_panels),
    list(figure_id = "figure3", panel_ids = figure3_panel_ids, render_fn = render_figure3_panels),
    list(figure_id = "extended_figure2", panel_ids = extended_figure2_panel_ids, render_fn = render_extended_figure2_panels),
    list(figure_id = "extended_figure4", panel_ids = extended_figure4_panel_ids, render_fn = render_extended_figure4_panels),
    list(figure_id = "extended_figure5", panel_ids = extended_figure5_panel_ids, render_fn = render_extended_figure5_panels),
    list(figure_id = "extended_figure6", panel_ids = extended_figure6_panel_ids, render_fn = render_extended_figure6_panels),
    list(figure_id = "extended_figure7", panel_ids = extended_figure7_panel_ids, render_fn = render_extended_figure7_panels),
    list(figure_id = "extended_figure8", panel_ids = extended_figure8_panel_ids, render_fn = render_extended_figure8_panels)
  )
}

validate_rendered_panel_set <- function(spec, results) {
  if (!is.list(results) || length(results) == 0) {
    stop(
      sprintf("%s returned no standalone panel metadata.", spec$figure_id),
      call. = FALSE
    )
  }
  if (length(results) != length(spec$panel_ids)) {
    stop(
      sprintf(
        "%s returned %d panels, expected %d (%s).",
        spec$figure_id,
        length(results),
        length(spec$panel_ids),
        paste(spec$panel_ids, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  observed_figure_ids <- vapply(results, function(meta) meta$figure_id, character(1))
  if (!all(observed_figure_ids == spec$figure_id)) {
    stop(
      sprintf(
        "%s returned unexpected figure ids: %s.",
        spec$figure_id,
        paste(unique(observed_figure_ids), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  observed_panel_ids <- vapply(results, function(meta) meta$panel_id, character(1))
  if (anyDuplicated(observed_panel_ids)) {
    stop(
      sprintf(
        "%s returned duplicate panel ids: %s.",
        spec$figure_id,
        paste(unique(observed_panel_ids[duplicated(observed_panel_ids)]), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  missing_panel_ids <- setdiff(spec$panel_ids, observed_panel_ids)
  extra_panel_ids <- setdiff(observed_panel_ids, spec$panel_ids)
  if (length(missing_panel_ids) > 0 || length(extra_panel_ids) > 0) {
    detail <- c(
      if (length(missing_panel_ids) > 0) {
        sprintf("missing: %s", paste(missing_panel_ids, collapse = ", "))
      },
      if (length(extra_panel_ids) > 0) {
        sprintf("unexpected: %s", paste(extra_panel_ids, collapse = ", "))
      }
    )
    stop(
      sprintf("%s panel inventory mismatch (%s).", spec$figure_id, paste(detail, collapse = "; ")),
      call. = FALSE
    )
  }
  missing_outputs <- vapply(
    results,
    function(meta) !file.exists(meta$output_path),
    logical(1)
  )
  if (any(missing_outputs)) {
    stop(
      sprintf(
        "%s produced missing panel PDFs for %s.",
        spec$figure_id,
        paste(observed_panel_ids[missing_outputs], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  missing_legends <- vapply(
    results,
    function(meta) {
      legend_path <- meta$legend_output_path
      if (is.null(legend_path) || is.na(legend_path) || !nzchar(legend_path)) {
        return(FALSE)
      }
      !file.exists(legend_path)
    },
    logical(1)
  )
  if (any(missing_legends)) {
    stop(
      sprintf(
        "%s produced missing legend PDFs for %s.",
        spec$figure_id,
        paste(observed_panel_ids[missing_legends], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  ordered <- results[match(spec$panel_ids, observed_panel_ids)]
  if (any(vapply(ordered, is.null, logical(1)))) {
    stop(sprintf("%s panel ordering could not be resolved.", spec$figure_id), call. = FALSE)
  }
  ordered
}

render_panel_spec <- function(spec) {
  args <- list(export = TRUE)
  results <- tryCatch(
    do.call(spec$render_fn, args),
    error = function(err) {
      stop(
        sprintf("%s panel render failed: %s", spec$figure_id, conditionMessage(err)),
        call. = FALSE
      )
    }
  )
  validate_rendered_panel_set(spec, results)
}

write_panel_manifest <- function(results) {
  manifest_path <- file.path(figures_root(), "panel_export_manifest.tsv")
  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  if (!length(results)) {
    stop("Cannot write a panel export manifest with zero panel rows.", call. = FALSE)
  }
  analysis_root_value <- resolved_analysis_root()
  plot_ready_root_value <- resolved_plot_ready_root()
  figures_root_value <- resolved_figures_root()
  rows <- lapply(
    results,
    function(meta) {
      data.frame(
        figure_id = meta$figure_id,
        panel_id = meta$panel_id,
        output_path = meta$output_path,
        legend_output_path = ifelse(is.null(meta$legend_output_path), NA_character_, meta$legend_output_path),
        base_width_mm = ifelse(is.null(meta$base_width_mm), NA_real_, meta$base_width_mm),
        base_height_mm = ifelse(is.null(meta$base_height_mm), NA_real_, meta$base_height_mm),
        width_mm = meta$width_mm,
        height_mm = meta$height_mm,
        export_unit = ifelse(is.null(meta$export_unit), "cm", meta$export_unit),
        qc_policy = ifelse(is.null(meta$qc_policy), NA_character_, meta$qc_policy),
        qc_status = ifelse(is.null(meta$qc_status), NA_character_, meta$qc_status),
        overlap_count = ifelse(is.null(meta$overlap_count), NA_real_, meta$overlap_count),
        qc_detail = ifelse(is.null(meta$qc_detail), NA_character_, meta$qc_detail),
        rerender_count = ifelse(is.null(meta$rerender_count), 0, meta$rerender_count),
        legend_width_mm = ifelse(is.null(meta$legend_width_mm), NA_real_, meta$legend_width_mm),
        legend_height_mm = ifelse(is.null(meta$legend_height_mm), NA_real_, meta$legend_height_mm),
        device_used = meta$device_used,
        font_family = meta$font_family,
        analysis_root = analysis_root_value,
        plot_ready_root = plot_ready_root_value,
        figures_root = figures_root_value,
        stringsAsFactors = FALSE
      )
    }
  )
  manifest <- do.call(rbind, rows)
  utils::write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)
  manifest_path
}

run_spatial_revision_panels <- function(verbose = FALSE) {
  specs <- panel_render_specs()
  results <- list()
  for (spec in specs) {
    if (isTRUE(verbose)) {
      cat(paste0("render_start\t", spec$figure_id, "\n"))
    }
    results <- c(results, render_panel_spec(spec))
    if (isTRUE(verbose)) {
      cat(paste0("render_done\t", spec$figure_id, "\n"))
    }
  }
  manifest_path <- write_panel_manifest(results)
  list(
    manifest_path = manifest_path,
    results = results
  )
}

main <- function() {
  outcome <- run_spatial_revision_panels(verbose = TRUE)
  cat(paste0("panel_export_manifest\t", outcome$manifest_path, "\n"))
  for (meta in outcome$results) {
    cat(
      paste(
        meta$figure_id,
        meta$panel_id,
        meta$output_path,
        ifelse(is.null(meta$legend_output_path), "", meta$legend_output_path),
        meta$width_mm,
        meta$height_mm,
        ifelse(is.null(meta$qc_status), "", meta$qc_status),
        sep = "\t"
      ),
      "\n",
      sep = ""
    )
  }
  invisible(outcome)
}

if (sys.nframe() == 0) {
  main()
}
