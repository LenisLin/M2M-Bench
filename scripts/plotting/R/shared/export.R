figure_dimensions_config_path <- function() {
  file.path(plotting_config_root(), "figure_dimensions.csv")
}

figure_dimensions_table <- NULL

get_figure_dimensions_table <- function() {
  if (is.null(figure_dimensions_table)) {
    figure_dimensions_table <<- read_m2m_csv(
      figure_dimensions_config_path(),
      required_cols = c(
        "figure_id",
        "default_width_mm",
        "width_hard_limit_mm",
        "default_height_mm",
        "height_policy",
        "orientation",
        "split_before_render",
        "export_filename_pattern"
      )
    )
  }
  figure_dimensions_table
}

default_figure_output_path <- function(figure_id, output_dir = NULL) {
  table <- get_figure_dimensions_table()
  row <- table[table$figure_id == figure_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop(sprintf("Expected exactly one dimension row for %s.", figure_id), call. = FALSE)
  }
  relative_pattern <- row$export_filename_pattern[[1]]
  if (is.null(output_dir)) {
    return(file.path(figures_root(), relative_pattern))
  }
  file.path(output_dir, basename(relative_pattern))
}

default_panel_output_dir <- function(figure_id, output_dir = NULL) {
  if (!is.null(output_dir)) {
    return(output_dir)
  }
  file.path(figures_root(), figure_id, "panels")
}

default_panel_output_path <- function(figure_id, panel_id, output_dir = NULL) {
  file.path(
    default_panel_output_dir(figure_id, output_dir = output_dir),
    sprintf("%s_panel_%s.pdf", figure_id, tolower(panel_id))
  )
}

default_panel_legend_output_path <- function(figure_id, panel_id, output_dir = NULL) {
  file.path(
    default_panel_output_dir(figure_id, output_dir = output_dir),
    sprintf("%s_panel_%s_legend.pdf", figure_id, tolower(panel_id))
  )
}

resolve_figure_dimensions <- function(figure_id, width_mm = NULL, height_mm = NULL) {
  table <- get_figure_dimensions_table()
  row <- table[table$figure_id == figure_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop(sprintf("Expected exactly one dimension row for %s.", figure_id), call. = FALSE)
  }
  resolved_width <- if (is.null(width_mm)) row$default_width_mm[[1]] else width_mm
  resolved_height <- if (is.null(height_mm)) row$default_height_mm[[1]] else height_mm
  width_hard_limit_mm <- row$width_hard_limit_mm[[1]]
  if (!is.numeric(resolved_width) || resolved_width <= 0) {
    stop("width_mm must be a positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(resolved_height) || resolved_height <= 0) {
    stop("height_mm must be a positive numeric value.", call. = FALSE)
  }
  if (resolved_width > width_hard_limit_mm || resolved_width > 180) {
    stop(
      sprintf(
        "Requested width %.1f mm exceeds the frozen hard limit of %.1f mm.",
        resolved_width,
        min(width_hard_limit_mm, 180)
      ),
      call. = FALSE
    )
  }
  list(
    figure_id = figure_id,
    width_mm = resolved_width,
    height_mm = resolved_height,
    width_hard_limit_mm = width_hard_limit_mm,
    height_policy = row$height_policy[[1]],
    orientation = row$orientation[[1]],
    split_before_render = row$split_before_render[[1]],
    export_filename_pattern = row$export_filename_pattern[[1]]
  )
}

mm_to_inches <- function(value_mm) {
  value_mm / 25.4
}

cm_to_mm <- function(value_cm) {
  value_cm * 10
}

base_panel_dimensions_mm <- function() {
  list(
    width_mm = cm_to_mm(20),
    height_mm = cm_to_mm(12),
    unit = "cm"
  )
}

panel_qc_policy <- function(panel_id) {
  if (grepl("^1B", panel_id)) {
    return("both")
  }
  if (panel_id %in% c("2E", "2F", "3B", "3D", "EF4A", "EF4B", "EF6A", "EF6B", "EF8C")) {
    return("height_first")
  }
  if (panel_id %in% c("2B", "2C", "3E", "3F", "EF8A", "EF8B")) {
    return("width_first")
  }
  "standard"
}

next_qc_dimensions <- function(width_mm, height_mm, policy) {
  max_width_mm <- 260
  if (policy == "both") {
    return(list(width_mm = min(width_mm + 20, max_width_mm), height_mm = height_mm + 20))
  }
  if (policy == "height_first") {
    return(list(width_mm = width_mm, height_mm = height_mm + 20))
  }
  if (policy == "width_first") {
    next_width_mm <- min(width_mm + 20, max_width_mm)
    next_height_mm <- if (identical(next_width_mm, width_mm)) height_mm + 20 else height_mm
    return(list(width_mm = next_width_mm, height_mm = next_height_mm))
  }
  list(width_mm = min(width_mm + 10, max_width_mm), height_mm = height_mm + 10)
}

find_headless_chrome <- function() {
  candidates <- c(
    Sys.getenv("M2M_CHROME_BIN", unset = ""),
    Sys.which("google-chrome"),
    Sys.which("google-chrome-stable"),
    Sys.which("chromium"),
    Sys.which("chromium-browser")
  )
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates) == 0) {
    stop(
      "No Chrome/Chromium binary found for htmlwidget PDF export.",
      call. = FALSE
    )
  }
  candidates[[1]]
}

export_htmlwidget_pdf <- function(
  widget,
  output_path,
  width_mm,
  height_mm
) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    stop("The 'htmlwidgets' package is required for widget PDF export.", call. = FALSE)
  }
  chrome_bin <- find_headless_chrome()
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  html_path <- tempfile("m2m_panel_", fileext = ".html")
  chrome_profile_dir <- tempfile("m2m_chrome_profile_")
  on.exit(unlink(html_path), add = TRUE)
  on.exit(unlink(chrome_profile_dir, recursive = TRUE), add = TRUE)
  dir.create(chrome_profile_dir, recursive = TRUE, showWarnings = FALSE)
  htmlwidgets::saveWidget(widget, file = html_path, selfcontained = TRUE)
  width_px <- max(800, round(mm_to_inches(width_mm) * 96))
  height_px <- max(600, round(mm_to_inches(height_mm) * 96))
  args <- c(
    "--headless",
    "--no-sandbox",
    "--disable-gpu",
    "--disable-dev-shm-usage",
    "--disable-breakpad",
    "--disable-crash-reporter",
    "--run-all-compositor-stages-before-draw",
    "--virtual-time-budget=10000",
    sprintf("--user-data-dir=%s", normalizePath(chrome_profile_dir, winslash = "/", mustWork = TRUE)),
    sprintf("--window-size=%d,%d", width_px, height_px),
    "--no-pdf-header-footer",
    paste0("--print-to-pdf=", normalizePath(output_path, winslash = "/", mustWork = FALSE)),
    paste0("file://", normalizePath(html_path, winslash = "/", mustWork = TRUE))
  )
  status <- suppressWarnings(system2(chrome_bin, args = args, stdout = TRUE, stderr = TRUE))
  if (!file.exists(output_path)) {
    stop(
      sprintf(
        "Chrome PDF export failed for %s. Chrome output: %s",
        output_path,
        paste(status, collapse = " | ")
      ),
      call. = FALSE
    )
  }
  invisible(
    list(
      output_path = output_path,
      device_used = "google-chrome --headless",
      font_family = NA_character_,
      font_path = NA_character_,
      export_log = paste(status, collapse = "\n")
    )
  )
}

parse_pdf_word_boxes <- function(output_path) {
  pdftotext_bin <- Sys.which("pdftotext")
  if (!nzchar(pdftotext_bin)) {
    stop(
      "The 'pdftotext' binary is required for standalone panel QC.",
      call. = FALSE
    )
  }
  html <- suppressWarnings(
    system2(
      pdftotext_bin,
      args = c("-bbox-layout", output_path, "-"),
      stdout = TRUE,
      stderr = TRUE
    )
  )
  word_lines <- grep("<word ", html, value = TRUE)
  if (length(word_lines) == 0) {
    return(data.frame())
  }
  extract_attr <- function(line, key) {
    sub(sprintf(".*%s=\"([^\"]+)\".*", key), "\\1", line)
  }
  data.frame(
    xmin = suppressWarnings(as.numeric(vapply(word_lines, extract_attr, character(1), key = "xMin"))),
    ymin = suppressWarnings(as.numeric(vapply(word_lines, extract_attr, character(1), key = "yMin"))),
    xmax = suppressWarnings(as.numeric(vapply(word_lines, extract_attr, character(1), key = "xMax"))),
    ymax = suppressWarnings(as.numeric(vapply(word_lines, extract_attr, character(1), key = "yMax"))),
    text = gsub("^.*<word[^>]*>|</word>.*$", "", word_lines),
    stringsAsFactors = FALSE
  )
}

detect_pdf_text_overlap <- function(output_path) {
  if (!file.exists(output_path)) {
    return(list(status = "missing_pdf", overlap_count = NA_integer_, detail = "Output PDF missing."))
  }
  boxes <- parse_pdf_word_boxes(output_path)
  if (nrow(boxes) < 2) {
    return(list(status = "pass", overlap_count = 0L, detail = "No overlapping text boxes detected."))
  }
  overlap_count <- 0L
  for (idx in seq_len(nrow(boxes) - 1)) {
    dx <- pmin(boxes$xmax[[idx]], boxes$xmax[(idx + 1):nrow(boxes)]) -
      pmax(boxes$xmin[[idx]], boxes$xmin[(idx + 1):nrow(boxes)])
    dy <- pmin(boxes$ymax[[idx]], boxes$ymax[(idx + 1):nrow(boxes)]) -
      pmax(boxes$ymin[[idx]], boxes$ymin[(idx + 1):nrow(boxes)])
    overlap_count <- overlap_count + sum(dx > 0.75 & dy > 0.75, na.rm = TRUE)
    if (overlap_count > 0) {
      break
    }
  }
  if (overlap_count > 0) {
    return(list(status = "overlap_detected", overlap_count = overlap_count, detail = "PDF word boxes overlap."))
  }
  list(status = "pass", overlap_count = 0L, detail = "No overlapping text boxes detected.")
}

pdf_font_available <- function(font_family) {
  if (!requireNamespace("systemfonts", quietly = TRUE)) {
    return(FALSE)
  }
  font_match <- systemfonts::match_font(font_family)
  path <- font_match$path
  is.character(path) && length(path) == 1 && nzchar(path) && file.exists(path)
}

resolve_pdf_font <- function(font_family) {
  if (!requireNamespace("systemfonts", quietly = TRUE)) {
    stop(
      "The 'systemfonts' package is required for final PDF export.",
      call. = FALSE
    )
  }
  font_match <- systemfonts::match_font(font_family)
  path <- font_match$path
  if (!is.character(path) || length(path) != 1 || !nzchar(path) || !file.exists(path)) {
    stop(
      sprintf("Font family '%s' could not be resolved to an installed system font.", font_family),
      call. = FALSE
    )
  }
  list(
    family = font_family,
    path = normalizePath(path, winslash = "/", mustWork = TRUE)
  )
}

assert_pdf_font_family <- function(font_family) {
  if (isTRUE(capabilities("cairo"))) {
    return(invisible(resolve_pdf_font(font_family)))
  }
  stop(
    paste(
      "Cairo PDF support is unavailable in this R session.",
      "Final manuscript export requires Cairo-backed vector PDF output."
    ),
    call. = FALSE
  )
}

export_editable_pdf <- function(
  plot,
  output_path,
  width_mm = 178,
  height_mm,
  device = "pdf",
  transparent_bg = FALSE,
  font_family = "Arial"
) {
  if (!identical(device, "pdf")) {
    stop("The frozen manuscript export device is editable vector PDF only.", call. = FALSE)
  }
  if (missing(height_mm) || is.null(height_mm)) {
    stop("height_mm must be provided to export_editable_pdf().", call. = FALSE)
  }
  if (!is.character(output_path) || length(output_path) != 1 || !nzchar(output_path)) {
    stop("output_path must be a single non-empty PDF path.", call. = FALSE)
  }
  if (!grepl("\\.pdf$", output_path, ignore.case = TRUE)) {
    stop("output_path must end with .pdf for final manuscript export.", call. = FALSE)
  }
  if (width_mm > 260) {
    stop("width_mm must be <= 260 mm for final manuscript export.", call. = FALSE)
  }
  if (isTRUE(transparent_bg)) {
    warning(
      "Transparent PDF backgrounds are allowed only when explicitly required; default opaque white is safer for Adobe editing.",
      call. = FALSE
    )
  }
  resolved_font <- assert_pdf_font_family(font_family)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(
    filename = output_path,
    width = mm_to_inches(width_mm),
    height = mm_to_inches(height_mm),
    onefile = FALSE,
    family = resolved_font$family,
    bg = if (isTRUE(transparent_bg)) "transparent" else "white",
    fallback_resolution = 300
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  if (grid::is.grob(plot)) {
    grid::grid.newpage()
    grid::grid.draw(plot)
  } else {
    print(plot)
  }
  invisible(
    list(
      output_path = output_path,
      device_used = "grDevices::cairo_pdf",
      font_family = resolved_font$family,
      font_path = resolved_font$path
    )
  )
}

extract_plot_legend <- function(plot) {
  grob <- ggplot2::ggplotGrob(plot)
  legend <- gtable::gtable_filter(grob, "guide-box", fixed = TRUE)
  if (length(legend$grobs) == 0) {
    return(NULL)
  }
  keep <- which(!vapply(legend$grobs, inherits, logical(1), "zeroGrob"))
  if (length(keep) == 0) {
    return(NULL)
  }
  legend$grobs <- legend$grobs[keep]
  legend$layout <- legend$layout[keep, , drop = FALSE]
  legend
}

resolve_panel_legend_dimensions <- function(
  panel_id,
  panel_width_mm,
  legend_dimensions = NULL
) {
  if (!is.null(legend_dimensions)) {
    dims <- legend_dimensions[[panel_id]]
    if (!is.null(dims$width_mm) && !is.null(dims$height_mm)) {
      return(dims)
    }
  }
  list(
    width_mm = min(160, panel_width_mm),
    height_mm = 30
  )
}

render_composed_figure <- function(
  figure_id,
  plot,
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE,
  panel_ids = character(),
  remaining_placeholder_flags = character(),
  all_real_panels_included = NULL
) {
  dims <- resolve_figure_dimensions(figure_id, width_mm = width_mm, height_mm = height_mm)
  output_path <- default_figure_output_path(figure_id, output_dir = output_dir)
  export_result <- list(
    output_path = output_path,
    device_used = "grDevices::cairo_pdf",
    font_family = font_family,
    font_path = NULL
  )
  if (isTRUE(export)) {
    export_result <- export_editable_pdf(
      plot = plot,
      output_path = output_path,
      width_mm = dims$width_mm,
      height_mm = dims$height_mm,
      device = "pdf",
      transparent_bg = transparent_bg,
      font_family = font_family
    )
  }
  placeholder_flags <- unique(remaining_placeholder_flags[!is.na(remaining_placeholder_flags) & nzchar(remaining_placeholder_flags)])
  if (is.null(all_real_panels_included)) {
    all_real_panels_included <- length(placeholder_flags) == 0
  }
  invisible(
    list(
      figure_id = figure_id,
      output_path = export_result$output_path,
      width_mm = dims$width_mm,
      height_mm = dims$height_mm,
      device_used = export_result$device_used,
      font_family = export_result$font_family,
      font_path = export_result$font_path,
      panel_ids = panel_ids,
      all_real_panels_included = isTRUE(all_real_panels_included),
      remaining_placeholder_flags = placeholder_flags,
      exported = isTRUE(export),
      plot = plot
    )
  )
}

render_individual_panel <- function(
  figure_id,
  panel_id,
  plot,
  output_dir = NULL,
  width_mm = 178,
  height_mm,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  if (missing(height_mm) || is.null(height_mm)) {
    stop("height_mm must be provided to render_individual_panel().", call. = FALSE)
  }
  if (!is.numeric(width_mm) || width_mm <= 0 || width_mm > 260) {
    stop("render_individual_panel() requires width_mm in the range (0, 260].", call. = FALSE)
  }
  if (!is.numeric(height_mm) || height_mm <= 0) {
    stop("render_individual_panel() requires positive height_mm.", call. = FALSE)
  }
  output_path <- default_panel_output_path(figure_id, panel_id, output_dir = output_dir)
  qc_policy <- panel_qc_policy(panel_id)
  base_width_mm <- width_mm
  base_height_mm <- height_mm
  rerender_count <- 0L
  export_result <- list(
    output_path = output_path,
    device_used = "grDevices::cairo_pdf",
    font_family = font_family,
    font_path = NULL
  )
  if (isTRUE(export)) {
    render_once <- function(width_value, height_value) {
      if (inherits(plot, "htmlwidget")) {
        return(
          export_htmlwidget_pdf(
            widget = plot,
            output_path = output_path,
            width_mm = width_value,
            height_mm = height_value
          )
        )
      }
      export_editable_pdf(
        plot = plot,
        output_path = output_path,
        width_mm = width_value,
        height_mm = height_value,
        device = "pdf",
        transparent_bg = transparent_bg,
        font_family = font_family
      )
    }
    export_result <- render_once(width_mm, height_mm)
    qc_result <- detect_pdf_text_overlap(output_path)
    while (identical(qc_result$status, "overlap_detected") && rerender_count < 3L) {
      next_dims <- next_qc_dimensions(width_mm, height_mm, qc_policy)
      width_mm <- next_dims$width_mm
      height_mm <- next_dims$height_mm
      rerender_count <- rerender_count + 1L
      export_result <- render_once(width_mm, height_mm)
      qc_result <- detect_pdf_text_overlap(output_path)
    }
  } else {
    qc_result <- list(status = "not_exported", overlap_count = NA_integer_, detail = "Panel not exported.")
  }
  invisible(
    list(
      figure_id = figure_id,
      panel_id = panel_id,
      output_path = export_result$output_path,
      base_width_mm = base_width_mm,
      base_height_mm = base_height_mm,
      width_mm = width_mm,
      height_mm = height_mm,
      export_unit = "cm",
      device_used = export_result$device_used,
      font_family = export_result$font_family,
      font_path = export_result$font_path,
      qc_policy = qc_policy,
      qc_status = qc_result$status,
      overlap_count = qc_result$overlap_count,
      qc_detail = qc_result$detail,
      rerender_count = rerender_count,
      exported = isTRUE(export),
      plot = plot
    )
  )
}

render_individual_panel_legend <- function(
  figure_id,
  panel_id,
  legend_plot,
  output_dir = NULL,
  width_mm = 140,
  height_mm = 30,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  output_path <- default_panel_legend_output_path(figure_id, panel_id, output_dir = output_dir)
  export_result <- list(
    output_path = output_path,
    device_used = "grDevices::cairo_pdf",
    font_family = font_family,
    font_path = NULL
  )
  if (isTRUE(export)) {
    export_result <- export_editable_pdf(
      plot = legend_plot,
      output_path = output_path,
      width_mm = width_mm,
      height_mm = height_mm,
      device = "pdf",
      transparent_bg = transparent_bg,
      font_family = font_family
    )
  }
  invisible(
    list(
      figure_id = figure_id,
      panel_id = panel_id,
      legend_output_path = export_result$output_path,
      legend_width_mm = width_mm,
      legend_height_mm = height_mm,
      legend_device_used = export_result$device_used,
      legend_font_family = export_result$font_family,
      legend_font_path = export_result$font_path,
      legend_exported = isTRUE(export)
    )
  )
}

render_panel_set <- function(
  figure_id,
  panel_ids,
  panel_builders,
  panel_dimensions,
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE,
  separate_legends = FALSE,
  legend_builders = NULL,
  legend_dimensions = NULL,
  legend_output_dir = NULL
) {
  target_dir <- default_panel_output_dir(figure_id, output_dir = output_dir)
  lapply(
    panel_ids,
    function(panel_id) {
      builder <- panel_builders[[panel_id]]
      dims <- panel_dimensions[[panel_id]]
      base_dims <- base_panel_dimensions_mm()
      if (!is.function(builder)) {
        stop(sprintf("No panel builder registered for %s %s.", figure_id, panel_id), call. = FALSE)
      }
      if (is.null(dims$width_mm) || is.na(dims$width_mm)) {
        dims$width_mm <- base_dims$width_mm
      }
      if (is.null(dims$height_mm) || is.na(dims$height_mm)) {
        dims$height_mm <- base_dims$height_mm
      }
      plot <- builder()
      panel_plot <- if (isTRUE(separate_legends) && !grid::is.grob(plot) && !inherits(plot, "htmlwidget")) {
        plot + ggplot2::theme(legend.position = "none")
      } else {
        plot
      }
      panel_meta <- render_individual_panel(
        figure_id = figure_id,
        panel_id = panel_id,
        plot = panel_plot,
        output_dir = target_dir,
        width_mm = dims$width_mm,
        height_mm = dims$height_mm,
        transparent_bg = transparent_bg,
        font_family = font_family,
        export = export
      )
      panel_meta$legend_output_path <- NA_character_
      panel_meta$legend_width_mm <- NA_real_
      panel_meta$legend_height_mm <- NA_real_
      panel_meta$legend_exported <- FALSE
      if (isTRUE(separate_legends)) {
        legend_plot <- NULL
        if (!is.null(legend_builders) && is.function(legend_builders[[panel_id]])) {
          legend_plot <- legend_builders[[panel_id]]()
        } else if (!grid::is.grob(plot) && !inherits(plot, "htmlwidget")) {
          legend_plot <- extract_plot_legend(plot)
        }
        if (!is.null(legend_plot)) {
          legend_dims <- resolve_panel_legend_dimensions(
            panel_id = panel_id,
            panel_width_mm = dims$width_mm,
            legend_dimensions = legend_dimensions
          )
          legend_meta <- render_individual_panel_legend(
            figure_id = figure_id,
            panel_id = panel_id,
            legend_plot = legend_plot,
            output_dir = default_panel_output_dir(figure_id, output_dir = legend_output_dir),
            width_mm = legend_dims$width_mm,
            height_mm = legend_dims$height_mm,
            transparent_bg = transparent_bg,
            font_family = font_family,
            export = export
          )
          panel_meta$legend_output_path <- legend_meta$legend_output_path
          panel_meta$legend_width_mm <- legend_meta$legend_width_mm
          panel_meta$legend_height_mm <- legend_meta$legend_height_mm
          panel_meta$legend_exported <- legend_meta$legend_exported
        }
      }
      panel_meta
    }
  )
}
