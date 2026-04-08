representation_palette_config_path <- function() {
  file.path(plotting_config_root(), "representation_palette.csv")
}

palette_values_config_path <- function() {
  file.path(plotting_config_root(), "palette_values.csv")
}

validate_hex_column <- function(values) {
  grepl("^#[0-9A-Fa-f]{6}$", values)
}

named_palette_values <- function(config, group_name) {
  rows <- config[config$group == group_name, , drop = FALSE]
  if (nrow(rows) == 0) {
    stop(sprintf("palette_values.csv is missing group '%s'.", group_name), call. = FALSE)
  }
  stats::setNames(rows$hex, rows$name)
}

load_palette_values <- function(config_path = palette_values_config_path()) {
  config <- read_m2m_csv(
    config_path,
    required_cols = c("group", "name", "hex")
  )
  if (!all(validate_hex_column(config$hex))) {
    stop(
      "palette_values.csv contains invalid hex values.",
      call. = FALSE
    )
  }
  list(
    section_band = named_palette_values(config, "section_band"),
    auxiliary_count = named_palette_values(config, "auxiliary_count"),
    direction = named_palette_values(config, "direction"),
    perturbation = named_palette_values(config, "perturbation"),
    status = named_palette_values(config, "status"),
    heatmap = named_palette_values(config, "heatmap")
  )
}

load_m2m_palette <- function(
  representation_config_path = representation_palette_config_path(),
  palette_values_path = palette_values_config_path()
) {
  representation_config <- read_m2m_csv(
    representation_config_path,
    required_cols = c(
      "representation",
      "hex",
      "shape",
      "linetype",
      "point_fill",
      "point_colour",
      "fill_behavior",
      "status"
    )
  )
  if (!all(validate_hex_column(representation_config$hex))) {
    stop(
      "representation_palette.csv contains placeholder or invalid representation hex values.",
      call. = FALSE
    )
  }
  if (!all(validate_hex_column(representation_config$point_fill))) {
    stop(
      "representation_palette.csv contains invalid point_fill hex values.",
      call. = FALSE
    )
  }
  if (!all(validate_hex_column(representation_config$point_colour))) {
    stop(
      "representation_palette.csv contains invalid point_colour hex values.",
      call. = FALSE
    )
  }
  aux <- load_palette_values(palette_values_path)
  list(
    section_band = aux$section_band,
    auxiliary_count = aux$auxiliary_count,
    representation = stats::setNames(representation_config$hex, representation_config$representation),
    representation_shape = stats::setNames(representation_config$shape, representation_config$representation),
    representation_linetype = stats::setNames(representation_config$linetype, representation_config$representation),
    representation_point_fill = stats::setNames(
      representation_config$point_fill,
      representation_config$representation
    ),
    representation_point_colour = stats::setNames(
      representation_config$point_colour,
      representation_config$representation
    ),
    representation_fill_behavior = stats::setNames(
      representation_config$fill_behavior,
      representation_config$representation
    ),
    direction = aux$direction,
    perturbation = aux$perturbation,
    status = aux$status,
    status_fill = aux$status,
    heatmap = aux$heatmap
  )
}

m2m_palette <- NULL

get_m2m_palette <- function() {
  if (is.null(m2m_palette)) {
    m2m_palette <<- load_m2m_palette()
  }
  m2m_palette
}

scale_color_representation <- function(...) {
  ggplot2::scale_color_manual(values = get_m2m_palette()$representation, ...)
}

scale_fill_representation <- function(...) {
  ggplot2::scale_fill_manual(values = get_m2m_palette()$representation, ...)
}

scale_shape_representation <- function(...) {
  ggplot2::scale_shape_manual(values = get_m2m_palette()$representation_shape, ...)
}

scale_linetype_representation <- function(...) {
  ggplot2::scale_linetype_manual(values = get_m2m_palette()$representation_linetype, ...)
}

scale_color_direction <- function(...) {
  ggplot2::scale_color_manual(values = get_m2m_palette()$direction, ...)
}

scale_fill_direction <- function(...) {
  ggplot2::scale_fill_manual(values = get_m2m_palette()$direction, ...)
}

scale_color_perturbation <- function(...) {
  ggplot2::scale_color_manual(values = get_m2m_palette()$perturbation, ...)
}

scale_fill_perturbation <- function(...) {
  ggplot2::scale_fill_manual(values = get_m2m_palette()$perturbation, ...)
}

scale_fill_status <- function(...) {
  ggplot2::scale_fill_manual(values = get_m2m_palette()$status_fill, ...)
}

compact_colorbar_guide <- function(
  width_mm = 30,
  height_mm = 3.2
) {
  ggplot2::guide_colorbar(
    title.position = "top",
    barwidth = grid::unit(width_mm, "mm"),
    barheight = grid::unit(height_mm, "mm")
  )
}

scale_fill_percentile_heatmap <- function(name = "Within-metric percentile", ...) {
  palette <- get_m2m_palette()$heatmap
  ggplot2::scale_fill_gradientn(
    colours = c("#FAFBFD", "#E5EDF6", palette[["mid"]], palette[["high"]]),
    values = c(0, 0.3, 0.6, 1),
    guide = compact_colorbar_guide(),
    name = name,
    ...
  )
}

scale_fill_reference_diverging <- function(name, midpoint = 0, ...) {
  ggplot2::scale_fill_gradient2(
    low = "#C76898",
    mid = "#F7F8FA",
    high = "#5E88C1",
    midpoint = midpoint,
    guide = compact_colorbar_guide(),
    name = name,
    ...
  )
}

scale_fill_heatmap_gradient <- function(...) {
  palette <- get_m2m_palette()$heatmap
  ggplot2::scale_fill_gradient2(
    low = palette[["low"]],
    mid = palette[["mid"]],
    high = palette[["high"]],
    midpoint = 0,
    guide = compact_colorbar_guide(),
    ...
  )
}
