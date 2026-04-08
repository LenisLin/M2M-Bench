suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source(file.path("plotting", "R", "shared", "io.R"))
source(file.path("plotting", "R", "shared", "palette.R"))
source(file.path("plotting", "R", "shared", "labels.R"))
source(file.path("plotting", "R", "shared", "theme.R"))
source(file.path("plotting", "R", "shared", "assembly.R"))
source(file.path("plotting", "R", "shared", "export.R"))

figure1_panel_ids <- c("1B_LINCS", "1B_SCPERTURB")

figure1_support_cache <- NULL

figure1_cell_line_rank_palette <- c(
  "#1E5A88",
  "#C05B3E",
  "#4C7F60",
  "#7B5EA7",
  "#AD7A32",
  "#2C8C99",
  "#B14763",
  "#7A6F42",
  "#4E73A8",
  "#466A57",
  "#99627A",
  "#A86E46",
  "#5E8A74",
  "#8B6A9F",
  "#6C82B5",
  "#B65C55",
  "#5F748B",
  "#8B8B45",
  "#3F8E8B",
  "#A87038"
)

figure1_others_fill <- function() {
  "#E7EDF3"
}

figure1_outer_others_fill <- function() {
  "#F2F5F8"
}

figure1_inner_ring_fill <- function() {
  palette <- get_m2m_palette()$perturbation
  c(
    Chemical = palette[["chemical"]],
    Genetic = palette[["genetic"]]
  )
}

figure1_task_ring_fill <- function() {
  c(
    Task1 = "#6C97C4",
    Task2 = "#D6866B"
  )
}

mix_with_white <- function(fill, amount = 0.06) {
  rgb_matrix <- grDevices::col2rgb(fill) / 255
  mixed <- rgb_matrix + (1 - rgb_matrix) * amount
  grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ])
}

segment_text_fill <- function(fill) {
  rgb_matrix <- grDevices::col2rgb(fill)
  luminance <- 0.2126 * rgb_matrix[1, ] + 0.7152 * rgb_matrix[2, ] + 0.0722 * rgb_matrix[3, ]
  ifelse(luminance < 145, "white", "#213547")
}

figure1_compact_count <- function(values) {
  numeric_values <- suppressWarnings(as.numeric(values))
  vapply(
    numeric_values,
    function(value) {
      if (is.na(value)) {
        return("")
      }
      if (abs(value) >= 1e6) {
        return(sprintf("%.1fM", value / 1e6))
      }
      if (abs(value) >= 1e3) {
        return(sprintf("%.1fk", value / 1e3))
      }
      format(round(value), big.mark = ",")
    },
    character(1)
  )
}

load_figure1_registry <- function() {
  read_plot_ready_csv(
    "figure1/figure1_panel_1b_dataset_context_coverage.csv",
    required_cols = c(
      "task",
      "dataset",
      "perturbation_type",
      "cell_line",
      "context_summary_label",
      "wedge_value",
      "dataset_cell_line_count",
      "dataset_target_count",
      "dataset_leaf_weight",
      "context_rank_within_parent",
      "context_folded_bool"
    )
  )
}

load_figure1_support_leaves <- function() {
  if (!is.null(figure1_support_cache)) {
    return(figure1_support_cache)
  }

  data <- load_figure1_registry()
  data$task <- trimws(as.character(data$task))
  data$dataset <- trimws(as.character(data$dataset))
  data$perturbation_type <- trimws(as.character(data$perturbation_type))
  data$cell_line <- trimws(as.character(data$cell_line))
  data$context_summary_label <- trimws(as.character(data$context_summary_label))
  data$weight <- suppressWarnings(as.numeric(data$wedge_value))
  data$dataset_cell_line_count <- suppressWarnings(as.numeric(data$dataset_cell_line_count))
  data$dataset_target_count <- suppressWarnings(as.numeric(data$dataset_target_count))
  data$dataset_leaf_weight <- suppressWarnings(as.numeric(data$dataset_leaf_weight))
  data$context_rank_within_parent <- suppressWarnings(as.integer(data$context_rank_within_parent))
  data$context_folded_bool <- tolower(as.character(data$context_folded_bool)) == "true"
  data <- data[
    nzchar(data$task) &
    nzchar(data$dataset) &
      nzchar(data$perturbation_type) &
      nzchar(data$cell_line) &
      nzchar(data$context_summary_label) &
      !is.na(data$weight) &
      data$weight > 0,
    c(
      "task",
      "dataset",
      "perturbation_type",
      "cell_line",
      "context_summary_label",
      "weight",
      "dataset_cell_line_count",
      "dataset_target_count",
      "dataset_leaf_weight",
      "context_rank_within_parent",
      "context_folded_bool"
    ),
    drop = FALSE
  ]
  data$task <- factor(data$task, levels = c("Task1", "Task2"))
  data$dataset <- factor(data$dataset, levels = c("LINCS", "scPerturb"))
  data$perturbation_type <- factor(data$perturbation_type, levels = c("Chemical", "Genetic"))
  figure1_support_cache <<- data
  data
}

figure1_rank_cell_lines <- function(dataset_data, top_n = 10) {
  cell_totals <- stats::aggregate(weight ~ cell_line, data = dataset_data, FUN = sum)
  cell_totals <- cell_totals[cell_totals$cell_line != "Others", , drop = FALSE]
  cell_totals <- cell_totals[order(-cell_totals$weight, cell_totals$cell_line), , drop = FALSE]
  head(cell_totals$cell_line, top_n)
}

figure1_cross_dataset_cell_line_palette <- function(top_n = 10) {
  leaves <- load_figure1_support_leaves()
  datasets <- c("LINCS", "scPerturb")
  top_by_dataset <- stats::setNames(
    lapply(
      datasets,
      function(dataset_name) {
        dataset_data <- leaves[as.character(leaves$dataset) == dataset_name, , drop = FALSE]
        figure1_rank_cell_lines(dataset_data, top_n = top_n)
      }
    ),
    datasets
  )
  shared_top <- Reduce(intersect, top_by_dataset)
  if (length(shared_top) > 0) {
    shared_weights <- vapply(
      shared_top,
      function(cell_line) {
        sum(leaves$weight[leaves$cell_line == cell_line], na.rm = TRUE)
      },
      numeric(1)
    )
    shared_top <- shared_top[order(-shared_weights, shared_top)]
  }
  dataset_specific <- unlist(
    lapply(
      datasets,
      function(dataset_name) {
        dataset_cells <- setdiff(top_by_dataset[[dataset_name]], shared_top)
        if (length(dataset_cells) == 0) {
          return(character(0))
        }
        dataset_data <- leaves[as.character(leaves$dataset) == dataset_name, , drop = FALSE]
        dataset_weights <- vapply(
          dataset_cells,
          function(cell_line) {
            sum(dataset_data$weight[dataset_data$cell_line == cell_line], na.rm = TRUE)
          },
          numeric(1)
        )
        dataset_cells[order(-dataset_weights, dataset_cells)]
      }
    ),
    use.names = FALSE
  )
  palette_cells <- c(shared_top, dataset_specific)
  if (length(palette_cells) > length(figure1_cell_line_rank_palette)) {
    stop(
      sprintf(
        "Figure 1B requires %d cell-line colors but only %d are defined.",
        length(palette_cells),
        length(figure1_cell_line_rank_palette)
      ),
      call. = FALSE
    )
  }
  list(
    palette_map = c(
      stats::setNames(figure1_cell_line_rank_palette[seq_along(palette_cells)], palette_cells),
      Others = figure1_others_fill()
    ),
    top_by_dataset = top_by_dataset,
    shared_top = shared_top
  )
}

build_figure1_dataset_payload <- function(dataset_name, top_n = 10) {
  leaves <- load_figure1_support_leaves()
  dataset_data <- leaves[as.character(leaves$dataset) == dataset_name, , drop = FALSE]
  if (nrow(dataset_data) == 0) {
    stop(sprintf("Figure 1 has no support rows for dataset %s.", dataset_name), call. = FALSE)
  }

  palette_bundle <- figure1_cross_dataset_cell_line_palette(top_n = top_n)
  palette_map <- palette_bundle$palette_map
  inner <- stats::aggregate(weight ~ task, data = dataset_data, FUN = sum)
  inner <- inner[match(c("Task1", "Task2"), inner$task), , drop = FALSE]
  inner <- inner[!is.na(inner$task), , drop = FALSE]

  middle <- stats::aggregate(weight ~ task + perturbation_type, data = dataset_data, FUN = sum)
  middle <- middle[order(match(middle$task, c("Task1", "Task2")), match(middle$perturbation_type, c("Chemical", "Genetic"))), , drop = FALSE]

  outer <- dataset_data[
    ,
    c(
      "task",
      "perturbation_type",
      "cell_line",
      "context_summary_label",
      "weight",
      "context_rank_within_parent",
      "context_folded_bool"
    ),
    drop = FALSE
  ]
  outer <- outer[
    order(
      match(outer$task, c("Task1", "Task2")),
      match(outer$perturbation_type, c("Chemical", "Genetic")),
      outer$context_rank_within_parent,
      outer$context_summary_label
    ),
    ,
    drop = FALSE
  ]

  list(
    dataset = dataset_name,
    inner = inner,
    middle = middle,
    outer = outer,
    palette_map = palette_map,
    summary = list(
      n_cell_lines = max(dataset_data$dataset_cell_line_count, na.rm = TRUE),
      n_targets = max(dataset_data$dataset_target_count, na.rm = TRUE),
      total_weight = max(dataset_data$dataset_leaf_weight, na.rm = TRUE)
    )
  )
}

compute_child_angles <- function(weights, start_angle, end_angle) {
  total_weight <- sum(weights)
  shares <- c(0, cumsum(weights) / total_weight)
  data.frame(
    start = start_angle + (end_angle - start_angle) * shares[-length(shares)],
    end = start_angle + (end_angle - start_angle) * shares[-1]
  )
}

build_figure1_segments <- function(payload) {
  task_fill <- figure1_task_ring_fill()
  perturb_fill <- figure1_inner_ring_fill()
  segment_cols <- c("segment_id", "ring", "label", "weight", "start", "end", "inner_r", "outer_r", "fill")
  segments <- list()

  inner_angles <- compute_child_angles(
    payload$inner$weight,
    start_angle = pi / 2,
    end_angle = pi / 2 - 2 * pi
  )
  inner <- cbind(payload$inner, inner_angles)
  inner$segment_id <- paste(payload$dataset, "task", inner$task, sep = "::")
  inner$ring <- "task"
  inner$label <- as.character(inner$task)
  inner$inner_r <- 0.82
  inner$outer_r <- 1.42
  inner$fill <- task_fill[inner$label]
  segments[[length(segments) + 1]] <- inner[, segment_cols, drop = FALSE]

  middle_segments <- list()
  outer_segments <- list()

  for (idx in seq_len(nrow(inner))) {
    task_label <- as.character(inner$task[[idx]])
    middle_sub <- payload$middle[payload$middle$task == task_label, , drop = FALSE]
    if (nrow(middle_sub) == 0) {
      next
    }
    middle_angles <- compute_child_angles(middle_sub$weight, inner$start[[idx]], inner$end[[idx]])
    middle_sub <- cbind(middle_sub, middle_angles)
    middle_sub$segment_id <- paste(
      payload$dataset,
      "perturbation",
      task_label,
      middle_sub$perturbation_type,
      sep = "::"
    )
    middle_sub$ring <- "perturbation"
    middle_sub$label <- as.character(middle_sub$perturbation_type)
    middle_sub$inner_r <- 1.50
    middle_sub$outer_r <- 2.06
    middle_sub$fill <- perturb_fill[middle_sub$label]
    middle_segments[[length(middle_segments) + 1]] <- middle_sub[, segment_cols, drop = FALSE]

    for (jdx in seq_len(nrow(middle_sub))) {
      perturbation_label <- as.character(middle_sub$perturbation_type[[jdx]])
      outer_sub <- payload$outer[
        payload$outer$task == task_label &
          payload$outer$perturbation_type == perturbation_label,
        ,
        drop = FALSE
      ]
      if (nrow(outer_sub) == 0) {
        next
      }
      outer_angles <- compute_child_angles(outer_sub$weight, middle_sub$start[[jdx]], middle_sub$end[[jdx]])
      outer_sub <- cbind(outer_sub, outer_angles)
      outer_sub$segment_id <- paste(
        payload$dataset,
        "context",
        task_label,
        perturbation_label,
        outer_sub$cell_line,
        sep = "::"
      )
      outer_sub$ring <- "context"
      outer_sub$label <- outer_sub$context_summary_label
      outer_sub$inner_r <- 2.14
      outer_sub$outer_r <- 3.42
      base_fill <- unname(payload$palette_map[outer_sub$cell_line])
      base_fill[is.na(base_fill)] <- figure1_others_fill()
      outer_sub$fill <- ifelse(
        outer_sub$cell_line == "Others" | outer_sub$context_folded_bool,
        figure1_outer_others_fill(),
        mix_with_white(base_fill, amount = 0.04)
      )
      outer_segments[[length(outer_segments) + 1]] <- outer_sub[, segment_cols, drop = FALSE]
    }
  }

  if (length(middle_segments) > 0) {
    segments[[length(segments) + 1]] <- do.call(rbind, middle_segments)
  }
  if (length(outer_segments) > 0) {
    segments[[length(segments) + 1]] <- do.call(rbind, outer_segments)
  }
  do.call(rbind, segments)
}

segment_polygon <- function(segment_row) {
  angle_span <- abs(segment_row[["end"]] - segment_row[["start"]])
  n_points <- max(16, ceiling(angle_span / (2 * pi) * 220))
  theta <- seq(segment_row[["start"]], segment_row[["end"]], length.out = n_points)
  data.frame(
    segment_id = segment_row[["segment_id"]],
    ring = segment_row[["ring"]],
    fill = segment_row[["fill"]],
    x = c(segment_row[["outer_r"]] * cos(theta), rev(segment_row[["inner_r"]] * cos(theta))),
    y = c(segment_row[["outer_r"]] * sin(theta), rev(segment_row[["inner_r"]] * sin(theta))),
    stringsAsFactors = FALSE
  )
}

build_figure1_polygons <- function(segments) {
  do.call(
    rbind,
    lapply(seq_len(nrow(segments)), function(idx) segment_polygon(segments[idx, , drop = FALSE]))
  )
}

circle_polygon <- function(radius, n = 360) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = radius * cos(theta),
    y = radius * sin(theta)
  )
}

figure1_label_data <- function(segments) {
  labels <- segments
  labels$radius_mid <- (labels$inner_r + labels$outer_r) / 2
  labels$angle_mid <- (labels$start + labels$end) / 2
  labels$arc_length <- abs(labels$end - labels$start) * labels$radius_mid
  labels$label_keep <- FALSE
  labels$label_keep[labels$ring == "task"] <- labels$arc_length[labels$ring == "task"] >= 0.74
  labels$label_keep[labels$ring == "perturbation"] <- labels$arc_length[labels$ring == "perturbation"] >= 0.74
  labels$label_keep[labels$ring == "context"] <- (
    labels$arc_length[labels$ring == "context"] >= 0.82 &
      labels$label[labels$ring == "context"] != "Others"
  )
  labels <- labels[labels$label_keep, , drop = FALSE]
  if (nrow(labels) > 0) {
    target_labels <- labels[labels$ring == "context", , drop = FALSE]
    if (nrow(target_labels) > 12) {
      keep_target_ids <- head(
        target_labels[order(-target_labels$weight, target_labels$label), "segment_id"],
        12
      )
      labels <- labels[
        labels$ring != "context" | labels$segment_id %in% keep_target_ids,
        ,
        drop = FALSE
      ]
    }
  }
  if (nrow(labels) == 0) {
    return(labels)
  }

  labels$x <- labels$radius_mid * cos(labels$angle_mid)
  labels$y <- labels$radius_mid * sin(labels$angle_mid)
  labels$text_angle <- labels$angle_mid * 180 / pi - 90
  flip <- labels$text_angle < -90 | labels$text_angle > 90
  labels$text_angle[flip] <- labels$text_angle[flip] + 180
  labels$text_fill <- segment_text_fill(labels$fill)
  labels
}

finish_sunburst_panel <- function(plot, panel_id) {
  panel_labels(panel_id)
  plot +
    labs(tag = "1B") +
    theme_void(base_family = "Arial") +
    theme(
      plot.tag = element_text(face = "bold", size = 11, color = "#111111"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(6, 6, 6, 6),
      legend.position = "none"
    )
}

build_panel_1B_dataset <- function(dataset_name, panel_id) {
  payload <- build_figure1_dataset_payload(dataset_name)
  segments <- build_figure1_segments(payload)
  polygons <- build_figure1_polygons(segments)
  labels <- figure1_label_data(segments)
  center_disc <- circle_polygon(radius = 0.66)
  dataset_label_size <- if (identical(dataset_name, "scPerturb")) 4.7 else 5.0

  summary_label <- paste(
    sprintf("%s cell lines", format(payload$summary$n_cell_lines, big.mark = ",")),
    sprintf("%s targets", format(payload$summary$n_targets, big.mark = ",")),
    sprintf("%s rows", figure1_compact_count(payload$summary$total_weight)),
    sep = "\n"
  )

  plot <- ggplot() +
    geom_polygon(
      data = polygons,
      aes(x = x, y = y, group = segment_id, fill = fill),
      colour = "white",
      linewidth = 0.42,
      show.legend = FALSE
    ) +
    geom_polygon(
      data = center_disc,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      fill = "#F8FAFC",
      colour = "#E6EBF2",
      linewidth = 0.5
    ) +
    geom_text(
      data = labels[labels$ring == "perturbation", , drop = FALSE],
      aes(x = x, y = y, label = label, angle = text_angle, colour = text_fill),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "bold",
      size = 3.8,
      lineheight = 0.9,
      show.legend = FALSE
    ) +
    geom_text(
      data = labels[labels$ring == "task", , drop = FALSE],
      aes(x = x, y = y, label = label, angle = text_angle, colour = text_fill),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "bold",
      size = 3.0,
      lineheight = 0.9,
      show.legend = FALSE
    ) +
    geom_text(
      data = labels[labels$ring == "context", , drop = FALSE],
      aes(x = x, y = y, label = label, angle = text_angle, colour = text_fill),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "plain",
      size = 2.15,
      lineheight = 0.9,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    annotate(
      "text",
      x = 0,
      y = 0.34,
      label = dataset_name,
      family = "Arial",
      fontface = "bold",
      size = dataset_label_size,
      colour = "#1F2D3A"
    ) +
    annotate(
      "text",
      x = 0,
      y = -0.08,
      label = "Center = dataset support summary\nInner / middle / outer = task / perturbation / context",
      family = "Arial",
      fontface = "plain",
      size = 2.18,
      lineheight = 0.98,
      colour = "#5B6B7A"
    ) +
    annotate(
      "text",
      x = 0,
      y = -0.44,
      label = summary_label,
      family = "Arial",
      fontface = "plain",
      size = 2.55,
      lineheight = 1.0,
      colour = "#1F2D3A"
    ) +
    coord_equal(xlim = c(-3.62, 3.62), ylim = c(-3.62, 3.62), clip = "off")

  finish_sunburst_panel(plot, panel_id)
}

build_panel_1B_LINCS <- function() {
  build_panel_1B_dataset("LINCS", "1B_LINCS")
}

build_panel_1B_SCPERTURB <- function() {
  build_panel_1B_dataset("scPerturb", "1B_SCPERTURB")
}

figure1_color_legend_data <- function(dataset_name) {
  payload <- build_figure1_dataset_payload(dataset_name)
  palette_bundle <- figure1_cross_dataset_cell_line_palette()
  top_cell_lines <- palette_bundle$top_by_dataset[[dataset_name]]
  task_rows <- data.frame(
    section = "Task ring",
    label = c("Task1", "Task2"),
    fill = unname(figure1_task_ring_fill()[c("Task1", "Task2")]),
    note = c("Inner ring segment.", "Inner ring segment."),
    stringsAsFactors = FALSE
  )
  cell_line_rows <- data.frame(
    section = "Context hue",
    label = c(top_cell_lines, "Others"),
    fill = unname(payload$palette_map[c(top_cell_lines, "Others")]),
    note = c(
      ifelse(
        top_cell_lines %in% palette_bundle$shared_top,
        "Shared top-10 hue across LINCS and scPerturb.",
        "Dataset-specific top-10 hue."
      ),
      "Outside this dataset's top-10 cell lines."
    ),
    stringsAsFactors = FALSE
  )
  perturbation_rows <- data.frame(
    section = "Perturbation ring",
    label = c("Chemical", "Genetic"),
    fill = unname(figure1_inner_ring_fill()[c("Chemical", "Genetic")]),
    note = c("Middle ring segment.", "Middle ring segment."),
    stringsAsFactors = FALSE
  )
  legend_data <- rbind(task_rows, perturbation_rows, cell_line_rows)
  legend_data$row_id <- factor(rev(seq_len(nrow(legend_data))), levels = rev(seq_len(nrow(legend_data))))
  legend_data
}

build_panel_1B_legend_dataset <- function(dataset_name) {
  data <- figure1_color_legend_data(dataset_name)
  ggplot(data, aes(x = 1, y = row_id)) +
    geom_tile(aes(fill = fill), width = 0.18, height = 0.78, colour = "#D8E1EA", linewidth = 0.28) +
    geom_text(
      aes(x = 1.18, label = label),
      hjust = 0,
      family = "Arial",
      fontface = "bold",
      size = 2.9,
      colour = "#1F2D3A"
    ) +
    geom_text(
      aes(x = 2.35, label = note),
      hjust = 0,
      family = "Arial",
      size = 2.35,
      colour = "#5B6B7A"
    ) +
    scale_fill_identity() +
    facet_grid(section ~ ., scales = "free_y", space = "free_y", switch = "y") +
    coord_cartesian(xlim = c(0.88, 3.18), clip = "off") +
    labs(
      title = paste(dataset_name, "Figure 1B color legend"),
      subtitle = "Shared high-support contexts keep the same hue across LINCS and scPerturb; Others stays gray."
    ) +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = 10.2, colour = "#1F2D3A"),
      plot.subtitle = element_text(size = 8.2, colour = "#526273"),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5, colour = "#1F2D3A"),
      plot.margin = margin(8, 16, 8, 24)
    )
}

figure1_legend_builders <- function() {
  list(
    `1B_LINCS` = function() build_panel_1B_legend_dataset("LINCS"),
    `1B_SCPERTURB` = function() build_panel_1B_legend_dataset("scPerturb")
  )
}

figure1_legend_dimensions <- list(
  `1B_LINCS` = list(width_mm = 156, height_mm = 120),
  `1B_SCPERTURB` = list(width_mm = 156, height_mm = 120)
)

compose_figure1 <- function() {
  compose_manuscript_figure(
    "figure1",
    panels = list(build_panel_1B_LINCS(), build_panel_1B_SCPERTURB()),
    panel_ids = figure1_panel_ids,
    ncol = 2
  )
}

figure1_panel_builders <- function() {
  list(
    `1B_LINCS` = build_panel_1B_LINCS,
    `1B_SCPERTURB` = build_panel_1B_SCPERTURB
  )
}

figure1_panel_dimensions <- list(
  `1B_LINCS` = list(width_mm = 132, height_mm = 132),
  `1B_SCPERTURB` = list(width_mm = 132, height_mm = 132)
)

render_figure1 <- function(
  output_dir = NULL,
  width_mm = NULL,
  height_mm = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE
) {
  plot <- compose_figure1()
  render_composed_figure(
    "figure1",
    plot,
    output_dir = output_dir,
    width_mm = width_mm,
    height_mm = height_mm,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export
  )
}

render_figure1_panels <- function(
  output_dir = NULL,
  transparent_bg = FALSE,
  font_family = "Arial",
  export = TRUE,
  separate_legends = TRUE
) {
  render_panel_set(
    figure_id = "figure1",
    panel_ids = figure1_panel_ids,
    panel_builders = figure1_panel_builders(),
    panel_dimensions = figure1_panel_dimensions,
    output_dir = output_dir,
    transparent_bg = transparent_bg,
    font_family = font_family,
    export = export,
    separate_legends = separate_legends,
    legend_builders = figure1_legend_builders(),
    legend_dimensions = figure1_legend_dimensions
  )
}

if (sys.nframe() == 0) {
  render_figure1()
}
