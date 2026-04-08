metric_display_names <- c(
  mean_cosine_centroid = "Cosine",
  mean_pcc_centroid = "Pearson",
  mean_edist_biascorr = "eDist",
  cosine_centroid = "Cosine",
  pcc_centroid = "Pearson",
  edist_biascorr = "eDist",
  cosine = "Cosine",
  pcc = "Pearson",
  edist = "eDist",
  mean_mrr_corrected = "MRR",
  mean_hit1_corrected = "Hit@1",
  mean_hit5_corrected = "Hit@5",
  mean_hit10_corrected = "Hit@10",
  mrr_corrected = "MRR",
  hit1_corrected = "Hit@1",
  hit5_corrected = "Hit@5",
  hit10_corrected = "Hit@10",
  mrr = "MRR",
  hit1 = "Hit@1",
  hit5 = "Hit@5",
  hit10 = "Hit@10",
  internal_mean = "Internal mean",
  cross_mean = "Cross mean",
  delta_mean = "Delta mean",
  delta_median = "Delta median"
)

analysis_family_display_names <- c(
  group_concordance = "Group concordance",
  retrieval = "Retrieval",
  cross_eligibility = "Cross eligibility"
)

dataset_display_names <- c(
  LINCS = "LINCS",
  scPerturb = "scPerturb",
  LINCS_to_scPerturb = "LINCS -> scPerturb",
  scPerturb_to_LINCS = "scPerturb -> LINCS",
  cross_genetic = "Cross genetic"
)

scope_display_names <- c(
  internal = "Internal",
  cross = "Cross",
  dataset_cell_line = "Dataset / cell line"
)

direction_display_names <- c(
  ALL = "All",
  C2G = "C2G",
  G2C = "G2C"
)

status_display_names <- c(
  materialized = "Materialized",
  available = "Available",
  support_only = "Support only",
  canonical = "Canonical",
  historical_only = "Historical",
  excluded_by_support_gate = "Excluded",
  not_applicable_scope = "Not applicable"
)

fm_model_display_names <- c(
  geneformer = "Geneformer",
  scbert = "scBERT",
  scfoundation = "scFoundation",
  scgpt = "scGPT",
  state = "STATE",
  `tahoe-x1` = "Tahoe-X1",
  tahoex1_3b = "Tahoe-X1 3B",
  uce = "UCE",
  `FM:geneformer` = "Geneformer",
  `FM:scbert` = "scBERT",
  `FM:scfoundation` = "scFoundation",
  `FM:scgpt` = "scGPT",
  `FM:state` = "STATE",
  `FM:tahoe-x1` = "Tahoe-X1",
  `FM:tahoex1_3b` = "Tahoe-X1 3B",
  `FM:uce` = "UCE",
  `FM:pca50` = "FM PCA50",
  `FM:pca100` = "FM PCA100",
  `FM:pca200` = "FM PCA200"
)

fm_model_shape_values <- c(
  Gene = 16,
  Pathway = 17,
  FM = 21,
  Geneformer = 21,
  scBERT = 22,
  scFoundation = 23,
  scGPT = 24,
  STATE = 25,
  `Tahoe-X1` = 7,
  `Tahoe-X1 3B` = 8,
  UCE = 9,
  `FM PCA50` = 10,
  `FM PCA100` = 11,
  `FM PCA200` = 12
)

task2_representation_levels <- c(
  "Gene",
  "Pathway",
  "Geneformer",
  "scBERT",
  "scFoundation",
  "scGPT",
  "STATE",
  "Tahoe-X1",
  "Tahoe-X1 3B",
  "UCE",
  "FM PCA50",
  "FM PCA100",
  "FM PCA200"
)

lookup_display <- function(values, mapping, default = NULL) {
  text <- as.character(values)
  mapped <- unname(mapping[text])
  missing <- is.na(mapped) | !nzchar(mapped)
  mapped[missing] <- if (is.null(default)) text[missing] else default
  mapped
}

metric_display_name <- function(values) {
  lookup_display(values, metric_display_names)
}

analysis_family_display_name <- function(values) {
  lookup_display(values, analysis_family_display_names)
}

dataset_display_name <- function(values) {
  lookup_display(values, dataset_display_names)
}

scope_display_name <- function(values) {
  lookup_display(values, scope_display_names)
}

direction_display_name <- function(values) {
  lookup_display(values, direction_display_names)
}

status_display_name <- function(values) {
  lookup_display(values, status_display_names)
}

representation_detail_name <- function(values) {
  text <- as.character(values)
  out <- lookup_display(text, fm_model_display_names)
  out[text %in% c("Gene", "Pathway", "FM")] <- text[text %in% c("Gene", "Pathway", "FM")]
  out
}

representation_class_name <- function(values) {
  ifelse(values %in% c("Gene", "Pathway", "FM"), values, "FM")
}

compact_value_label <- function(values, digits = 2) {
  numeric_values <- suppressWarnings(as.numeric(values))
  out <- rep("", length(numeric_values))
  keep <- !is.na(numeric_values)
  out[keep & abs(numeric_values) >= 1000] <- format(round(numeric_values[keep & abs(numeric_values) >= 1000]), big.mark = ",")
  out[keep & abs(numeric_values) < 1000] <- sprintf(paste0("%.", digits, "f"), numeric_values[keep & abs(numeric_values) < 1000])
  out
}

truncate_label <- function(values, width = 32) {
  vapply(
    as.character(values),
    function(value) {
      if (is.na(value) || nchar(value) <= width) {
        return(value)
      }
      paste0(substr(value, 1, width - 3), "...")
    },
    character(1)
  )
}

percent_rank_within <- function(data, group_cols, value_col, out_col = "fill_value") {
  out <- data
  if (nrow(out) == 0) {
    out[[out_col]] <- numeric(0)
    return(out)
  }
  keys <- interaction(out[group_cols], drop = TRUE, lex.order = TRUE)
  values <- out[[value_col]]
  out[[out_col]] <- stats::ave(
    values,
    keys,
    FUN = function(group_values) {
      if (all(is.na(group_values))) {
        return(rep(NA_real_, length(group_values)))
      }
      valid_n <- sum(!is.na(group_values))
      base::rank(group_values, ties.method = "average", na.last = "keep") / valid_n
    }
  )
  out
}

signed_effect_value <- function(effect_direction, median_delta) {
  magnitude <- suppressWarnings(as.numeric(median_delta))
  sign_value <- ifelse(
    effect_direction %in% c("group_a_gt_group_b", "positive_trend", "any_level_effect"),
    1,
    ifelse(effect_direction %in% c("group_b_gt_group_a", "negative_trend"), -1, 0)
  )
  sign_value * ifelse(is.na(magnitude), 0, magnitude)
}

significance_marker <- function(bh_q, test_status = NULL) {
  q_values <- suppressWarnings(as.numeric(bh_q))
  out <- rep("", length(q_values))
  if (!is.null(test_status)) {
    out[test_status == "not_tested_low_n"] <- "n<5"
    out[test_status == "not_tested_underpowered"] <- "underpowered"
  }
  out[is.na(q_values)] <- ifelse(nzchar(out[is.na(q_values)]), out[is.na(q_values)], "")
  out[!is.na(q_values) & q_values < 0.001] <- "***"
  out[!is.na(q_values) & q_values >= 0.001 & q_values < 0.01] <- "**"
  out[!is.na(q_values) & q_values >= 0.01 & q_values < 0.05] <- "*"
  out
}

format_q_value <- function(values, digits = 2) {
  numeric_values <- suppressWarnings(as.numeric(values))
  vapply(
    numeric_values,
    function(value) {
      if (is.na(value)) {
        return("")
      }
      if (value < 0.001) {
        return(format(value, scientific = TRUE, digits = digits, trim = TRUE))
      }
      format(round(value, digits), nsmall = digits, trim = TRUE)
    },
    character(1)
  )
}

comparison_annotation_label <- function(bh_q, test_status = NULL, n_test_units = NULL, digits = 2) {
  q_values <- suppressWarnings(as.numeric(bh_q))
  status <- if (is.null(test_status)) rep(NA_character_, length(q_values)) else as.character(test_status)
  n_values <- suppressWarnings(as.numeric(n_test_units))
  markers <- significance_marker(q_values, status)
  q_labels <- format_q_value(q_values, digits = digits)
  out <- rep("", length(q_values))

  tested <- is.na(status) | status == "tested"
  out[tested] <- ifelse(nzchar(markers[tested]), markers[tested], "ns")
  out[tested & nzchar(q_labels)] <- paste0(out[tested & nzchar(q_labels)], "\nq=", q_labels[tested & nzchar(q_labels)])

  low_n <- status == "not_tested_low_n"
  out[low_n] <- ifelse(
    !is.na(n_values[low_n]),
    paste0("low n\nn=", format(round(n_values[low_n]), trim = TRUE)),
    "low n"
  )

  underpowered <- status == "not_tested_underpowered"
  out[underpowered] <- ifelse(
    !is.na(n_values[underpowered]),
    paste0("underpowered\nn=", format(round(n_values[underpowered]), trim = TRUE)),
    "underpowered"
  )

  out
}

format_signed_delta <- function(values, digits = 2) {
  numeric_values <- suppressWarnings(as.numeric(values))
  vapply(
    numeric_values,
    function(value) {
      if (is.na(value)) {
        return("")
      }
      sprintf(paste0("%+.", digits, "f"), value)
    },
    character(1)
  )
}

representation_short_name <- function(values) {
  labels <- representation_detail_name(values)
  out <- labels
  out[labels == "Gene"] <- "G"
  out[labels == "Pathway"] <- "P"
  out
}

comparison_stats_label <- function(
  bh_q,
  test_status = NULL,
  n_test_units = NULL,
  median_delta = NULL,
  digits = 2
) {
  q_values <- suppressWarnings(as.numeric(bh_q))
  status <- if (is.null(test_status)) rep(NA_character_, length(q_values)) else as.character(test_status)
  n_values <- suppressWarnings(as.numeric(n_test_units))
  delta_values <- if (is.null(median_delta)) rep(NA_real_, length(q_values)) else suppressWarnings(as.numeric(median_delta))
  q_labels <- format_q_value(q_values, digits = digits)
  delta_labels <- format_signed_delta(delta_values, digits = digits)
  out <- rep("", length(q_values))

  low_n <- status == "not_tested_low_n"
  out[low_n] <- ifelse(
    !is.na(n_values[low_n]),
    paste0("low n\nn=", format(round(n_values[low_n]), trim = TRUE)),
    "low n"
  )

  underpowered <- status == "not_tested_underpowered"
  out[underpowered] <- ifelse(
    !is.na(n_values[underpowered]),
    paste0("underpowered\nn=", format(round(n_values[underpowered]), trim = TRUE)),
    "underpowered"
  )

  tested <- is.na(status) | status == "tested"
  for (index in seq_along(out)) {
    if (!tested[index]) {
      next
    }
    parts <- character()
    if (!is.na(q_values[index])) {
      parts <- c(parts, paste0("q=", q_labels[index]))
    }
    if (!is.na(n_values[index])) {
      parts <- c(parts, paste0("n=", format(round(n_values[index]), trim = TRUE)))
    }
    if (!is.na(delta_values[index])) {
      parts <- c(parts, paste0("Δ=", delta_labels[index]))
    }
    if (!length(parts)) {
      next
    }
    out[index] <- paste(parts, collapse = "\n")
  }

  out
}

comparison_stats_inline_label <- function(
  bh_q,
  test_status = NULL,
  n_test_units = NULL,
  median_delta = NULL,
  prefix = NULL,
  digits = 2
) {
  q_values <- suppressWarnings(as.numeric(bh_q))
  status <- if (is.null(test_status)) rep(NA_character_, length(q_values)) else as.character(test_status)
  n_values <- suppressWarnings(as.numeric(n_test_units))
  delta_values <- if (is.null(median_delta)) rep(NA_real_, length(q_values)) else suppressWarnings(as.numeric(median_delta))
  q_labels <- format_q_value(q_values, digits = digits)
  delta_labels <- format_signed_delta(delta_values, digits = digits)
  prefix_values <- if (is.null(prefix)) rep("", length(q_values)) else as.character(prefix)
  out <- rep("", length(q_values))

  for (index in seq_along(out)) {
    header <- if (nzchar(prefix_values[index])) paste0(prefix_values[index], ": ") else ""
    if (!is.na(status[index]) && status[index] == "not_tested_low_n") {
      count_text <- if (!is.na(n_values[index])) paste0(" (n=", format(round(n_values[index]), trim = TRUE), ")") else ""
      out[index] <- paste0(header, "low n", count_text)
      next
    }
    if (!is.na(status[index]) && status[index] == "not_tested_underpowered") {
      count_text <- if (!is.na(n_values[index])) paste0(" (n=", format(round(n_values[index]), trim = TRUE), ")") else ""
      out[index] <- paste0(header, "underpowered", count_text)
      next
    }

    parts <- character()
    if (!is.na(q_values[index])) {
      parts <- c(parts, paste0("q=", q_labels[index]))
    }
    if (!is.na(n_values[index])) {
      parts <- c(parts, paste0("n=", format(round(n_values[index]), trim = TRUE)))
    }
    if (!is.na(delta_values[index])) {
      parts <- c(parts, paste0("d=", delta_labels[index]))
    }
    if (length(parts)) {
      out[index] <- paste0(header, paste(parts, collapse = " | "))
    }
  }

  out
}

metric_uses_signed_log <- function(metric_name) {
  grepl("edist", as.character(metric_name), ignore.case = TRUE)
}

signed_log10_1p <- function(values) {
  numeric_values <- suppressWarnings(as.numeric(values))
  sign(numeric_values) * log10(1 + abs(numeric_values))
}

metric_plot_value <- function(metric_name, values) {
  ifelse(metric_uses_signed_log(metric_name), signed_log10_1p(values), suppressWarnings(as.numeric(values)))
}

metric_value_label <- function(metric_name, values, digits = 2) {
  numeric_values <- suppressWarnings(as.numeric(values))
  out <- compact_value_label(numeric_values, digits = digits)
  out
}

dataset_local_levels <- function(data, level_col, dataset_col = "dataset") {
  datasets <- unique(data[[dataset_col]])
  local_levels <- lapply(datasets, function(ds) {
    unique(data[[level_col]][data[[dataset_col]] == ds])
  })
  unique(unlist(local_levels))
}

drop_empty_facet_levels <- function(data, facet_col) {
  data[[facet_col]] <- droplevels(factor(data[[facet_col]]))
  data
}

abbreviate_metric_header <- function(values, max_width = 18) {
  vapply(as.character(values), function(value) {
    if (is.na(value) || nchar(value) <= max_width) return(value)
    value <- gsub("Group concordance", "Group", value)
    value <- gsub("C2G retrieval", "C2G retr.", value)
    value <- gsub("Retrieval", "Retr.", value)
    if (nchar(value) <= max_width) return(value)
    paste0(substr(value, 1, max_width - 2), "..")
  }, character(1))
}

density_adjusted_height <- function(
  n_rows,
  base_height_mm = 120,
  row_height_mm = 5.5,
  min_height_mm = 90,
  max_height_mm = 400
) {
  estimated <- base_height_mm + max(0, n_rows - 12) * row_height_mm
  max(min_height_mm, min(max_height_mm, estimated))
}

strip_hierarchy_theme <- function(
  outer_size = 9.2,
  inner_size = 8.2,
  y_angle = 0
) {
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(
      size = outer_size,
      face = "bold",
      margin = ggplot2::margin(3, 3, 4, 3)
    ),
    strip.text.y = ggplot2::element_text(
      size = inner_size,
      face = "bold",
      angle = y_angle,
      margin = ggplot2::margin(3, 4, 3, 3)
    )
  )
}

paired_enrichment_exemplar_plot <- function(
  data,
  panel_id,
  x_label = "log2(observed success rate / background success rate)",
  annotation_primary = NULL,
  annotation_secondary = NULL,
  facet_row = NULL,
  facet_col = NULL
) {
  if (!nrow(data)) {
    stop(sprintf("No paired enrichment rows available for panel %s", panel_id), call. = FALSE)
  }
  panel_key <- if (!is.null(facet_row) && !is.null(facet_col) && facet_row %in% colnames(data) && facet_col %in% colnames(data)) {
    paste(data[[facet_row]], data[[facet_col]], sep = "||")
  } else if (!is.null(facet_col) && facet_col %in% colnames(data)) {
    as.character(data[[facet_col]])
  } else {
    as.character(data$facet_value)
  }
  data$panel_key <- panel_key
  shared_column <- if ("shared_flag" %in% colnames(data)) "shared_flag" else "shared_across_modalities_bool"
  point_data <- rbind(
    data.frame(
      panel_key = data$panel_key,
      facet_value = if ("facet_value" %in% colnames(data)) data$facet_value else NA_character_,
      facet_row = if (!is.null(facet_row) && facet_row %in% colnames(data)) data[[facet_row]] else NA_character_,
      facet_col = if (!is.null(facet_col) && facet_col %in% colnames(data)) data[[facet_col]] else NA_character_,
      row_factor = data$row_factor,
      value = data$gene_enrichment_score,
      representation = "Gene",
      marker = data$gene_significance_label,
      shared_flag = data[[shared_column]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      panel_key = data$panel_key,
      facet_value = if ("facet_value" %in% colnames(data)) data$facet_value else NA_character_,
      facet_row = if (!is.null(facet_row) && facet_row %in% colnames(data)) data[[facet_row]] else NA_character_,
      facet_col = if (!is.null(facet_col) && facet_col %in% colnames(data)) data[[facet_col]] else NA_character_,
      row_factor = data$row_factor,
      value = data$pathway_enrichment_score,
      representation = "Pathway",
      marker = data$pathway_significance_label,
      shared_flag = data[[shared_column]],
      stringsAsFactors = FALSE
    )
  )
  point_data$representation <- factor(point_data$representation, levels = c("Gene", "Pathway"))
  support_cols <- c("panel_key", "row_factor", "support_label")
  if ("facet_value" %in% colnames(data)) support_cols <- c("facet_value", support_cols)
  if (!is.null(facet_row) && facet_row %in% colnames(data)) support_cols <- c(facet_row, support_cols)
  if (!is.null(facet_col) && facet_col %in% colnames(data)) support_cols <- c(facet_col, support_cols)
  support_labels <- unique(data[, support_cols, drop = FALSE])
  row_labels <- unique(data[, c("row_factor", "entity_label"), drop = FALSE])
  label_lookup <- setNames(as.character(row_labels$entity_label), as.character(row_labels$row_factor))

  annotation_primary_data <- NULL
  if (!is.null(annotation_primary) && annotation_primary %in% colnames(data)) {
    primary_cols <- c("panel_key", "row_factor", annotation_primary)
    if ("facet_value" %in% colnames(data)) primary_cols <- c("facet_value", primary_cols)
    if (!is.null(facet_row) && facet_row %in% colnames(data)) primary_cols <- c(facet_row, primary_cols)
    if (!is.null(facet_col) && facet_col %in% colnames(data)) primary_cols <- c(facet_col, primary_cols)
    annotation_primary_data <- unique(data[, primary_cols, drop = FALSE])
    colnames(annotation_primary_data)[length(colnames(annotation_primary_data))] <- "label_value"
    annotation_primary_data <- annotation_primary_data[nzchar(annotation_primary_data$label_value), , drop = FALSE]
  }

  annotation_secondary_data <- NULL
  if (!is.null(annotation_secondary) && annotation_secondary %in% colnames(data)) {
    secondary_cols <- c("panel_key", "row_factor", annotation_secondary)
    if ("facet_value" %in% colnames(data)) secondary_cols <- c("facet_value", secondary_cols)
    if (!is.null(facet_row) && facet_row %in% colnames(data)) secondary_cols <- c(facet_row, secondary_cols)
    if (!is.null(facet_col) && facet_col %in% colnames(data)) secondary_cols <- c(facet_col, secondary_cols)
    annotation_secondary_data <- unique(data[, secondary_cols, drop = FALSE])
    colnames(annotation_secondary_data)[length(colnames(annotation_secondary_data))] <- "label_value"
    annotation_secondary_data <- annotation_secondary_data[nzchar(annotation_secondary_data$label_value), , drop = FALSE]
  }

  range_data <- do.call(
    rbind,
    lapply(
      split(data, data$panel_key, drop = TRUE),
      function(panel_frame) {
        has_primary <- !is.null(annotation_primary) && annotation_primary %in% colnames(panel_frame) &&
          any(nzchar(as.character(panel_frame[[annotation_primary]])))
        has_secondary <- !is.null(annotation_secondary) && annotation_secondary %in% colnames(panel_frame) &&
          any(nzchar(as.character(panel_frame[[annotation_secondary]])))
        x_values <- c(
          suppressWarnings(as.numeric(panel_frame$gene_enrichment_score)),
          suppressWarnings(as.numeric(panel_frame$pathway_enrichment_score))
        )
        x_values <- x_values[is.finite(x_values)]
        if (!length(x_values)) {
          x_values <- c(-1, 1)
        }
        x_min <- min(x_values)
        x_max <- max(x_values)
        x_span <- max(x_max - x_min, 0.8)
        left_pad <- if (has_primary && has_secondary) 0.42 else if (has_primary || has_secondary) 0.30 else 0.18
        support_pad <- 0.16
        data.frame(
          panel_key = unique(panel_frame$panel_key)[1],
          x_low = x_min - x_span * left_pad,
          x_high = x_max + x_span * support_pad,
          primary_x = x_min - x_span * max(left_pad - 0.07, 0.12),
          secondary_x = x_min - x_span * max(left_pad - 0.20, 0.06),
          support_x = x_max + x_span * 0.10,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  segment_data <- merge(data, range_data, by = "panel_key", all.x = TRUE, sort = FALSE)
  point_data <- merge(point_data, range_data, by = "panel_key", all.x = TRUE, sort = FALSE)
  support_labels <- merge(support_labels, range_data, by = "panel_key", all.x = TRUE, sort = FALSE)
  if (!is.null(annotation_primary_data) && nrow(annotation_primary_data)) {
    annotation_primary_data <- merge(annotation_primary_data, range_data, by = "panel_key", all.x = TRUE, sort = FALSE)
  }
  if (!is.null(annotation_secondary_data) && nrow(annotation_secondary_data)) {
    annotation_secondary_data <- merge(annotation_secondary_data, range_data, by = "panel_key", all.x = TRUE, sort = FALSE)
  }
  bounds_data <- unique(support_labels[, c("panel_key", "row_factor", "x_low", "x_high"), drop = FALSE])

  plot <- ggplot() +
    geom_blank(
      data = bounds_data,
      aes(x = x_low, y = row_factor)
    ) +
    geom_blank(
      data = bounds_data,
      aes(x = x_high, y = row_factor)
    ) +
    geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "#A7B3BF") +
    geom_segment(
      data = segment_data,
      aes(x = gene_enrichment_score, xend = pathway_enrichment_score, y = row_factor, yend = row_factor),
      linewidth = 0.48,
      color = "#D8E0E8"
    )

  if (!is.null(annotation_primary_data) && nrow(annotation_primary_data)) {
    plot <- plot +
      geom_text(
        data = annotation_primary_data,
        aes(x = primary_x, y = row_factor, label = label_value),
        inherit.aes = FALSE,
        hjust = 0,
        size = 1.72,
        color = "#526273",
        fontface = "bold"
      )
  }

  if (!is.null(annotation_secondary_data) && nrow(annotation_secondary_data)) {
    plot <- plot +
      geom_text(
        data = annotation_secondary_data,
        aes(x = secondary_x, y = row_factor, label = label_value),
        inherit.aes = FALSE,
        hjust = 0,
        size = 1.64,
        color = "#8A4B08"
      )
  }

  plot <- plot +
    geom_point(
      data = point_data,
      aes(x = value, y = row_factor, color = representation, alpha = shared_flag, shape = representation),
      size = 2.28,
      stroke = 0.34
    ) +
    geom_text(
      data = point_data[!is.na(point_data$marker) & nzchar(point_data$marker), , drop = FALSE],
      aes(x = value, y = row_factor, label = marker, color = representation),
      nudge_y = 0.15,
      size = 2.16,
      show.legend = FALSE
    ) +
    geom_text(
      data = support_labels,
      aes(x = support_x, y = row_factor, label = support_label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 1.72,
      color = "#6B7D8E"
    ) +
    scale_color_representation(drop = FALSE) +
    scale_shape_manual(values = c(Gene = 16, Pathway = 17), drop = FALSE) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.78), guide = "none") +
    scale_y_discrete(labels = function(values) unname(label_lookup[as.character(values)])) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(clip = "off") +
    labs(x = x_label, y = NULL, color = "Representation", shape = "Representation") +
    theme_m2m() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 6.8, lineheight = 0.9),
      strip.text.x = element_text(face = "bold", size = 8.8),
      strip.text.y = element_text(face = "bold", size = 7.8, angle = 0),
      panel.spacing.x = grid::unit(6.5, "mm"),
      panel.spacing.y = grid::unit(5, "mm"),
      legend.position = "bottom",
      plot.margin = margin(8, 24, 6, 48)
    )

  if (!is.null(facet_row) && !is.null(facet_col) && facet_row %in% colnames(data) && facet_col %in% colnames(data)) {
    plot <- plot + facet_grid(stats::as.formula(paste(facet_row, "~", facet_col)), scales = "free", space = "free_y")
  } else if (!is.null(facet_col) && facet_col %in% colnames(data)) {
    plot <- plot + facet_grid(stats::as.formula(paste(". ~", facet_col)), scales = "free", space = "free_y")
  } else {
    plot <- plot + facet_grid(. ~ facet_value, scales = "free", space = "free_y")
  }

  finish_panel_plot(plot, panel_id)
}

panel_caption_text <- function(panel_id, caption_lines = character()) {
  invisible(panel_id)
  invisible(caption_lines)
  NULL
}

finish_panel_plot <- function(
  plot,
  panel_id,
  caption_lines = character(),
  legend_position = "bottom"
) {
  invisible(caption_lines)
  panel_labels(panel_id)
  plot +
    ggplot2::labs(
      tag = panel_id
    ) +
    theme_m2m() +
    ggplot2::theme(
      legend.position = legend_position,
      axis.title.x = ggplot2::element_text(
        color = "#1F2D3A", face = "bold", size = 9.5
      ),
      axis.title.y = ggplot2::element_text(
        color = "#1F2D3A", face = "bold", size = 9.5
      ),
      axis.text.x = ggplot2::element_text(
        color = "#3B5068", margin = ggplot2::margin(t = 2)
      ),
      axis.text.y = ggplot2::element_text(color = "#3B5068"),
      plot.margin = ggplot2::margin(10, 10, 7, 8)
    )
}

fm_shape_scale <- function(...) {
  ggplot2::scale_shape_manual(values = fm_model_shape_values, ...)
}

striped_vline <- function(xintercept = 0) {
  ggplot2::geom_vline(
    xintercept = xintercept,
    linewidth = 0.3,
    linetype = "dashed",
    colour = "#8F8F8F"
  )
}
