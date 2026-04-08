source(file.path("plotting", "R", "submission_refresh", "common.R"))

build_submission_refresh_figure1 <- function() {
  composition <- submission_refresh_csv(
    file.path("Data_Description", "Fig1B1_SourcePanels_StackedBars_PairCount.csv"),
    required_cols = c("SourceDB", "Bar", "Segment", "Count", "CountUnit")
  ) %>%
    mutate(
      SourceDB = factor(SourceDB, levels = c("LINCS", "scPerturb")),
      Bar = factor(Bar, levels = c("Modality", "Level"))
    )

  composition_totals <- composition %>%
    group_by(SourceDB, Bar) %>%
    summarise(total_count = sum(Count), .groups = "drop")

  panel_a <- ggplot(composition, aes(x = SourceDB, y = Count, fill = Segment)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    geom_text(
      data = composition_totals,
      aes(x = SourceDB, y = total_count, label = comma(total_count)),
      inherit.aes = FALSE,
      nudge_y = max(composition_totals$total_count) * 0.03,
      size = 2.8,
      fontface = "bold",
      colour = submission_refresh_palette()$dark_text
    ) +
    facet_wrap(~Bar, nrow = 1) +
    scale_fill_manual(
      values = c(
        Chemical = submission_refresh_palette()$standard,
        Genetic = submission_refresh_palette()$protocol_sensitive,
        Level1 = submission_refresh_palette()$systema,
        Level3_only = submission_refresh_palette()$robust_high
      )
    ) +
    scale_y_continuous(labels = comma) +
    labs(x = NULL, y = "Eligible pair count", fill = NULL) +
    theme_submission_refresh() +
    theme(panel.grid.major.x = element_blank())

  hierarchy <- submission_refresh_csv(
    file.path("Data_Description", "Hierarchy_Polar_PairCount.csv"),
    required_cols = c("Modality", "tissue", "cell_std", "Count")
  ) %>%
    mutate(
      Modality = factor(Modality, levels = c("Chemical", "Genetic")),
      tissue_label = clean_label(tissue)
    )

  tissue_summary <- hierarchy %>%
    group_by(Modality, tissue_label) %>%
    summarise(
      pair_count = sum(Count),
      cell_line_count = n_distinct(cell_std),
      .groups = "drop"
    ) %>%
    arrange(Modality, desc(pair_count), tissue_label) %>%
    group_by(Modality) %>%
    slice_head(n = 9) %>%
    ungroup() %>%
    mutate(
      tissue_label = factor(
        tissue_label,
        levels = rev(unique(tissue_label[order(Modality, pair_count)]))
      ),
      label_text = paste0(comma(pair_count), " pairs | ", cell_line_count, " cell lines")
    )

  panel_b <- ggplot(tissue_summary, aes(x = pair_count, y = tissue_label, fill = Modality)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    geom_text(
      aes(label = label_text),
      hjust = 0,
      nudge_x = max(tissue_summary$pair_count) * 0.02,
      size = 2.4,
      colour = submission_refresh_palette()$dark_text
    ) +
    coord_cartesian(clip = "off") +
    facet_wrap(~Modality, scales = "free_y") +
    scale_fill_manual(
      values = c(
        Chemical = submission_refresh_palette()$standard,
        Genetic = submission_refresh_palette()$protocol_sensitive
      ),
      guide = "none"
    ) +
    scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.35))) +
    labs(x = "Pair count", y = NULL) +
    theme_submission_refresh() +
    theme(panel.grid.major.y = element_blank())

  pairwise <- submission_refresh_csv(
    file.path("Task1_Metrics", "Task1_Pairwise_Metrics_Long.csv"),
    required_cols = c("Track", "cell_std", "target_std", "View", "Value", "Metric")
  ) %>%
    filter(Metric == "cosine") %>%
    group_by(Track, View, cell_std, target_std) %>%
    summarise(pairwise_cosine = mean(Value, na.rm = TRUE), .groups = "drop")

  retrieval <- submission_refresh_csv(
    file.path("Task1_Retrieval_MultiScenario", "Task1_Retrieval_MultiScenario_PerQuery.csv"),
    required_cols = c("Track", "View", "Cell", "Target", "Success_Score")
  ) %>%
    group_by(Track, View, Cell, Target) %>%
    summarise(mean_success = mean(Success_Score, na.rm = TRUE), .groups = "drop")

  scatter_data <- inner_join(
    pairwise,
    retrieval,
    by = c("Track" = "Track", "View" = "View", "cell_std" = "Cell", "target_std" = "Target")
  ) %>%
    mutate(
      Track = factor(Track, levels = c("Gene", "Pathway")),
      View = factor(View, levels = c("Standard", "Systema"))
    )

  panel_stats <- scatter_data %>%
    group_by(Track, View) %>%
    summarise(
      rho = suppressWarnings(cor(pairwise_cosine, mean_success, method = "spearman", use = "complete.obs")),
      n = n(),
      x = quantile(pairwise_cosine, 0.06, na.rm = TRUE),
      y = quantile(mean_success, 0.94, na.rm = TRUE),
      label = paste0("rho = ", number(rho, accuracy = 0.01), "\n", "n = ", comma(n)),
      .groups = "drop"
    )

  panel_c <- ggplot(scatter_data, aes(x = pairwise_cosine, y = mean_success, colour = View)) +
    geom_point(alpha = 0.33, size = 1.1) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.7,
      span = 0.85
    ) +
    geom_text(
      data = panel_stats,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.6,
      colour = submission_refresh_palette()$dark_text
    ) +
    facet_grid(Track ~ View) +
    scale_colour_manual(values = submission_view_palette(), guide = "none") +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    labs(
      x = "Mean pairwise cosine across source pairs",
      y = "Mean retrieval success"
    ) +
    theme_submission_refresh()

  tag_submission_panels((panel_a | panel_b) / panel_c)
}

render_submission_refresh_figure1 <- function() {
  plot <- build_submission_refresh_figure1()
  save_submission_refresh_pdf(plot, "figure1_submission_refresh.pdf", width_mm = 178, height_mm = 170)
}
