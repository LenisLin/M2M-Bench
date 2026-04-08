source(file.path("plotting", "R", "submission_refresh", "common.R"))

build_submission_refresh_figure3 <- function() {
  target_gain <- submission_refresh_csv(
    file.path("Task3_Unified", "viz_target_gain.csv"),
    required_cols = c("Metric", "View", "Target", "Target_Tier", "Track", "Delta_Better")
  ) %>%
    filter(Metric %in% c("CentroidCosine", "Success_CRISPR2Drug")) %>%
    mutate(
      metric_label = metric_label_submission(Metric),
      View = factor(View, levels = c("Standard", "Systema")),
      Target_Tier = factor(
        Target_Tier,
        levels = c("The Cleanest Hits", "The Family Hits", "The Promiscuous Hits")
      )
    ) %>%
    group_by(metric_label, View, Target_Tier, Target, Track) %>%
    summarise(delta_better = mean(Delta_Better, na.rm = TRUE), .groups = "drop") %>%
    group_by(metric_label, View, Target_Tier, Track) %>%
    mutate(track_rank = rank(-delta_better, ties.method = "first")) %>%
    ungroup() %>%
    mutate(
      Track = factor(Track, levels = rev(unique(Track))),
      target_label = factor(
        paste(Target_Tier, Target, sep = "___"),
        levels = unique(paste(Target_Tier, Target, sep = "___"))
      )
    )

  panel_a <- ggplot(target_gain, aes(x = target_label, y = Track, fill = delta_better)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    submission_signed_scale("Improvement vs gene baseline") +
    facet_grid(metric_label ~ View, scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    labs(x = NULL, y = NULL) +
    theme_submission_refresh(base_size = 8.8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  scoreboard <- submission_refresh_csv(
    file.path("Task3_Unified", "viz_scoreboard_long.csv"),
    required_cols = c("Track", "View", "Panel", "Metric", "Scaled_Best", "Rank_Best")
  ) %>%
    mutate(
      metric_label = metric_label_submission(Metric),
      panel_metric = paste(Panel, metric_label, sep = " | "),
      View = factor(View, levels = c("Standard", "Systema")),
      panel_metric = factor(panel_metric, levels = unique(panel_metric)),
      Track = factor(Track, levels = rev(unique(Track)))
    )

  panel_b <- ggplot(scoreboard, aes(x = panel_metric, y = Track)) +
    geom_tile(aes(fill = Scaled_Best), colour = "white", linewidth = 0.25) +
    geom_text(aes(label = Rank_Best), size = 2.3, colour = submission_refresh_palette()$dark_text) +
    scale_fill_gradient(
      low = "#F7F8FA",
      high = submission_refresh_palette()$standard,
      name = "Scaled best"
    ) +
    facet_wrap(~View, ncol = 1) +
    labs(x = NULL, y = NULL) +
    theme_submission_refresh(base_size = 8.7) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  lollipop <- submission_refresh_csv(
    file.path("Task3_Unified", "viz_fig2_lollipop.csv"),
    required_cols = c("ComparePanel", "Metric", "Track", "MeanDiff", "StdErrMean", "NegLog10P")
  ) %>%
    mutate(
      metric_label = metric_label_submission(Metric),
      xmin = MeanDiff - 1.96 * StdErrMean,
      xmax = MeanDiff + 1.96 * StdErrMean,
      ComparePanel = factor(ComparePanel, levels = unique(ComparePanel)),
      Track = factor(Track, levels = rev(unique(Track)))
    )

  panel_c <- ggplot(lollipop, aes(x = MeanDiff, y = Track, colour = NegLog10P)) +
    geom_vline(xintercept = 0, colour = submission_refresh_palette()$neutral_line, linewidth = 0.4, linetype = "dashed") +
    geom_segment(aes(x = xmin, xend = xmax, yend = Track), linewidth = 0.6) +
    geom_point(size = 2.3) +
    scale_colour_gradient(
      low = submission_refresh_palette()$neutral_line,
      high = submission_refresh_palette()$systema,
      name = "-log10 P"
    ) +
    facet_grid(ComparePanel ~ metric_label, scales = "free_x") +
    labs(x = "Mean difference", y = NULL, colour = "-log10 P") +
    theme_submission_refresh() +
    theme(panel.grid.major.y = element_blank())

  tag_submission_panels((panel_a / panel_b) / panel_c)
}

render_submission_refresh_figure3 <- function() {
  plot <- build_submission_refresh_figure3()
  save_submission_refresh_pdf(plot, "figure3_submission_refresh.pdf", width_mm = 178, height_mm = 250)
}
