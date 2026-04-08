source(file.path("plotting", "R", "submission_refresh", "common.R"))

build_submission_refresh_figure2 <- function() {
  bubble <- submission_refresh_csv(
    file.path("Figures_Submission", "Fig2A_Bubble_Effects.csv"),
    required_cols = c("Metric", "MetricPretty", "CompType", "CompShort", "RelChangeCap", "MinusLog10FDR")
  ) %>%
    mutate(
      CompType = factor(CompType, levels = c("Methods", "Modality", "Domain")),
      MetricPretty = factor(
        MetricPretty,
        levels = rev(unique(MetricPretty[order(match(Metric, c("Cosine", "DEG_PCC", "JaccAbs", "Mean_MRR", "Mean_Success", "NegLogEDist")))]))
      ),
      sig_size = pmin(MinusLog10FDR, 4)
    )

  panel_a <- ggplot(bubble, aes(x = CompShort, y = MetricPretty)) +
    geom_point(aes(size = sig_size, fill = RelChangeCap), shape = 21, colour = "#1F2D3A", stroke = 0.3) +
    submission_signed_scale("Relative change") +
    scale_size_continuous(name = "-log10 BH FDR", range = c(1.5, 7), limits = c(0, 4)) +
    facet_wrap(~CompType, nrow = 1, scales = "free_x") +
    labs(x = NULL, y = NULL) +
    theme_submission_refresh() +
    theme(panel.grid.major = element_blank())

  tracer <- submission_refresh_csv(
    file.path("Task2_Unified", "Step4_CaseStudy_Tracer.csv"),
    required_cols = c("Performance_Class", "Scenario", "Track", "View", "MeanSuccess")
  ) %>%
    mutate(
      Performance_Class = factor(
        Performance_Class,
        levels = c("Robust_High", "Protocol_Sensitive", "Robust_Low")
      ),
      Scenario = factor(Scenario, levels = c("A_L1Base", "B_CellDep", "C_TgtSpec")),
      Track = factor(Track, levels = c("Gene", "Pathway")),
      View = factor(View, levels = c("Standard", "Systema"))
    )

  tracer_summary <- tracer %>%
    group_by(Performance_Class, Scenario, Track, View) %>%
    summarise(
      median_success = median(MeanSuccess, na.rm = TRUE),
      q25 = quantile(MeanSuccess, 0.25, na.rm = TRUE),
      q75 = quantile(MeanSuccess, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  panel_b <- ggplot(tracer_summary, aes(x = Scenario, y = median_success, group = View, colour = View, fill = View)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 1.9) +
    facet_grid(Track ~ Performance_Class) +
    scale_colour_manual(values = submission_view_palette()) +
    scale_fill_manual(values = submission_view_palette()) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    scale_x_discrete(labels = c(A_L1Base = "L1", B_CellDep = "L2", C_TgtSpec = "L3")) +
    labs(x = "Evaluation level", y = "Median retrieval success", colour = "View", fill = "View") +
    theme_submission_refresh()

  enrich_paths <- c(
    file.path("Task2_Unified", "Step5_Enrichment_Targets_RobustLow.csv"),
    file.path("Task2_Unified", "Step5_Enrichment_Targets_ProtocolSensitive.csv"),
    file.path("Task2_Unified", "Step5_Enrichment_Targets_RobustHigh.csv")
  )
  enrichment <- bind_rows(lapply(enrich_paths, submission_refresh_csv)) %>%
    mutate(
      Target_Class = factor(
        Target_Class,
        levels = c("Robust_High", "Protocol_Sensitive", "Robust_Low")
      )
    ) %>%
    group_by(Target_Class) %>%
    mutate(effect_direction = ifelse(Log2_OR_Stable >= 0, "Positive", "Negative")) %>%
    arrange(Target_Class, effect_direction, desc(abs(Log2_OR_Stable)), Feature, .by_group = TRUE) %>%
    group_by(Target_Class, effect_direction) %>%
    slice_head(n = 4) %>%
    ungroup() %>%
    mutate(
      Feature = reorder(paste(Feature, Target_Class, sep = "___"), Log2_OR_Stable)
    )

  panel_c <- ggplot(enrichment, aes(x = Log2_OR_Stable, y = Feature)) +
    geom_segment(aes(x = 0, xend = Log2_OR_Stable, yend = Feature), colour = submission_refresh_palette()$neutral_line, linewidth = 0.5) +
    geom_point(aes(size = pmin(NegLog10_P, 6), fill = Target_Class), shape = 21, colour = "#1F2D3A", stroke = 0.3) +
    geom_vline(xintercept = 0, colour = submission_refresh_palette()$neutral_line, linewidth = 0.4, linetype = "dashed") +
    scale_fill_manual(values = submission_class_palette(), guide = "none") +
    scale_size_continuous(name = "-log10 P", range = c(1.7, 6.2)) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    facet_wrap(~Target_Class, scales = "free_y") +
    labs(x = "Stabilized log2 odds ratio", y = NULL) +
    theme_submission_refresh() +
    theme(panel.grid.major.y = element_blank())

  protocol <- submission_refresh_csv(
    file.path("Task2_Unified", "Step5_Protocol_Correlations_ProtocolSensitive.csv"),
    required_cols = c("Cell", "Target", "Dose_Corr", "Time_Corr", "Dose_FDR", "Time_FDR")
  ) %>%
    pivot_longer(
      cols = c(Dose_Corr, Time_Corr),
      names_to = "protocol_metric",
      values_to = "corr_value"
    ) %>%
    mutate(
      fdr_value = ifelse(protocol_metric == "Dose_Corr", Dose_FDR, Time_FDR),
      protocol_type = protocol_label_submission(protocol_metric)
    ) %>%
    filter(!is.na(corr_value), !is.na(fdr_value)) %>%
    mutate(
      signed_score = sign(corr_value) * pmin(-log10(pmax(fdr_value, 1e-12)), 5),
      context_label = paste(Target, "@", Cell)
    ) %>%
    group_by(protocol_type) %>%
    arrange(desc(abs(signed_score)), context_label, .by_group = TRUE) %>%
    slice_head(n = 18) %>%
    ungroup() %>%
    mutate(
      context_label = factor(
        context_label,
        levels = rev(unique(context_label[order(protocol_type, signed_score)]))
      )
    )

  panel_d <- ggplot(protocol, aes(x = signed_score, y = context_label, colour = signed_score)) +
    geom_vline(xintercept = 0, colour = submission_refresh_palette()$neutral_line, linewidth = 0.4, linetype = "dashed") +
    geom_segment(aes(x = 0, xend = signed_score, yend = context_label), colour = submission_refresh_palette()$neutral_line, linewidth = 0.45) +
    geom_point(size = 2.3) +
    submission_signed_colour_scale("Signed\nsignal") +
    facet_wrap(~protocol_type, scales = "free_y") +
    labs(x = "Signed -log10 FDR", y = NULL, colour = "Signed\nsignal") +
    theme_submission_refresh() +
    theme(panel.grid.major.y = element_blank())

  tag_submission_panels((panel_a / panel_b) / (panel_c / panel_d))
}

render_submission_refresh_figure2 <- function() {
  plot <- build_submission_refresh_figure2()
  save_submission_refresh_pdf(plot, "figure2_submission_refresh.pdf", width_mm = 178, height_mm = 220)
}
