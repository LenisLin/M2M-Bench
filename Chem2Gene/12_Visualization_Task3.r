#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ============================================================
# Chem2Gen-Bench — Task3 Visualization (R)
# Inputs: outputs from 11_Task3_Visualization.py
#   - Step1_Task3_Scoreboard_Long.csv
#   - Step0_Task3_Annotated_Pairwise_PerTarget.csv
#   - Step0_Task3_Annotated_Retrieval_PerQuery.csv
#
# Figures:
#   Fig1_Task3_Scoreboard.pdf
#   Fig2_Task3_Improvements_Stats.pdf
#   Fig3_Task3_TargetLevel_Improvements_Heatmap.pdf
#
# Run:
#   Rscript 11_Task3_Visualization.R \
#     --in_dir  /mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Unified \
#     --out_dir /mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Figures \
#     --direction Drug->CRISPR
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggsci)
  library(ggpubr)
  library(ggthemes)
  library(ggrepel)
  library(gghalves)

  library(tidyverse)
  library(dplyr)

  library(scales)
  library(ComplexHeatmap)
  library(grid)
  library(circlize)
})

# =======================
# Paths
# =======================
BASE_RVIS <- "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready"
IN_DIR <- "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Unified"
OUT_DIR <- "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Figures_Submission"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

pick_file <- function(...) {
  cands <- c(...)
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) {
    stop("Missing input file. Tried:\n  - ", paste(cands, collapse = "\n  - "))
  }
  hit
}

# compatible with both old Step* naming and new viz_* naming
FILE_SCORE_LONG <- pick_file(
  file.path(IN_DIR, "Step1_Task3_Scoreboard_Long.csv"),
  file.path(IN_DIR, "viz_scoreboard_long.csv")
)

FILE_PAIR_ANN <- pick_file(
  file.path(IN_DIR, "Step0_Task3_Annotated_Pairwise_PerTarget.csv"),
  file.path(IN_DIR, "annotated_pairwise_per_target.csv")
)

FILE_RET_ANN <- pick_file(
  file.path(IN_DIR, "Step0_Task3_Annotated_Retrieval_PerQuery.csv"),
  file.path(IN_DIR, "annotated_retrieval_per_query.csv")
)

# NEW: produced by updated Script11
FILE_FIG2_LOLLIPOP <- pick_file(
  file.path(IN_DIR, "viz_fig2_lollipop.csv")
)

# =======================
# Helpers
# =======================
save_pdf <- function(path, w = 10, h = 8, expr) {
  has_cairo <- requireNamespace("Cairo", quietly = TRUE)
  if (has_cairo) Cairo::cairo_pdf(path, width = w, height = h) else pdf(path, width = w, height = h)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

scale01_clip <- function(x, q = c(0.05, 0.95), higher_is_better = TRUE) {
  x <- suppressWarnings(as.numeric(x))
  if (all(!is.finite(x))) {
    return(rep(NA_real_, length(x)))
  }
  lo <- quantile(x, probs = q[1], na.rm = TRUE)
  hi <- quantile(x, probs = q[2], na.rm = TRUE)
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
    s <- rep(0.5, length(x))
  } else {
    x2 <- pmin(pmax(x, lo), hi)
    s <- (x2 - lo) / (hi - lo)
  }
  if (!higher_is_better) s <- 1 - s
  s
}

canon_panel_levels <- function(x) {
  # 与 Script11 默认保持一致
  lv <- c("All", "The Cleanest Hits", "The Family Hits", "The Promiscuous Hits")
  factor(x, levels = lv)
}

canon_direction_levels <- function(x) {
  factor(x, levels = c("Drug<->CRISPR", "Drug->CRISPR", "CRISPR->Drug"))
}

short_metric <- function(m) {
  recode(m,
    "Mean_CentroidCosine" = "Cosine",
    "Mean_NegEDist" = "-EDist",
    "Mean_MRR" = "MRR",
    "Mean_Success" = "Success",
    "Median_Rank" = "Rank",
    .default = m
  )
}

tier_order <- c("The Cleanest Hits", "The Family Hits", "The Promiscuous Hits", "Control", "Unknown")


# =======================
# Aesthetics (Mature benchmark palette; unified across scripts)
# =======================
PAL_VIZ <- list(
  view = c(Standard = "#6F99AD", Systema = "#D9B44A"),
  view_bg = c(Standard = "#EEF6FD", Systema = "#FFF7E3"),
  domain = c(SameDomain = "#4C78A8", CrossDomain = "#E45756"),
  tier = c(
    "All" = "#E6E6E6",
    "The Cleanest Hits"    = "#DCEAF7",
    "The Family Hits"      = "#E6E0F2",
    "The Promiscuous Hits" = "#F7DCE0",
    "Control"              = "#ECECEC",
    "Unknown"              = "#ECECEC"
  ),
  direction = c(
    "Drug<->CRISPR" = "#4D4D4D",
    "Drug->CRISPR"  = "#6F99AD",
    "CRISPR->Drug"  = "#8E79B9"
  ),
  heat_seq = c("#F7F7FB", "#C6C4E0", "#5E4FA2"),
  heat_div = c("#4C78A8", "#F7F7F7", "#E45756")
)

# Muted categorical palette for arbitrary Track lists (keeps Gene/Pathway stable when present)
make_muted_tracks <- function(levels_vec) {
  lv <- as.character(levels_vec)
  lv <- unique(lv)
  out <- rep(NA_character_, length(lv))
  names(out) <- lv

  # Prefer semantic anchors if present
  if ("Gene" %in% lv) out["Gene"] <- "#54A24B"
  if ("Pathway" %in% lv) out["Pathway"] <- "#B279A2"

  missing <- names(out)[is.na(out)]
  if (length(missing) > 0) {
    # Pastel HCL palette (stable order)
    n <- length(missing)
    hues <- seq(15, 375, length.out = n + 1)[1:n]
    cols <- grDevices::hcl(h = hues, c = 38, l = 70)
    names(cols) <- missing
    out[missing] <- cols
  }

  out
}


# ============================================================
# Fig1 — Scoreboard heatmap (Scaled_Best)
# ============================================================
message(">>> Fig1: Scoreboard")

# IMPORTANT: canon_direction_levels MUST include Drug<->CRISPR

# tier colors (you can adjust)
tier_cols <- PAL_VIZ$tier

score_long <- readr::read_csv(FILE_SCORE_LONG, show_col_types = FALSE) %>%
  mutate(
    # make Pairwise direction explicit & stable
    Direction = if_else(is.na(Direction) | Direction == "NA", "Drug<->CRISPR", Direction)
  )

df_plot <- score_long %>%
  mutate(
    Panel = canon_panel_levels(Panel),
    MetricGroup = factor(MetricGroup, levels = c("Pairwise", "Retrieval")),
    Direction = canon_direction_levels(Direction),
    View = factor(View, levels = c("Standard", "Systema")),
    Value = suppressWarnings(as.numeric(Scaled_Best))
  ) %>%
  filter(is.finite(Value)) %>%
  mutate(
    ColKey = paste(as.character(Panel), as.character(MetricGroup), as.character(Direction), Metric, sep = " | ")
  ) %>%
  select(View, Track, ColKey, Value) %>%
  distinct()

# track levels: Gene first
track_levels <- df_plot %>%
  distinct(Track) %>%
  pull(Track) %>%
  {
    c("Gene", sort(setdiff(., "Gene")))
  } %>%
  unique()

df_plot <- df_plot %>%
  mutate(Track = factor(Track, levels = track_levels)) %>%
  arrange(View, Track)

df_matrix <- df_plot %>%
  pivot_wider(names_from = ColKey, values_from = Value)

col_keys <- setdiff(colnames(df_matrix), c("View", "Track"))
mat <- as.matrix(df_matrix[, col_keys, drop = FALSE])
rownames(mat) <- as.character(df_matrix$Track)

# ---- column parse + deterministic ordering
col_info <- tibble(ColKey = colnames(mat)) %>%
  separate(ColKey,
    into = c("Panel", "MetricGroup", "Direction", "Metric"),
    sep = " \\| ", fill = "right", remove = FALSE
  ) %>%
  mutate(
    Panel = canon_panel_levels(Panel),
    MetricGroup = factor(MetricGroup, levels = c("Pairwise", "Retrieval")),
    Direction = if_else(is.na(Direction) | Direction == "", "Drug<->CRISPR", Direction),
    Direction = canon_direction_levels(Direction),
    Metric = factor(Metric, levels = c("Mean_CentroidCosine", "Mean_NegEDist", "Mean_MRR", "Mean_Success", "Median_Rank"))
  ) %>%
  arrange(Panel, MetricGroup, Direction, Metric)

mat <- mat[, col_info$ColKey, drop = FALSE]

# labels: short metric
col_labels_vec <- short_metric(as.character(col_info$Metric))

# ---- splits
row_split_vec <- factor(as.character(df_matrix$View), levels = c("Standard", "Systema"))
col_split_vec <- col_info$Panel

# ---- SIMPLE ROW ORDERING (within each view) by All-panel mean score
cols_for_sort <- col_info %>%
  filter(as.character(Panel) == "All") %>%
  pull(ColKey)
if (length(cols_for_sort) < 1) cols_for_sort <- col_info$ColKey

row_score <- rowMeans(mat[, cols_for_sort, drop = FALSE], na.rm = TRUE)

row_df <- tibble(
  idx = seq_len(nrow(mat)),
  View = as.character(df_matrix$View),
  Track = as.character(df_matrix$Track),
  score = row_score
) %>%
  group_by(View) %>%
  arrange(desc(score), Track, .by_group = TRUE) %>%
  ungroup()

ord <- row_df$idx
mat <- mat[ord, , drop = FALSE]
df_matrix <- df_matrix[ord, , drop = FALSE]
row_split_vec <- factor(as.character(df_matrix$View), levels = c("Standard", "Systema"))

# ---- colors
track_cols <- make_muted_tracks(track_levels)

dir_cols <- PAL_VIZ$direction

# ---- TOP ANNOTATION: Tier block with labels + Direction stripe
tier_levels <- levels(col_split_vec)
tier_fill <- tier_cols[tier_levels]
tier_fill[is.na(tier_fill)] <- "#E6E6E6"

dir_vec <- as.character(col_info$Direction)
dir_vec[is.na(dir_vec) | dir_vec == ""] <- "Drug<->CRISPR"

top_ha <- HeatmapAnnotation(
  Tier = anno_block(
    gp = gpar(fill = tier_fill, col = "white"),
    labels = tier_levels,
    labels_gp = gpar(fontsize = 9, fontface = "bold")
  ),
  Direction = dir_vec,
  col = list(Direction = dir_cols),
  show_annotation_name = FALSE,
  annotation_height = unit.c(unit(5, "mm"), unit(2.2, "mm"))
)

left_ha <- rowAnnotation(
  Track = as.character(df_matrix$Track),
  col = list(Track = track_cols),
  show_annotation_name = FALSE,
  width = unit(2.8, "mm")
)

# ---- bar cell_fun (same as yours)
cell_bar_fun <- function(j, i, x, y, width, height, fill) {
  v <- mat[i, j]
  grid.rect(x, y,
    width = width * 0.92, height = height * 0.78,
    gp = gpar(fill = "#F2F2F2", col = NA)
  )
  if (!is.na(v)) {
    v_clamp <- max(min(v, 1.0), 0.03)
    bar_w <- (width * 0.92) * v_clamp
    this_track <- rownames(mat)[i]
    bar_fill <- track_cols[[this_track]]
    grid.rect(x - (width * 0.92) / 2 + bar_w / 2, y,
      width = bar_w, height = height * 0.78,
      gp = gpar(fill = bar_fill, col = NA, alpha = 0.75)
    )
  }
}

ht <- Heatmap(
  mat,
  name = "Scaled",
  col = NULL,
  rect_gp = gpar(type = "none"),
  cell_fun = cell_bar_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split_vec,
  column_split = col_split_vec,
  column_title = rep("", length(levels(col_split_vec))), # avoid duplicated tier text
  left_annotation = left_ha,
  top_annotation = top_ha,
  row_names_gp = gpar(fontsize = 8, fontfamily = "Helvetica"),
  column_labels = col_labels_vec,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 8, fontfamily = "Helvetica"),
  column_names_rot = 90,
  row_gap = unit(3, "mm"),
  column_gap = unit(1.2, "mm"),
  border = FALSE,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold")
)

out1 <- file.path(OUT_DIR, "FigT3_1_Scoreboard_Bars.pdf")
save_pdf(out1, w = 10, h = 4, {
  draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
})
message("Saved: ", out1)


# ============================================================
# Fig2 — Lollipop improvements (Mean; metric-column panels)
#   - Columns: 3 metrics (Cosine / ED-score / Success)
#   - Rows: methods (Track)
#   - Two big panels (top/bottom):
#       (A) View(Systema-Standard)
#       (B) Standard(Track-Gene)  [drop Gene]
#   - Length = MeanDiff; CI = MeanDiff ± 1.96*StdErrMean
#   - Dot size = -log10(p) (capped)
#   - Colors: reuse track_cols from Fig1
# ============================================================
message(">>> Fig2: Lollipop improvements (MeanDiff; metric-column panels)")

DIR_KEEP <- "CRISPR->Drug"  # retrieval direction to display for Success

fig2 <- readr::read_csv(FILE_FIG2_LOLLIPOP, show_col_types = FALSE) %>%
  mutate(
    ComparePanel = as.character(ComparePanel),
    Metric = as.character(Metric),
    Track = as.character(Track),
    Direction = as.character(Direction),
    MeanDiff = suppressWarnings(as.numeric(MeanDiff)),
    StdErrMean = suppressWarnings(as.numeric(StdErrMean)),
    NegLog10P = suppressWarnings(as.numeric(NegLog10P))
  ) %>%
  # keep only 3 metrics
  filter(Metric %in% c("CentroidCosine", "ED_score", "Success_CRISPR2Drug")) %>%
  # IMPORTANT: only filter Direction for Success; keep Pairwise rows where Direction == "NA"
  filter(Metric != "Success_CRISPR2Drug" | Direction == DIR_KEEP) %>%
  # keep only valid MeanDiff
  filter(is.finite(MeanDiff)) %>%
  mutate(
    ComparePanel = factor(
      ComparePanel,
      levels = c("View(Systema-Standard)", "Standard(Track-Gene)")
    ),
    Metric = factor(
      Metric,
      levels = c("CentroidCosine", "ED_score", "Success_CRISPR2Drug"),
      labels = c("Cosine", "ED-score", "Success")
    ),
    NegLog10P = pmin(NegLog10P, 30),  # cap for readability
    # CI only if StdErrMean is finite
    CI_lo = ifelse(is.finite(StdErrMean), MeanDiff - 1.96 * StdErrMean, NA_real_),
    CI_hi = ifelse(is.finite(StdErrMean), MeanDiff + 1.96 * StdErrMean, NA_real_)
  )

# ---- Ensure we reuse Fig1 track order/colors
if (!exists("track_levels")) {
  track_levels <- fig2 %>% distinct(Track) %>% pull(Track) %>% { c("Gene", sort(setdiff(., "Gene"))) } %>% unique()
}
fig2 <- fig2 %>% mutate(Track = factor(Track, levels = track_levels))

if (!exists("track_cols")) {
  npg <- ggsci::pal_npg("nrc")(10)
  track_cols <- make_muted_tracks(track_levels)
}

plot_metric_columns <- function(panel_name, drop_gene = FALSE) {
  dat <- fig2 %>% filter(as.character(ComparePanel) == panel_name)
  if (drop_gene) dat <- dat %>% filter(as.character(Track) != "Gene")

  ggplot(dat, aes(x = MeanDiff, y = Track)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey40") +
    # CI (95%)
    geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.30, linewidth = 0.55, alpha = 0.85, na.rm = TRUE) +
    # lollipop stick
    geom_segment(aes(x = 0, xend = MeanDiff, yend = Track, color = Track),
                 linewidth = 0.75, alpha = 0.9, show.legend = FALSE) +
    # bubble: significance (size only; alpha suppressed)
    geom_point(aes(size = NegLog10P, alpha = NegLog10P, color = Track),
               stroke = 0, show.legend = TRUE) +
    scale_color_manual(values = track_cols, guide = "none") +
    scale_size_continuous(name = expression(-log[10](p)), range = c(1.6, 5.0)) +
    scale_alpha_continuous(range = c(0.35, 1.0), guide = "none") +
    facet_grid(. ~ Metric, scales = "free_x") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = 10, face = "bold"),
      axis.title = element_blank(),
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    ) +
    ggtitle(panel_name)
}

p2_top <- plot_metric_columns("View(Systema-Standard)", drop_gene = FALSE)
p2_bot <- plot_metric_columns("Standard(Track-Gene)", drop_gene = TRUE)

out2 <- file.path(OUT_DIR, "FigT3_2_Lollipop_Improvements_Mean_ByMetricPanels.pdf")
pdf(out2, width = 8, height = 6)
ggpubr::ggarrange(p2_top, p2_bot, ncol = 1, heights = c(1, 1.05),
                   common.legend = TRUE, legend = "right")
dev.off()
message("Saved: ", out2)


# ============================================================
# Fig3 — Target-level improvements vs Gene across 17 targets
#   Fix: use metric-specific color ranges (Pairwise vs Retrieval)
# ============================================================
message(">>> Fig3: Target-level improvements heatmap (17 targets)")

pair_pt <- readr::read_csv(FILE_PAIR_ANN, show_col_types = FALSE)
ret_pq <- readr::read_csv(FILE_RET_ANN, show_col_types = FALSE) %>%
  mutate(specificity_tier = factor(as.character(specificity_tier), levels = tier_order))

# Target order by tier then by baseline Gene(Standard) centroid cosine (desc)
target_order_df <- pair_pt %>%
  filter(Track == "Gene", View == "Standard") %>%
  select(Target, specificity_tier, centroid_cosine) %>%
  mutate(
    specificity_tier = factor(as.character(specificity_tier), levels = tier_order),
    centroid_cosine = as.numeric(centroid_cosine)
  ) %>%
  arrange(specificity_tier, desc(centroid_cosine)) %>%
  distinct(Target, .keep_all = TRUE)

target_levels <- target_order_df$Target

# --- Pairwise per-target delta vs Gene ---
pw_delta <- pair_pt %>%
  select(Track, View, Target, specificity_tier, centroid_cosine) %>%
  mutate(Metric = "Pairwise_CentroidCosine") %>%
  group_by(View, Target) %>%
  mutate(GeneBase = centroid_cosine[Track == "Gene"][1]) %>%
  ungroup() %>%
  mutate(Delta = as.numeric(centroid_cosine) - as.numeric(GeneBase)) %>%
  select(Metric, View, Track, Target, specificity_tier, Delta)

# --- Retrieval per-target (mean MRR) delta vs Gene ---
ret_target <- ret_pq %>%
  filter(Direction == DIR_KEEP) %>%
  group_by(Track, View, True_Target) %>%
  summarise(
    mean_mrr = mean(as.numeric(MRR), na.rm = TRUE),
    tier = first(as.character(specificity_tier)),
    .groups = "drop"
  ) %>%
  rename(Target = True_Target, specificity_tier = tier) %>%
  mutate(Metric = "Retrieval_MeanMRR")

rt_delta <- ret_target %>%
  group_by(View, Target) %>%
  mutate(GeneBase = mean_mrr[Track == "Gene"][1]) %>%
  ungroup() %>%
  mutate(Delta = as.numeric(mean_mrr) - as.numeric(GeneBase)) %>%
  select(Metric, View, Track, Target, specificity_tier, Delta)

delta_all <- bind_rows(pw_delta, rt_delta) %>%
  mutate(
    specificity_tier = factor(as.character(specificity_tier), levels = tier_order),
    Target = factor(Target, levels = target_levels),
    Track = factor(Track, levels = track_levels),
    View = factor(View, levels = c("Standard", "Systema")),
    Metric = factor(Metric, levels = c("Pairwise_CentroidCosine", "Retrieval_MeanMRR"))
  ) %>%
  filter(Target %in% target_levels) %>%
  filter(as.character(specificity_tier) %in% c("The Cleanest Hits","The Family Hits","The Promiscuous Hits","Unknown","Control"))

# --- helper: symmetric, robust color limit per metric (avoid outliers dominating)
sym_limit_q <- function(x, q = 0.95, min_lim = 1e-6) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(min_lim)
  lim <- as.numeric(stats::quantile(abs(x), probs = q, na.rm = TRUE))
  lim <- max(lim, min_lim)
  lim
}

# metric-specific limits
lim_pw <- sym_limit_q(delta_all %>% filter(Metric == "Pairwise_CentroidCosine") %>% pull(Delta), q = 0.95)
lim_rt <- sym_limit_q(delta_all %>% filter(Metric == "Retrieval_MeanMRR") %>% pull(Delta), q = 0.95)

# ---- plot builder (keeps your tier block layout; just different limits)
plot_metric_heat <- function(dat, metric_title, lim) {
  ggplot(dat, aes(x = Target, y = Track, fill = Delta)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(View ~ specificity_tier, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
      limits = c(-lim, lim),
      oob = scales::squish,
      name = paste0("Emb vs Gene\n(", metric_title, ")")
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 7),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA)
    ) +
    labs(x = NULL, y = NULL)
}

p_pw <- plot_metric_heat(
  delta_all %>% filter(Metric == "Pairwise_CentroidCosine"),
  metric_title = "Cosine",
  lim = lim_pw
)

p_rt <- plot_metric_heat(
  delta_all %>% filter(Metric == "Retrieval_MeanMRR"),
  metric_title = paste0("MRR (", DIR_KEEP, ")"),
  lim = lim_rt
)

# ---- combine (two metric blocks stacked), keep layout otherwise unchanged
out3 <- file.path(OUT_DIR, "Fig3_Task3_TargetLevel_Improvements_Heatmap.pdf")
pdf(out3, width = 8, height = 8)
ggpubr::ggarrange(p_pw, p_rt,ncol = 1,heights = c(1, 1),common.legend = FALSE,align = "v")
dev.off()
message("Saved: ", out3)
