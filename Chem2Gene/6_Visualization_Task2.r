#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# Chem2Gen-Bench - Figure 2 Visualization (v5)
#
# v5 changes vs v4:
#   (1) Fig 2A: move "Methods / Modality / Domain" into top facet strips
#       and shorten x-axis labels
#   (2) Restore boxplot code:
#       - Fig2A_Boxplot_KeyDeltas.pdf (Cosine & Mean success, delta distributions)
#       - Fig2A_Boxplot_Domain.pdf    (Same vs Cross, rank-sum + BH FDR)
#   (3) Fix Fig 2B bug: geom_errorbar missing x aesthetic; make summary layers explicit
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

# =========================
# Config
# =========================
BASE_DIR <- "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready"
OUT_DIR  <- file.path(BASE_DIR, "Figures_Submission")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FILE_TASK1_PAIRWISE_WIDE <- file.path(BASE_DIR, "Task1_Metrics", "Task1_Pairwise_Metrics_Wide.csv")
FILE_STEP1_CTX           <- file.path(BASE_DIR, "Task2_Unified", "Step1_L1_Context_Aggregated.csv")
FILE_STEP2_LABELS        <- file.path(BASE_DIR, "Task2_Unified", "Step2_Context_Labels.csv")
FILE_STEP4_TRACER        <- file.path(BASE_DIR, "Task2_Unified", "Step4_CaseStudy_Tracer.csv")

FILE_ENR_HIGH            <- file.path(BASE_DIR, "Task2_Unified", "Step5_Enrichment_Targets_RobustHigh.csv")
FILE_ENR_SENS            <- file.path(BASE_DIR, "Task2_Unified", "Step5_Enrichment_Targets_ProtocolSensitive.csv")
FILE_ENR_LOW             <- file.path(BASE_DIR, "Task2_Unified", "Step5_Enrichment_Targets_RobustLow.csv")

FILE_PROTOCOL            <- file.path(BASE_DIR, "Task2_Unified", "Step5_Protocol_Correlations_ProtocolSensitive.csv")


# =========================
# Aesthetics (Mature benchmark palette; unified across scripts)
# =========================
PAL <- list(
  view  = c(Standard = "#6F99AD", Systema = "#D9B44A"),
  track = c(Gene = "#54A24B", Pathway = "#B279A2"),
  domain = c(SameDomain = "#4C78A8", CrossDomain = "#E45756"),
  perf = c(
    Robust_High        = "#4FAF9A",
    Protocol_Sensitive = "#D9B44A",
    Robust_Low         = "#E45756",
    Other              = "grey75"
  ),
  sig3 = c(
    "Not significant"      = "grey80",
    "Significant positive" = "#4C78A8",
    "Significant negative" = "#E45756"
  ),
  heat_seq = c("#F7F7FB", "#C6C4E0", "#5E4FA2"),
  heat_div = c("#4C78A8", "#F7F7F7", "#E45756")
)

theme_set(

  theme_classic(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      axis.text = element_text(size = 9, color = "black"),
      strip.background = element_rect(fill = "#ECECEC", color = "grey85"),
      strip.text = element_text(face = "bold", size = 10)
    )
)

# =========================
# Helpers
# =========================
stop_if_missing <- function(path, label = "file") {
  if (is.na(path) || !file.exists(path)) stop("CRITICAL missing ", label, ": ", path, call. = FALSE)
}
read_csv_safe <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}
must_have_cols <- function(df, cols, label = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(sprintf("[%s] Missing required columns: %s", label, paste(miss, collapse = ", ")), call. = FALSE)
}
save_gg <- function(p, filename, w = 7.5, h = 5.5) {
  out <- file.path(OUT_DIR, filename)
  ggsave(out, p, width = w, height = h, device = "pdf", useDingbats = FALSE)
  message("Saved: ", out)
}

paired_wilcox <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 5) return(list(n = length(a), p = NA_real_, med_diff = NA_real_))
  d <- b - a
  p <- tryCatch(suppressWarnings(wilcox.test(d, mu = 0, exact = FALSE)$p.value), error = function(e) NA_real_)
  list(n = length(d), p = p, med_diff = median(d, na.rm = TRUE))
}

unpaired_wilcox <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 5 || length(y) < 5) return(list(nx = length(x), ny = length(y), p = NA_real_, med_diff = NA_real_))
  p <- tryCatch(suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value), error = function(e) NA_real_)
  list(nx = length(x), ny = length(y), p = p, med_diff = median(x, na.rm = TRUE) - median(y, na.rm = TRUE))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
}

# =========================
# Load data
# =========================
message(">>> Loading inputs ...")
stop_if_missing(FILE_TASK1_PAIRWISE_WIDE, "Task1_Pairwise_Metrics_Wide")
stop_if_missing(FILE_STEP1_CTX, "Step1_L1_Context_Aggregated")
stop_if_missing(FILE_STEP2_LABELS, "Step2_Context_Labels")
stop_if_missing(FILE_STEP4_TRACER, "Step4_CaseStudy_Tracer")
stop_if_missing(FILE_ENR_HIGH, "Enrichment_RobustHigh")
stop_if_missing(FILE_ENR_SENS, "Enrichment_ProtocolSensitive")
stop_if_missing(FILE_ENR_LOW,  "Enrichment_RobustLow")
stop_if_missing(FILE_PROTOCOL, "Protocol_Correlations_ProtocolSensitive")

pw_raw     <- read_csv_safe(FILE_TASK1_PAIRWISE_WIDE)
ctx_raw    <- read_csv_safe(FILE_STEP1_CTX)
labels_raw <- read_csv_safe(FILE_STEP2_LABELS)
tracer_raw <- read_csv_safe(FILE_STEP4_TRACER)
en_high    <- read_csv_safe(FILE_ENR_HIGH)
en_sens    <- read_csv_safe(FILE_ENR_SENS)
en_low     <- read_csv_safe(FILE_ENR_LOW)
prot_raw   <- read_csv_safe(FILE_PROTOCOL)

must_have_cols(pw_raw, c("Track","cell_std","target_std","domain_type",
                         "source_db_chem","source_db_gene",
                         "cosine_std","cosine_sys",
                         "deg_pcc_std","deg_pcc_sys",
                         "jaccard_abs_std","jaccard_abs_sys",
                         "log10_edist_std","log10_edist_sys"),
               "Task1_Pairwise_Wide")

must_have_cols(ctx_raw, c("Track","View","DomainType","Cell","Target","Source_Chem","Source_Gene",
                          "Mean_Success","Mean_MRR"),
               "Step1_Context_Aggregated")

must_have_cols(labels_raw, c("Cell","Target","Mean_Success","Peak_Success","Performance_Class"),
               "Step2_Labels")

must_have_cols(tracer_raw, c("Cell","Target","Scenario","Track","View","MeanSuccess","MeanMRR","N"),
               "Step4_Tracer")

must_have_cols(en_high, c("Feature","Target_Class","Log2_OR_Stable","FDR_BH"), "Enrichment_High")
must_have_cols(en_sens, c("Feature","Target_Class","Log2_OR_Stable","FDR_BH"), "Enrichment_Sensitive")
must_have_cols(en_low,  c("Feature","Target_Class","Log2_OR_Stable","FDR_BH"), "Enrichment_Low")

must_have_cols(prot_raw, c("Cell","Target","N","Dose_Corr","Dose_FDR","Time_Corr","Time_FDR"),
               "Protocol")

# ==============================================================================
# Prepare unified metric table for Fig 2A summary (7 metrics)
# ==============================================================================
message(">>> Preparing metrics table for Fig 2A summary ...")

eps <- 1e-9
cap_fdr <- 4
cap_rel <- 2

pw_base <- pw_raw %>%
  transmute(
    Track = as.character(Track),
    Cell = as.character(cell_std),
    Target = as.character(target_std),
    DomainType = as.character(domain_type),
    Source_Chem = as.character(source_db_chem),
    Source_Gene = as.character(source_db_gene),

    Cosine_Std      = as.numeric(cosine_std),
    Cosine_Sys      = as.numeric(cosine_sys),
    DEG_PCC_Std     = as.numeric(deg_pcc_std),
    DEG_PCC_Sys     = as.numeric(deg_pcc_sys),
    JaccAbs_Std     = as.numeric(jaccard_abs_std),
    JaccAbs_Sys     = as.numeric(jaccard_abs_sys),
    NegLogEDist_Std = -as.numeric(log10_edist_std),
    NegLogEDist_Sys = -as.numeric(log10_edist_sys)
  )

pw_std_long <- pw_base %>%
  transmute(Track, View = "Standard", DomainType, Cell, Target, Source_Chem, Source_Gene,
            Cosine = Cosine_Std, DEG_PCC = DEG_PCC_Std, JaccAbs = JaccAbs_Std, NegLogEDist = NegLogEDist_Std)
pw_sys_long <- pw_base %>%
  transmute(Track, View = "Systema", DomainType, Cell, Target, Source_Chem, Source_Gene,
            Cosine = Cosine_Sys, DEG_PCC = DEG_PCC_Sys, JaccAbs = JaccAbs_Sys, NegLogEDist = NegLogEDist_Sys)

ctx_metrics <- ctx_raw %>%
  transmute(
    Track = as.character(Track),
    View  = as.character(View),
    DomainType = as.character(DomainType),
    Cell = as.character(Cell),
    Target = as.character(Target),
    Source_Chem = as.character(Source_Chem),
    Source_Gene = as.character(Source_Gene),
    Mean_Success = as.numeric(Mean_Success),
    Mean_MRR     = as.numeric(Mean_MRR)
  )

metrics_all <- bind_rows(pw_std_long, pw_sys_long) %>%
  left_join(ctx_metrics, by = c("Track","View","DomainType","Cell","Target","Source_Chem","Source_Gene"))

metric_pretty <- c(
  Cosine       = "Pairwise cosine",
  DEG_PCC      = "DEG PCC",
  JaccAbs      = "Jaccard abs",
  NegLogEDist  = "Minus log10 E-dist",
  Mean_Success = "Mean success",
  Mean_MRR     = "Mean MRR"
)

metric_cols <- intersect(names(metric_pretty), names(metrics_all))

df_long <- metrics_all %>%
  pivot_longer(cols = all_of(metric_cols), names_to = "Metric", values_to = "Value") %>%
  mutate(MetricPretty = metric_pretty[Metric] %>% unname()) %>%
  filter(is.finite(Value))

# ==============================================================================
# Fig 2A - Bubble plot (facet top strips)
# ==============================================================================
message(">>> Fig 2A overall - bubble plot with top facets ...")

calc_method <- function(track_name) {
  df_long %>%
    filter(Track == track_name) %>%
    select(Metric, MetricPretty, DomainType, Cell, Target, Source_Chem, Source_Gene, View, Value) %>%
    pivot_wider(names_from = View, values_from = Value) %>%
    filter(is.finite(Standard), is.finite(Systema)) %>%
    group_by(Metric, MetricPretty) %>%
    summarise(
      CompType = "Methods",
      CompShort = track_name,  # Gene / Pathway
      N = n(),
      Med_Denom = median(Standard, na.rm = TRUE),
      Med_Numer = median(Systema,  na.rm = TRUE),
      RelChange = (Med_Numer - Med_Denom) / (abs(Med_Denom) + eps),
      P = paired_wilcox(Standard, Systema)$p,
      .groups = "drop"
    )
}

calc_modality <- function(view_name) {
  df_long %>%
    filter(View == view_name) %>%
    select(Metric, MetricPretty, DomainType, Cell, Target, Source_Chem, Source_Gene, Track, Value) %>%
    pivot_wider(names_from = Track, values_from = Value) %>%
    filter(is.finite(Gene), is.finite(Pathway)) %>%
    group_by(Metric, MetricPretty) %>%
    summarise(
      CompType = "Modality",
      CompShort = view_name,   # Standard / Systema
      N = n(),
      Med_Denom = median(Gene,    na.rm = TRUE),
      Med_Numer = median(Pathway, na.rm = TRUE),
      RelChange = (Med_Numer - Med_Denom) / (abs(Med_Denom) + eps),
      P = paired_wilcox(Gene, Pathway)$p,
      .groups = "drop"
    )
}

calc_domain <- function(track_name, view_name, short_label) {
  df_long %>%
    filter(Track == track_name, View == view_name) %>%
    transmute(Metric, MetricPretty, DomainType, Cell, Target, Value) %>%
    group_by(Metric, MetricPretty, DomainType, Cell, Target) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(Value), DomainType %in% c("SameDomain","CrossDomain")) %>%
    group_by(Metric, MetricPretty) %>%
    summarise(
      CompType = "Domain",
      CompShort = short_label,  # short x tick
      N = n(),
      Med_Denom = median(Value[DomainType == "CrossDomain"], na.rm = TRUE),
      Med_Numer = median(Value[DomainType == "SameDomain"],  na.rm = TRUE),
      RelChange = (Med_Numer - Med_Denom) / (abs(Med_Denom) + eps),
      P = unpaired_wilcox(Value[DomainType == "SameDomain"],
                          Value[DomainType == "CrossDomain"])$p,
      .groups = "drop"
    )
}

res_heat <- bind_rows(
  calc_method("Gene"),
  calc_method("Pathway"),
  calc_modality("Standard"),
  calc_modality("Systema"),
  calc_domain("Gene",   "Standard", "Gene Std"),
  calc_domain("Pathway","Standard", "Pathway Std"),
  calc_domain("Gene",   "Systema",  "Gene Sys")
) %>%
  mutate(
    CompType = factor(CompType, levels = c("Methods","Modality","Domain")),
    CompShort = factor(CompShort, levels = c("Gene","Pathway","Standard","Systema","Gene Std","Pathway Std","Gene Sys")),
    FDR = p.adjust(P, method = "BH"),
    MinusLog10FDR = pmin(-log10(FDR + 1e-300), cap_fdr),
    RelChangeCap = pmax(pmin(RelChange, cap_rel), -cap_rel)
  )

readr::write_csv(res_heat, file.path(OUT_DIR, "Fig2A_Bubble_Effects.csv"))

p2a <- ggplot(res_heat, aes(x = CompShort, y = MetricPretty)) +
  geom_point(aes(size = MinusLog10FDR, color = RelChangeCap), alpha = 0.90) +
  scale_color_gradient2(
    midpoint = 0, low = "#2166AC", mid = "white", high = "#B2182B",
    name = "Relative change (cap 2)"
  ) +
  scale_size_continuous(
    range = c(1.6, 7.2),
    name = "Minus log10 BH FDR (cap 4)"
  ) +
  facet_grid(. ~ CompType, scales = "free_x", space = "free_x") +
  scale_x_discrete(expand = expansion(mult = c(0.10, 0.10))) +
  scale_y_discrete(expand = expansion(mult = c(0.06, 0.06))) +
  labs(
    title = "Fig 2A - Summary across metrics and comparisons",
    subtitle = "Color shows relative change; size shows significance (BH FDR).",
    x = NULL, y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 9),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

save_gg(p2a, "Fig2A_SummaryBubble_FacetTop.pdf", w = 12, h = 6)

# ==============================================================================
# Fig 2A boxplots (restored): key metrics only
#   (A) delta distributions for paired comparisons (Methods + Modality)
# ==============================================================================
message(">>> Fig 2A boxplot supplement - key deltas ...")

KEY_METRICS <- c("Cosine", "Mean_Success")
key_pretty  <- c(Cosine = "Pairwise cosine", Mean_Success = "Mean success")

df_key <- metrics_all %>%
  transmute(
    Track, View, DomainType, Cell, Target, Source_Chem, Source_Gene,
    Cosine = as.numeric(Cosine),
    Mean_Success = as.numeric(Mean_Success)
  )

# Paired deltas: Methods (Sys-Std) within Track; Modality (Pathway-Gene) within View
delta_method <- df_key %>%
  filter(Track %in% c("Gene","Pathway")) %>%
  pivot_longer(cols = all_of(KEY_METRICS), names_to = "Metric", values_to = "Value") %>%
  group_by(Track, DomainType, Cell, Target, Source_Chem, Source_Gene, Metric) %>%
  pivot_wider(names_from = View, values_from = Value) %>%
  filter(is.finite(Standard), is.finite(Systema)) %>%
  transmute(
    CompType = "Methods",
    Comp = Track,
    Metric,
    Delta = Systema - Standard
  )

delta_modality <- df_key %>%
  filter(Track %in% c("Gene","Pathway")) %>%
  pivot_longer(cols = all_of(KEY_METRICS), names_to = "Metric", values_to = "Value") %>%
  group_by(View, DomainType, Cell, Target, Source_Chem, Source_Gene, Metric) %>%
  pivot_wider(names_from = Track, values_from = Value) %>%
  filter(is.finite(Gene), is.finite(Pathway)) %>%
  transmute(
    CompType = "Modality",
    Comp = View,
    Metric,
    Delta = Pathway - Gene
  )

delta_all <- bind_rows(delta_method, delta_modality) %>%
  mutate(
    CompType = factor(CompType, levels = c("Methods","Modality")),
    Comp = factor(Comp, levels = c("Gene","Pathway","Standard","Systema")),
    MetricPretty = key_pretty[Metric] %>% unname()
  ) %>%
  filter(is.finite(Delta))

# stats for annotation
delta_stats <- delta_all %>%
  group_by(CompType, Comp, Metric, MetricPretty) %>%
  summarise(
    n = n(),
    p = tryCatch(suppressWarnings(wilcox.test(Delta, mu = 0, exact = FALSE)$p.value), error = function(e) NA_real_),
    med = median(Delta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Metric) %>%
  mutate(fdr = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(
    label = paste0("n=", n, "  med=", sprintf("%.3f", med), "  FDR=", formatC(fdr, format = "e", digits = 2))
  )

p2a_box_delta <- ggplot(delta_all, aes(x = Comp, y = Delta)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.4) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, fill = "white", color = "grey30") +
  geom_jitter(width = 0.16, height = 0, alpha = 0.35, size = 0.7, color = "grey35") +
  geom_text(
    data = delta_stats,
    aes(x = Comp, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.15, size = 2.7, color = "grey15"
  ) +
  facet_grid(MetricPretty ~ CompType, scales = "free_x", space = "free_x") +
  labs(
    title = "Fig 2A supplement - Key delta distributions",
    subtitle = "Methods: delta = Systema minus Standard. Modality: delta = Pathway minus Gene.",
    x = NULL, y = "Delta"
  ) +
  theme(
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
    plot.margin = margin(8, 8, 8, 8)
  )

save_gg(p2a_box_delta, "Fig2A_Boxplot_KeyDeltas.pdf", w = 9.4, h = 5.8)

# ==============================================================================
# Fig 2A (boxplots) - keep v2 behavior (only Cosine + Mean success + BH FDR)
# ==============================================================================
message(">>> Fig 2A boxplots - key metrics only ...")

KEY_METRICS <- c("Cosine", "Mean_Success")

df_key <- metrics_all %>%
  select(Track, View, DomainType, Cell, Target, Source_Chem, Source_Gene, all_of(KEY_METRICS)) %>%
  mutate(
    Track = as.character(Track),
    View  = as.character(View),
    DomainType = as.character(DomainType),
    PairKey = paste(Track, Cell, Target, DomainType, Source_Chem, Source_Gene, sep = "||")
  ) %>%
  pivot_longer(cols = all_of(KEY_METRICS), names_to = "Metric", values_to = "Value") %>%
  mutate(
    MetricPretty = recode(Metric,
                          Cosine = "Pairwise cosine",
                          Mean_Success = "Mean success")
  ) %>%
  filter(is.finite(Value))

# (1) Method: Systema vs Standard (paired), facet by Metric x Track
meth_dat <- df_key %>%
  select(MetricPretty, Track, DomainType, Cell, Target, Source_Chem, Source_Gene, PairKey, View, Value) %>%
  pivot_wider(names_from = View, values_from = Value) %>%
  filter(is.finite(Standard), is.finite(Systema))

meth_stat <- meth_dat %>%
  group_by(MetricPretty, Track) %>%
  summarise(
    n = n(),
    p = paired_wilcox(Standard, Systema)$p,
    .groups = "drop"
  ) %>%
  mutate(fdr = p.adjust(p, method = "BH"),
         label = paste0("BH FDR=", fmt_p(fdr), "; n=", n))

meth_plot <- meth_dat %>%
  pivot_longer(cols = c(Standard, Systema), names_to = "Group", values_to = "Value") %>%
  mutate(Group = factor(Group, levels = c("Standard","Systema")))

p_meth <- ggplot(meth_plot, aes(x = Group, y = Value)) +
  geom_boxplot(aes(fill = Group), width = 0.55, outlier.shape = NA, alpha = 0.45) +
  geom_jitter(width = 0.15, size = 0.25, alpha = 0.12) +
  scale_fill_manual(values = PAL$view, name = "View") +
  facet_grid(MetricPretty ~ Track, scales = "free_y") +
  geom_text(
    data = meth_stat,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.25, size = 3.0
  ) +
  labs(
    title = "Fig 2A (Key) - Method contrast",
    subtitle = "Paired test used; pairing lines omitted to avoid overplotting.",
    x = NULL, y = NULL
  ) +
  theme(legend.position = "bottom")

save_gg(p_meth, "Fig2A_Box_Method_Cosine_MeanSuccess.pdf", w = 8, h = 6)

# (2) Modality: Pathway vs Gene (paired), facet by Metric x View
mod_dat <- df_key %>%
  select(MetricPretty, View, DomainType, Cell, Target, Source_Chem, Source_Gene, Track, Value) %>%
  pivot_wider(names_from = Track, values_from = Value) %>%
  filter(is.finite(Gene), is.finite(Pathway))

mod_stat <- mod_dat %>%
  group_by(MetricPretty, View) %>%
  summarise(
    n = n(),
    p = paired_wilcox(Gene, Pathway)$p,
    .groups = "drop"
  ) %>%
  mutate(fdr = p.adjust(p, method = "BH"),
         label = paste0("BH FDR=", fmt_p(fdr), "; n=", n))

mod_plot <- mod_dat %>%
  pivot_longer(cols = c(Gene, Pathway), names_to = "Group", values_to = "Value") %>%
  mutate(Group = factor(Group, levels = c("Gene","Pathway")))

p_mod <- ggplot(mod_plot, aes(x = Group, y = Value)) +
  geom_boxplot(aes(fill = Group), width = 0.55, outlier.shape = NA, alpha = 0.45) +
  geom_jitter(width = 0.15, size = 0.25, alpha = 0.12) +
  scale_fill_manual(values = PAL$track, name = "Track") +
  facet_grid(MetricPretty ~ View, scales = "free_y") +
  geom_text(
    data = mod_stat,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.25, size = 3.0
  ) +
  labs(
    title = "Fig 2A (Key) - Modality contrast",
    subtitle = "Paired test used; pairing lines omitted to avoid overplotting.",
    x = NULL, y = NULL
  ) +
  theme(legend.position = "bottom")

save_gg(p_mod, "Fig2A_Box_Modality_Cosine_MeanSuccess.pdf", w = 8, h = 6)

# (3) Domain: SameDomain vs CrossDomain (unpaired), unit Cell+Target
dom_dat <- df_key %>%
  group_by(MetricPretty, Track, View, DomainType, Cell, Target) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(Value)) %>%
  mutate(
    DomainType = factor(DomainType, levels = c("CrossDomain","SameDomain")),
    Slice = paste0(Track, " / ", View)
  )

dom_stat <- dom_dat %>%
  group_by(MetricPretty, Slice) %>%
  summarise(
    n_cross = sum(DomainType == "CrossDomain"),
    n_same  = sum(DomainType == "SameDomain"),
    p = unpaired_wilcox(Value[DomainType == "SameDomain"], Value[DomainType == "CrossDomain"])$p,
    .groups = "drop"
  ) %>%
  mutate(fdr = p.adjust(p, method = "BH"),
         label = paste0("BH FDR=", fmt_p(fdr), "; n_same=", n_same, "; n_cross=", n_cross))

p_dom <- ggplot(dom_dat, aes(x = DomainType, y = Value)) +
  geom_boxplot(aes(fill = DomainType), width = 0.55, outlier.shape = NA, alpha = 0.45) +
  geom_jitter(width = 0.15, size = 0.25, alpha = 0.12) +
  scale_fill_manual(values = c(CrossDomain = PAL$domain[["CrossDomain"]],
                               SameDomain  = PAL$domain[["SameDomain"]]),
                    name = "Domain") +
  facet_grid(MetricPretty ~ Slice, scales = "free_y") +
  geom_text(
    data = dom_stat,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.25, size = 2.9
  ) +
  labs(
    title = "Fig 2A (Key) - Domain contrast",
    subtitle = "Unpaired rank-sum test; unit is Cell and Target.",
    x = NULL, y = NULL
  ) +
  theme(legend.position = "bottom")

save_gg(p_dom, "Fig2A_Box_Domain_Cosine_MeanSuccess.pdf", w = 10, h = 6)

# ==============================================================================
# Fig 2B - Tracer: Full Density Background + Highlighted Examples
# ==============================================================================
message(">>> Fig 2B - Tracer (Full Density + Highlights) ...")
library(dplyr)
library(ggplot2)
library(ggrepel)

# --- 1. 数据准备 (与之前一致) ---
tr_plot <- tracer_raw %>%
  filter(Track == "Gene", View == "Standard") %>%
  mutate(
    ContextID = paste0(Cell, " | ", Target),
    ScenarioLabel = case_when(
      Scenario == "A_L1Base"  ~ "L1 (Exact)",
      Scenario == "B_CellDep" ~ "L2 (Cell)",
      Scenario == "C_TgtSpec" ~ "L3 (Target)",
      TRUE ~ Scenario
    ),
    ScenarioLabel = factor(ScenarioLabel, levels = c("L1 (Exact)", "L2 (Cell)", "L3 (Target)")),
    Performance_Class = factor(Performance_Class, 
                               levels = c("Robust_High", "Protocol_Sensitive", "Robust_Low"))
  ) %>%
  filter(!is.na(ScenarioLabel))

# --- 2. 统计计算 (均值 + Error Bar) ---
stats_df <- tr_plot %>%
  group_by(Performance_Class, ScenarioLabel) %>%
  summarise(
    MedianVal = median(MeanSuccess, na.rm = TRUE),
    MeanVal   = mean(MeanSuccess, na.rm = TRUE),
    Lower     = quantile(MeanSuccess, 0.25, na.rm = TRUE), # IQR 25%
    Upper     = quantile(MeanSuccess, 0.75, na.rm = TRUE), # IQR 75%
    .groups = "drop"
  )

# --- 3. 挑选高亮案例 (每组随机 3 个) ---
set.seed(2024) # 固定种子
highlight_ids <- tr_plot %>%
  select(ContextID, Performance_Class) %>%
  distinct() %>%
  group_by(Performance_Class) %>%
  slice_sample(n = 3) %>% # <--- 这里修改数量
  pull(ContextID)

# 提取高亮线条的数据
tr_highlight <- tr_plot %>% filter(ContextID %in% highlight_ids)

# --- 4. 绘图 ---
p_main <- ggplot() +
  
  # A. [背景层] 所有线条
  #    alpha 设得很低 (0.05 - 0.1)，形成“云雾”效果，展示整体分布密度
  geom_line(data = tr_plot, 
            aes(x = ScenarioLabel, y = MeanSuccess, group = ContextID, color = Performance_Class), 
            linewidth = 0.3, alpha = 0.1) + 
  
  # B. [统计层] Error Bars 和 均值线
  #    颜色设为深色或同色系，确保能压住背景
  geom_errorbar(data = stats_df, 
                aes(x = ScenarioLabel, ymin = Lower, ymax = Upper, color = Performance_Class),
                width = 0.15, linewidth = 0.8, alpha = 0.8) +
  
  geom_line(data = stats_df, 
            aes(x = ScenarioLabel, y = MeanVal, group = Performance_Class, color = Performance_Class), 
            linewidth = 1.5) +
  
  geom_point(data = stats_df, 
             aes(x = ScenarioLabel, y = MeanVal, fill = Performance_Class), 
             size = 3.5, shape = 21, color = "white", stroke = 1) +
  
  # C. [高亮层] 选中的案例重绘
  #    不透明 (alpha=1)，稍粗，为了看清这几条线的具体走势
  geom_line(data = tr_highlight, 
            aes(x = ScenarioLabel, y = MeanSuccess, group = ContextID, color = Performance_Class), 
            linewidth = 0.8, alpha = 1) + 
  
  # D. [标签层] 仅给高亮案例的终点 (L3) 加名字
  geom_text_repel(data = tr_highlight %>% filter(ScenarioLabel == "L3 (Target)"),
                  aes(x = ScenarioLabel, y = MeanSuccess, label = Target, color = Performance_Class),
                  size = 3.5, fontface = "bold", 
                  nudge_x = 0.2,       # 往右移一点
                  direction = "y",     # 尽量在Y轴方向避让
                  hjust = 0,           # 左对齐文本
                  segment.size = 0.3, 
                  show.legend = FALSE) +
  
  # E. 样式调整
  scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.25), limits = c(0, 1.05)) +
  # 统一配色
  scale_color_manual(values = c("Robust_High" = "#008B8B", "Protocol_Sensitive" = "#D55E00", "Robust_Low" = "#56B4E9")) +
  scale_fill_manual(values = c("Robust_High" = "#008B8B", "Protocol_Sensitive" = "#D55E00", "Robust_Low" = "#56B4E9")) +
  
  labs(
    title = "Retrieval Consistency & Generalization",
    subtitle = "Background: All pairs; Thick lines: Class Mean ± IQR; Highlighted: 3 random examples per class.",
    y = "Success Score", x = NULL, 
    color = "Class", fill = "Class"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold"),
    # 右侧留白给标签，防止切断
    plot.margin = margin(t = 10, r = 40, b = 10, l = 10, unit = "pt") 
  )

# --- 5. 保存 ---
save_gg(p_main, "Fig2B_Tracer_FullDensity_Highlighted.pdf", w = 10, h = 7.5)

# ==============================================================================
# Fig 2C - Enrichment
# ==============================================================================
message(">>> Fig 2C - Enrichment ...")

en_all <- bind_rows(
  en_high %>% mutate(Class = "Robust_High"),
  en_sens %>% mutate(Class = "Protocol_Sensitive"),
  en_low  %>% mutate(Class = "Robust_Low")
) %>%
  transmute(
    Class = as.character(Class),
    Feature = as.character(Feature),
    Log2OR = as.numeric(Log2_OR_Stable),
    FDR_BH = as.numeric(FDR_BH)
  ) %>%
  filter(is.finite(Log2OR), is.finite(FDR_BH)) %>%
  mutate(NegLog10FDR = -log10(FDR_BH + 1e-300))

TOPN_ENR <- 15
en_top <- en_all %>%
  group_by(Class) %>%
  arrange(FDR_BH, desc(abs(Log2OR))) %>%
  slice_head(n = TOPN_ENR) %>%
  ungroup() %>%
  mutate(
    Class = factor(Class, levels = c("Robust_Low","Protocol_Sensitive","Robust_High")),
    Feature = fct_reorder(Feature, Log2OR)
  )

p2c <- ggplot(en_top, aes(x = Log2OR, y = Feature)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_point(aes(size = NegLog10FDR, color = Class), alpha = 0.90) +
  facet_wrap(~Class, ncol = 1, scales = "free_y") +
  scale_color_manual(values = PAL$perf, name = "Target class") +
  scale_size_continuous(range = c(1.6, 5.2), name = "Minus log10 FDR") +
  labs(
    title = "Fig 2C - Target enrichment by performance class",
    subtitle = "Log2 odds ratio is reported as in Step5 outputs.",
    x = "Log2 odds ratio", y = NULL
  ) +
  theme(legend.position = "right")

save_gg(p2c, "Fig2C_Enrichment_Targets_Log2OR.pdf", w = 6, h = 8)

# ==============================================================================
# Fig 2D - Protocol tunability: violin + dot + group label box; pos/neg significant colors
# ==============================================================================
message(">>> Fig 2D - Protocol violin + dot ...")

# --- 1. Configuration ---
SIG_FDR  <- 0.05
LABEL_N  <- 10      
CAP_VAL  <- 5       # Strict cap for -log10(FDR)

# --- 2. Data Preparation ---
prot_long <- prot_raw %>%
  transmute(
    Cell   = as.character(Cell),
    Target = as.character(Target),
    Dose_Corr = as.numeric(Dose_Corr),
    Dose_FDR  = as.numeric(Dose_FDR),
    Time_Corr = as.numeric(Time_Corr),
    Time_FDR  = as.numeric(Time_FDR)
  ) %>%
  pivot_longer(
    cols = c(Dose_Corr, Dose_FDR, Time_Corr, Time_FDR),
    names_to = c("Factor", ".value"),
    names_pattern = "(Dose|Time)_(Corr|FDR)"
  ) %>%
  mutate(
    Factor = factor(Factor, levels = c("Dose", "Time")),
    # Strict cap at 5 for size calculations
    NegLog10FDR = pmin(-log10(FDR + 1e-300), CAP_VAL),
    # X coordinate for Miami-style layout
    X = ifelse(Factor == "Dose", -1, 1),
    # Signed significance capped between -5 and 5 for color scaling
    SignedSig = ifelse(FDR < SIG_FDR, sign(Corr) * NegLog10FDR, NA_real_),
    LabelTxt = paste0(Target, " (", Cell, ")")
  ) %>%
  filter(is.finite(Corr), is.finite(FDR))

# Separate layers
dat_insig <- prot_long %>% filter(FDR >= SIG_FDR | is.na(FDR))
dat_sig   <- prot_long %>% filter(FDR < SIG_FDR)

# Top hits for labeling
labs_top <- dat_sig %>%
  group_by(Factor) %>%
  slice_max(order_by = NegLog10FDR, n = LABEL_N) %>%
  ungroup()

# --- 3. Plotting ---
p_final <- ggplot(prot_long, aes(x = X, y = Corr)) +
  # 1. Background Distribution & Zero Line
  geom_hline(yintercept = 0, color = "grey60", size = 0.5) +
  geom_violin(aes(group = Factor), fill = "grey97", color = NA, width = 0.8) +
  
  # 2. Layer 1: Insignificant (Background)
  geom_jitter(data = dat_insig, aes(size = NegLog10FDR), 
              color = "grey88", alpha = 0.4, width = 0.25, seed = 42) +
  
  # 3. Layer 2: Significant (Foreground)
  geom_jitter(data = dat_sig, aes(fill = SignedSig, size = NegLog10FDR), 
              color = "white", stroke = 0.3, shape = 21, alpha = 0.9, 
              width = 0.25, seed = 42) +
  
  # 4. Central Annotation Boxes (Restored and Refined)
  # Dose Box
  annotate("rect", xmin = -1.4, xmax = -0.6, ymin = -0.06, ymax = 0.06, 
           fill = "white", color = "#2C7BB6", size = 0.8, alpha = 0.9) +
  annotate("text", x = -1, y = 0, label = "DOSE", fontface = "bold", 
           color = "#2C7BB6", size = 4.5) +
  # Time Box
  annotate("rect", xmin = 0.6, xmax = 1.4, ymin = -0.06, ymax = 0.06, 
           fill = "white", color = "#D7191C", size = 0.8, alpha = 0.9) +
  annotate("text", x = 1, y = 0, label = "TIME", fontface = "bold", 
           color = "#D7191C", size = 4.5) +
  
  # 5. High-Visibility Labels
  geom_text_repel(data = labs_top, aes(label = LabelTxt, color = SignedSig),
                  size = 3.2, fontface = "bold", 
                  force = 2, box.padding = 0.6,
                  bg.color = "white", bg.r = 0.12, max.overlaps = 50) +
  
  # --- Aesthetics & Theme ---
  scale_fill_gradient2(low = "#2C7BB6", mid = "grey90", high = "#D7191C", 
                       midpoint = 0, limits = c(-5, 5), 
                       name = expression(bold("Signed -log")[10]*"(FDR)")) +
  scale_color_gradient2(low = "#2C7BB6", mid = "grey30", high = "#D7191C", 
                        midpoint = 0, limits = c(-5, 5), guide = "none") +
  scale_size_continuous(range = c(1, 5), limits = c(0, 5), 
                        name = expression(bold("-log")[10]*"(FDR)")) +
  
  scale_x_continuous(limits = c(-1.8, 1.8)) +
  scale_y_continuous(limits = c(-1.1, 1.1), breaks = seq(-1, 1, 0.5)) +
  labs(y = "Spearman Correlation (r)", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10)
  )

save_gg(p_final, "Fig2D_Protocol_ViolinDot_Box_PosNeg_Cap4.pdf", w = 6, h = 8)

message(">>> Figure 2 (v5) completed.")
