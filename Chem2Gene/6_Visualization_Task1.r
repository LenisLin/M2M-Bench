#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# Chem2Gen-Bench - Visualization Panels (Figure 1 only)  v3.1
#
# This script generates Figure 1 panels aligned to your Task2 outputs:
#   - Fig 1B1: Source panels stacked bars (optional; Data_Description input)
#   - Fig 1B2: Sunburst Modality -> tissue -> cell line (optional; Data_Description input)
#   - Fig 1C : Context Scoreboard (bar-based; ONLY Standard/Systema views)
#   - Fig 1D : Retrieval phenotype fingerprints (purple heatmap + overall bubble)
#   - Fig 1E : Scatter panels (x panels = View Standard/Systema; y panels = Track Gene/Pathway),
#              includes Spearman r and p-value annotations (ASCII only)
#
# Inputs (required):
#   Task2_Unified/Step1_L1_Context_Aggregated.csv
#   Task2_Unified/Step2_Context_Labels.csv
#   Task2_Unified/Step1_L1_Instance_Tidy.csv   (required for Fig 1D)
#
# Inputs (optional, only for Fig 1B1/1B2):
#   Data_Description/Fig1B1_SourcePanels_StackedBars_PairCount.csv  (or *_RecordCount.csv)
#   Data_Description/Hierarchy_Polar.csv
#
# Outputs:
#   <BASE_DIR>/Figures_Submission/*.pdf
#
# Notes:
# - This script intentionally does NOT generate any UMAP panel (your Task2 has no Step3 UMAP).
# - No special characters (delta, rho, arrows) are used in plot labels/titles.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggsci)
  library(ggthemes)
  library(gghalves)

  library(tidyverse)
  library(dplyr)

  library(scales)
  library(ComplexHeatmap)
  library(grid)
  library(circlize)
})

has_cairo <- requireNamespace("Cairo", quietly = TRUE)

# ================= Configuration & Paths =================
BASE_DIR <- "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready"

OUT_DIR <- file.path(BASE_DIR, "Figures_Submission")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Required inputs (Task2 unified outputs) ----
FILE_STEP1_CTX    <- file.path(BASE_DIR, "Task2_Unified", "Step1_L1_Context_Aggregated.csv")
FILE_STEP2_LABELS <- file.path(BASE_DIR, "Task2_Unified", "Step2_Context_Labels.csv")
FILE_STEP1_INST   <- file.path(BASE_DIR, "Task2_Unified", "Step1_L1_Instance_Tidy.csv")

# ---- Task1 metric inputs (for Fig 1C null marker + extra pairwise metrics) ----
FILE_TASK1_WIDE     <- file.path(BASE_DIR, "Task1_Metrics", "Task1_Pairwise_Metrics_Wide.csv")
FILE_TASK1_NULL_COS <- file.path(BASE_DIR, "Task1_Metrics", "Task1_Pairwise_NullCosine_GlobalSummary.csv")
FILE_TASK1_RETR_NULL <- file.path(BASE_DIR, "Task1_Retrieval_MultiScenario", "Task1_Retrieval_MultiScenario_NullSummary.csv")

# ---- Optional (Data description) ----
FILE_FIG1B1_STACK      <- file.path(BASE_DIR, "Data_Description", "Fig1B1_SourcePanels_StackedBars_PairCount.csv")
FILE_FIG1B1_STACK_REC  <- file.path(BASE_DIR, "Data_Description", "Fig1B1_SourcePanels_StackedBars_RecordCount.csv")
FILE_HIER_POLAR        <- file.path(BASE_DIR, "Data_Description", "Hierarchy_Polar.csv")

# ================= Aesthetics (restore original palette semantics) =================
PAL <- list(
  modality = c(
    Chemical = "#7876b1ff",
    Genetic  = "#6f99adff"
  ),
  view = c(
    Standard = "#7aa6dcff",
    Systema  = "#efc000ff"
  ),
  track = c(
    Gene    = "#59A14F",
    Pathway = "#B07AA1"
  ),
  domain = c(
    SameDomain  = "#2C7BB6",
    CrossDomain = "#D7191C"
  ),
  cls = c(
    Robust_high         = "#1B9E77",
    Protocol_sensitive  = "#E6AB02",
    Intermediate        = "#7570B3",
    Robust_low          = "#D95F02"
  )
)

theme_set(
  theme_classic(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      axis.text = element_text(size = 8, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 9)
    )
)

# ================= Helpers =================
stop_if_missing <- function(path, label = "file") {
  if (is.na(path) || !file.exists(path)) stop("CRITICAL missing ", label, ": ", path, call. = FALSE)
}
read_csv_safe <- function(path) {
  if (is.na(path) || !file.exists(path)) return(tibble())
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}
must_have_cols <- function(df, cols, label = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(sprintf("[%s] Missing required columns: %s", label, paste(miss, collapse = ", ")), call. = FALSE)
}
save_gg <- function(p, filename, w = 7.5, h = 5.5, has_cairo = FALSE) {
  out <- file.path(OUT_DIR, filename)
  if (has_cairo) {
    ggsave(out, p, width = w, height = h, device = Cairo::cairo_pdf)
  } else {
    ggsave(out, p, width = w, height = h, device = "pdf", useDingbats = FALSE)
  }
  message("Saved: ", out)
}

scale01_clip <- function(x, q = c(0.05, 0.95)) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  lo <- as.numeric(quantile(x, q[1], na.rm = TRUE))
  hi <- as.numeric(quantile(x, q[2], na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) return(rep(0.5, length(x)))
  x2 <- pmin(pmax(x, lo), hi)
  (x2 - lo) / (hi - lo)
}

std_key <- function(x) toupper(trimws(as.character(x)))

map_perf_class4 <- function(x) {
  x0 <- as.character(x)
  x1 <- tolower(x0)
  x1 <- gsub("[^a-z0-9]+", "", x1)

  dplyr::case_when(
    is.na(x0) | x0 == "" ~ "Intermediate",
    x1 %in% c("robusthigh","high","consistentlygood","consistently_good") ~ "Robust_high",
    x1 %in% c("robustlow","low","consistentlybad","consistently_bad") ~ "Robust_low",
    x1 %in% c("protocolsensitive","protocol_sensitive","conditionallygood","conditionally_good") ~ "Protocol_sensitive",
    x1 %in% c("intermediate","mid") ~ "Intermediate",
    TRUE ~ "Intermediate"
  )
}

# ==============================================================================
# Load Data (required)
# ==============================================================================
message(">>> Loading inputs...")

stop_if_missing(FILE_STEP1_CTX, "Step1_L1_Context_Aggregated.csv")
stop_if_missing(FILE_STEP2_LABELS, "Step2_Context_Labels.csv")

ctx_raw    <- read_csv_safe(FILE_STEP1_CTX)
labels_raw <- read_csv_safe(FILE_STEP2_LABELS)

must_have_cols(ctx_raw, c(
  "Track","View","DomainType","Cell","Target",
  "Source_Chem","Source_Gene",
  "N_Instances","Mean_Success","Median_Success","Max_Success","Mean_MRR","Pairwise_Cosine"
), "Step1_L1_Context_Aggregated")

must_have_cols(labels_raw, c(
  "Cell","Target","Mean_Success","Peak_Success","N_Instances","Performance_Class"
), "Step2_Context_Labels")

message(sprintf("Loaded: ctx_raw=%s rows | labels_raw=%s rows", comma(nrow(ctx_raw)), comma(nrow(labels_raw))))

# Harmonize + context-level collapse (avoid Source_Chem/Source_Gene multiplicity)
ctx <- ctx_raw %>%
  transmute(
    Track = as.character(Track),
    View  = as.character(View),
    DomainType = as.character(DomainType),
    Cell = as.character(Cell),
    Target = as.character(Target),
    Source_Chem = as.character(Source_Chem),
    Source_Gene = as.character(Source_Gene),
    N_Instances = as.numeric(N_Instances),
    Mean_Success = as.numeric(Mean_Success),
    Median_Success = as.numeric(Median_Success),
    Max_Success = as.numeric(Max_Success),
    Mean_MRR = as.numeric(Mean_MRR),
    Pairwise_Cosine = as.numeric(Pairwise_Cosine),
    CellK = std_key(Cell),
    TargetK = std_key(Target)
  )



# ---- Enrich ctx with extra pairwise metrics for Fig 1C (Task1 Pairwise Wide) ----
stop_if_missing(FILE_TASK1_WIDE, "Task1_Pairwise_Metrics_Wide.csv")
pw_wide <- read_csv_safe(FILE_TASK1_WIDE)

must_have_cols(pw_wide, c(
  "Track","cell_std","target_std","domain_type","source_db_chem","source_db_gene",
  "deg_pcc_std","deg_pcc_sys",
  "jaccard_abs_std","jaccard_abs_sys",
  "log10_edist_std","log10_edist_sys"
), "Task1_Pairwise_Metrics_Wide")

pw_long <- pw_wide %>%
  transmute(
    Track = as.character(Track),
    DomainType = as.character(domain_type),
    CellK = std_key(cell_std),
    TargetK = std_key(target_std),
    Source_Chem = as.character(source_db_chem),
    Source_Gene = as.character(source_db_gene),
    deg_pcc_std = as.numeric(deg_pcc_std),
    deg_pcc_sys = as.numeric(deg_pcc_sys),
    jaccard_abs_std = as.numeric(jaccard_abs_std),
    jaccard_abs_sys = as.numeric(jaccard_abs_sys),
    log10_edist_std = as.numeric(log10_edist_std),
    log10_edist_sys = as.numeric(log10_edist_sys)
  ) %>%
  pivot_longer(
    cols = c(deg_pcc_std, deg_pcc_sys,
             jaccard_abs_std, jaccard_abs_sys,
             log10_edist_std, log10_edist_sys),
    names_to = c("Metric", "Suffix"),
    names_pattern = "^(.*)_(std|sys)$",
    values_to = "Value"
  ) %>%
  mutate(
    View = ifelse(Suffix == "std", "Standard", "Systema")
  ) %>%
  select(Track, View, DomainType, CellK, TargetK, Source_Chem, Source_Gene, Metric, Value) %>%
  pivot_wider(names_from = Metric, values_from = Value)

ctx <- ctx %>%
  left_join(pw_long, by = c("Track","View","DomainType","CellK","TargetK","Source_Chem","Source_Gene"))

labels <- labels_raw %>%
  transmute(
    Cell = as.character(Cell),
    Target = as.character(Target),
    CellK = std_key(Cell),
    TargetK = std_key(Target),
    Mean_Success_Lab = as.numeric(Mean_Success),
    Peak_Success_Lab = as.numeric(Peak_Success),
    N_Instances_Lab  = as.numeric(N_Instances),
    Context_Cosine_Lab = if ("Context_Cosine" %in% names(labels_raw)) as.numeric(Context_Cosine) else NA_real_,
    Performance_Class = as.character(Performance_Class),
    Class4 = map_perf_class4(Performance_Class)
  ) %>%
  distinct(CellK, TargetK, .keep_all = TRUE)

ctx_context <- ctx %>%
  group_by(Track, View, DomainType, CellK, TargetK) %>%
  summarise(
    Cell = dplyr::first(Cell),
    Target = dplyr::first(Target),
    N_Instances = sum(N_Instances, na.rm = TRUE),
    Mean_Success = weighted.mean(Mean_Success, w = pmax(N_Instances, 1), na.rm = TRUE),
    Max_Success  = max(Max_Success, na.rm = TRUE),
    Mean_MRR     = weighted.mean(Mean_MRR, w = pmax(N_Instances, 1), na.rm = TRUE),
    Pairwise_Cosine = weighted.mean(Pairwise_Cosine, w = pmax(N_Instances, 1), na.rm = TRUE),
    deg_pcc = weighted.mean(deg_pcc, w = pmax(N_Instances, 1), na.rm = TRUE),
    jaccard_abs = weighted.mean(jaccard_abs, w = pmax(N_Instances, 1), na.rm = TRUE),
    log10_edist = weighted.mean(log10_edist, w = pmax(N_Instances, 1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(labels %>% select(CellK, TargetK, Class4, Mean_Success_Lab, Peak_Success_Lab, Context_Cosine_Lab),
            by = c("CellK","TargetK")) %>%
  mutate(Class4 = ifelse(is.na(Class4) | Class4 == "", "Intermediate", Class4))

# ==============================================================================
# Fig 1B1 - Source panels stacked bars (optional)
# ==============================================================================
message(">>> Fig 1B1 (optional)...")

plot_fig1b1_stackbars <- function(df, title = "Fig 1B1 - Data composition", subtitle = NULL) {
  must_have_cols(df, c("SourceDB","Bar","Segment","Count"), "Fig1B1_StackedBars")

  df2 <- df %>%
    mutate(
      SourceDB = as.character(SourceDB),
      Bar = as.character(Bar),
      Segment = as.character(Segment),
      Count = as.numeric(Count)
    ) %>%
    filter(is.finite(Count), Count >= 0)

  df2 <- df2 %>%
    mutate(
      SourceDB = factor(SourceDB, levels = c("LINCS","scPerturb")),
      Bar = factor(Bar, levels = c("Modality","Level")),
      Segment = factor(Segment, levels = c("Chemical","Genetic","Level3_only","Level1"))
    )

  df2 <- df2 %>%
    group_by(SourceDB, Bar) %>%
    mutate(
      Total = sum(Count, na.rm = TRUE),
      Pct = ifelse(Total > 0, Count / Total, NA_real_),
      Lab = ifelse(is.finite(Pct) & Pct >= 0.06,
                   paste0(scales::comma(Count), "\n", sprintf("%.1f%%", 100 * Pct)),
                   NA_character_)
    ) %>%
    ungroup()

  pal_fill <- c(
    PAL$modality,
    Level1 = "#3182bdff",
    Level3_only = "grey85"
  )

  x_labs <- c(
    Modality = "Modality",
    Level = "Coverage"
  )

  ylab <- "Count"
  if ("CountUnit" %in% names(df2)) {
    unit <- unique(na.omit(as.character(df2$CountUnit)))
    if (length(unit) == 1) ylab <- paste0("Count (", unit, ")")
  }

  ggplot(df2, aes(x = Bar, y = Count, fill = Segment)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    geom_text(aes(label = Lab),
              position = position_stack(vjust = 0.5),
              size = 2.7, lineheight = 0.92, color = "black", na.rm = TRUE) +
    facet_wrap(~SourceDB, ncol = 1, scales = "free_y") +
    scale_x_discrete(labels = x_labs) +
    scale_fill_manual(values = pal_fill, drop = FALSE) +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_classic(base_size = 10.5) +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "grey30")
    )
}

df_fig1b1 <- read_csv_safe(FILE_FIG1B1_STACK)
if (nrow(df_fig1b1) == 0) df_fig1b1 <- read_csv_safe(FILE_FIG1B1_STACK_REC)

if (nrow(df_fig1b1) > 0 && all(c("SourceDB","Bar","Segment","Count") %in% names(df_fig1b1))) {
  subtxt <- "Panels split by source database; bar1 shows modality; bar2 shows Level1 vs Level3-only coverage"
  p1b1 <- plot_fig1b1_stackbars(df_fig1b1, subtitle = subtxt)
  save_gg(p1b1, "Fig1B1_SourcePanels_StackedBars.pdf", w = 7.2, h = 6.2, has_cairo = has_cairo)
} else {
  message("Skipping Fig 1B1: missing Data_Description input or schema mismatch.")
}

# ==============================================================================
# Fig 1B2 - Sunburst Modality -> tissue -> cell line (optional)
#   Requires Data_Description/Hierarchy_Polar.csv with columns:
#     Modality, tissue, cell_std, Count
# ==============================================================================
message(">>> Fig 1B2 (optional sunburst)...")

collapse_topk <- function(df, group_cols, key_col, count_col = "Count",
                          keep_n = 10, other_label = "Other") {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(.rk = dense_rank(desc(.data[[count_col]]))) %>%
    mutate("{key_col}" := ifelse(.rk <= keep_n, .data[[key_col]], other_label)) %>%
    ungroup() %>%
    group_by(across(all_of(group_cols)), across(all_of(key_col))) %>%
    summarise("{count_col}" := sum(.data[[count_col]], na.rm = TRUE), .groups = "drop")
}

compute_sunburst_rects <- function(df, levels, count_col = "Count",
                                   r0 = 1.2, ring_w = 0.8,
                                   sort_within_parent = TRUE) {
  must_have_cols(df, c(levels, count_col), "sunburst_input")

  df0 <- df %>%
    mutate(across(all_of(levels), ~ as.character(.x))) %>%
    group_by(across(all_of(levels))) %>%
    summarise(Count = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop")

  out_list <- list()
  parent_ranges <- tibble(.root = 1, ymin = 0, ymax = 1)

  for (d in seq_along(levels)) {
    cur_levels <- levels[1:d]
    parent_levels <- if (d == 1) character(0) else levels[1:(d - 1)]
    node_col <- levels[d]

    cur <- df0 %>%
      group_by(across(all_of(cur_levels))) %>%
      summarise(Count = sum(Count), .groups = "drop")

    if (d == 1) {
      cur <- cur %>%
        mutate(.root = 1) %>%
        left_join(parent_ranges, by = ".root") %>%
        group_by(.root) %>%
        mutate(parent_total = sum(Count), frac = Count / parent_total) %>%
        { if (sort_within_parent) arrange(., .root, .data[[node_col]]) else . } %>%
        mutate(
          cum = cumsum(frac),
          ymin_local = lag(cum, default = 0),
          ymax_local = cum,
          ymin_new = ymin + (ymax - ymin) * ymin_local,
          ymax_new = ymin + (ymax - ymin) * ymax_local
        ) %>%
        ungroup() %>%
        transmute(!!!syms(cur_levels), Count = Count, ymin = ymin_new, ymax = ymax_new)
    } else {
      parent_ranges2 <- out_list[[d - 1]] %>%
        select(all_of(parent_levels), ymin, ymax) %>%
        distinct()

      cur <- cur %>%
        left_join(parent_ranges2, by = parent_levels) %>%
        group_by(across(all_of(parent_levels))) %>%
        mutate(parent_total = sum(Count), frac = Count / parent_total) %>%
        { if (sort_within_parent) arrange(., across(all_of(parent_levels)), .data[[node_col]]) else . } %>%
        mutate(
          cum = cumsum(frac),
          ymin_local = lag(cum, default = 0),
          ymax_local = cum,
          ymin_new = ymin + (ymax - ymin) * ymin_local,
          ymax_new = ymin + (ymax - ymin) * ymax_local
        ) %>%
        ungroup() %>%
        transmute(!!!syms(cur_levels), Count = Count, ymin = ymin_new, ymax = ymax_new)
    }

    xmin <- r0 + (d - 1) * ring_w
    xmax <- xmin + ring_w

    out_list[[d]] <- cur %>%
      mutate(
        depth = d,
        xmin = xmin, xmax = xmax,
        y = (ymin + ymax) / 2,
        node = .data[[node_col]],
        frac_arc = (ymax - ymin)
      )
  }
  bind_rows(out_list)
}

build_color_map <- function(nodes) {
  nodes <- unique(as.character(nodes))
  fixed <- c(
    "Chemical" = PAL$modality[["Chemical"]],
    "Genetic"  = PAL$modality[["Genetic"]],
    "Other"    = "grey80"
  )
  dynamic_nodes <- setdiff(nodes, names(fixed))
  if (length(dynamic_nodes) > 0) {
    dyn_pal <- ggsci::pal_igv("default")(length(dynamic_nodes))
    names(dyn_pal) <- dynamic_nodes
    c(fixed, dyn_pal)
  } else {
    fixed
  }
}

plot_sunburst_multi <- function(df, levels, title, subtitle,
                                count_col = "Count",
                                r0 = 1.2, ring_w = 0.8,
                                min_label_frac = 0.03) {
  rects <- compute_sunburst_rects(df, levels = levels, count_col = count_col, r0 = r0, ring_w = ring_w)
  lab <- rects %>%
    mutate(label = paste0(node, "\n", scales::comma(Count))) %>%
    filter(frac_arc >= min_label_frac)
  cols <- build_color_map(rects$node)

  ggplot() +
    geom_rect(
      data = rects,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = node),
      color = "white", linewidth = 0.25
    ) +
    ggrepel::geom_text_repel(
      data = lab,
      aes(x = (xmin + xmax) / 2, y = y, label = label),
      size = 2.4, color = "black",
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.1,
      max.overlaps = Inf
    ) +
    coord_polar(theta = "y") +
    xlim(r0 - 0.2, r0 + ring_w * length(levels) + 0.2) +
    theme_void() +
    labs(title = title, subtitle = subtitle) +
    scale_fill_manual(values = cols, guide = "none")
}

df_hier <- read_csv_safe(FILE_HIER_POLAR)
if (nrow(df_hier) > 0 && all(c("Modality","tissue","cell_std","Count") %in% names(df_hier))) {

  df_tissue <- df_hier %>%
    group_by(Modality, tissue) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    group_by(Modality) %>%
    mutate(.rk = dense_rank(desc(Count)),
           tissue = ifelse(.rk <= 10, tissue, "Other")) %>%
    ungroup() %>%
    group_by(Modality, tissue) %>%
    summarise(Count = sum(Count), .groups = "drop")

  df_cell <- df_hier %>%
    mutate(tissue = as.character(tissue)) %>%
    inner_join(df_tissue %>% select(Modality, tissue), by = c("Modality","tissue")) %>%
    group_by(Modality, tissue, cell_std) %>%
    summarise(Count = sum(Count), .groups = "drop")

  df_cell <- collapse_topk(df_cell, group_cols = c("Modality","tissue"), key_col = "cell_std", keep_n = 8)

  p1b2 <- plot_sunburst_multi(
    df = df_cell,
    levels = c("Modality","tissue","cell_std"),
    title = "Fig 1B2 - Composition by tissue and cell line",
    subtitle = "Modality to tissue to cell line (top tissues and cells; rest collapsed)",
    min_label_frac = 0.02
  )
  save_gg(p1b2, "Fig1B2_Sunburst_Modality_Tissue_CellLine.pdf", w = 9.0, h = 5.2, has_cairo = has_cairo)

} else {
  message("Skipping Fig 1B2: Hierarchy_Polar.csv not found or schema mismatch.")
  message("If you want Fig 1B2, generate Data_Description/Hierarchy_Polar.csv with columns: Modality,tissue,cell_std,Count.")
}

# ==============================================================================
# Fig 1C - Context Scoreboard (bar-based; ONLY Standard/Systema)
# ==============================================================================
message(">>> Fig 1C - Context Scoreboard (bar-based; Std/Sys only)...")

select_contexts_base <- function(domain_name, n_total = 24, n_high = 6, n_low = 6) {
  pool <- ctx_context %>%
    filter(DomainType == domain_name, Track == "Gene", View == "Standard") %>%
    distinct(CellK, TargetK) %>%
    left_join(
      ctx_context %>%
        filter(DomainType == domain_name, Track == "Gene", View == "Standard") %>%
        group_by(CellK, TargetK) %>%
        summarise(N_Instances = sum(N_Instances, na.rm = TRUE), .groups = "drop"),
      by = c("CellK","TargetK")
    ) %>%
    left_join(labels %>% select(CellK, TargetK, Class4, Mean_Success_Lab, Peak_Success_Lab),
              by = c("CellK","TargetK"))

  high <- pool %>%
    filter(Class4 == "Robust_high") %>%
    arrange(desc(Mean_Success_Lab), desc(Peak_Success_Lab)) %>%
    slice_head(n = n_high) %>%
    mutate(Kind = "High")

  low <- pool %>%
    filter(Class4 == "Robust_low") %>%
    arrange(Mean_Success_Lab, Peak_Success_Lab) %>%
    slice_head(n = n_low) %>%
    mutate(Kind = "Low")

  support <- pool %>%
    anti_join(bind_rows(high, low) %>% select(CellK, TargetK), by = c("CellK","TargetK")) %>%
    arrange(desc(N_Instances)) %>%
    slice_head(n = max(n_total - nrow(high) - nrow(low), 0)) %>%
    mutate(Kind = "Common")

  bind_rows(high, low, support) %>%
    mutate(DomainType = domain_name) %>%
    select(CellK, TargetK, DomainType, Kind)
}

sel_base <- bind_rows(select_contexts_base("SameDomain"),
                      select_contexts_base("CrossDomain"))

metrics_to_plot <- c("Pairwise_Cosine","deg_pcc","jaccard_abs","log10_edist","N_Instances")

df_plot_long <- ctx_context %>%
  inner_join(sel_base, by = c("CellK","TargetK","DomainType")) %>%
  select(Cell, Target, DomainType, Kind, Track, View, all_of(metrics_to_plot)) %>%
  pivot_longer(cols = all_of(metrics_to_plot), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Track = factor(Track, levels = c("Gene","Pathway")),
    View  = factor(View, levels = c("Standard","Systema")),
    DomainType = factor(DomainType, levels = c("SameDomain","CrossDomain")),
    Kind = factor(Kind, levels = c("High","Low","Common"))
  )

df_matrix <- df_plot_long %>%
  mutate(ColKey = paste0(Metric, "|", as.character(View))) %>%
  select(Cell, Target, DomainType, Kind, Track, ColKey, Value) %>%
  pivot_wider(names_from = ColKey, values_from = Value) %>%
  arrange(Track, DomainType, Kind, Cell, Target)

views_ord <- c("Standard","Systema")
col_keys_ordered <- as.vector(outer(metrics_to_plot, views_ord, function(m,v) paste0(m,"|",v)))
col_keys_ordered <- col_keys_ordered[col_keys_ordered %in% names(df_matrix)]

mat_raw <- as.matrix(df_matrix[, col_keys_ordered, drop = FALSE])
rownames(mat_raw) <- paste0(df_matrix$Cell, " | ", df_matrix$Target)
mat_scaled <- apply(mat_raw, 2, scale01_clip)
mat_scaled <- as.matrix(mat_scaled)

row_split_vec <- factor(
  paste(df_matrix$Track, df_matrix$DomainType, sep = "\n"),
  levels = c("Gene\nSameDomain","Gene\nCrossDomain","Pathway\nSameDomain","Pathway\nCrossDomain")
)

view_vec <- sapply(strsplit(colnames(mat_scaled), "\\|"), `[`, 2)

col_ha <- HeatmapAnnotation(
  View = factor(view_vec, levels = c("Standard","Systema")),
  col = list(View = PAL$view),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  simple_anno_size = unit(3, "mm")
)

kind_col <- c(
  High = PAL$cls[["Robust_high"]],
  Low = PAL$cls[["Robust_low"]],
  Common = "grey80"
)
row_ha <- rowAnnotation(
  Kind = df_matrix$Kind,
  col = list(Kind = kind_col),
  show_annotation_name = FALSE,
  width = unit(3, "mm")
)

clean_metrics <- c(
  Pairwise_Cosine = "Pairwise Cosine",
  deg_pcc = "DEG PCC",
  jaccard_abs = "Jaccard abs",
  log10_edist = "log10 EDist",
  N_Instances  = "N"
)
col_labels_raw <- sapply(strsplit(colnames(mat_scaled), "\\|"), `[`, 1)
col_labels_view <- sapply(strsplit(colnames(mat_scaled), "\\|"), `[`, 2)
col_labels_vec <- paste0(clean_metrics[col_labels_raw], "\n", ifelse(col_labels_view=="Standard","Std","Sys"))


# ---- Baseline marker positions for cosine (null baseline; mapped into the same 0-1 scale) ----
null_cos <- read_csv_safe(FILE_TASK1_NULL_COS)

null_cos_lut <- tibble(Track = character(0), View = character(0), baseline = numeric(0))
if (nrow(null_cos) > 0 && all(c("Track","View") %in% names(null_cos))) {
  base_col <- if ("null_p50" %in% names(null_cos)) "null_p50" else if ("null_mu" %in% names(null_cos)) "null_mu" else NA_character_
  if (!is.na(base_col)) {
    null_cos_lut <- null_cos %>%
      transmute(
        Track = as.character(Track),
        View  = as.character(View),
        baseline = as.numeric(.data[[base_col]])
      ) %>%
      distinct()
  }
}

get_baseline <- function(track, view) {
  if (nrow(null_cos_lut) == 0) return(NA_real_)
  k <- paste0(track, "||", view)
  kk <- paste0(null_cos_lut$Track, "||", null_cos_lut$View)
  ii <- match(k, kk)
  if (is.na(ii)) return(NA_real_)
  null_cos_lut$baseline[ii]
}

scale01_point <- function(x, lo, hi) {
  if (!is.finite(x) || !is.finite(lo) || !is.finite(hi) || lo == hi) return(NA_real_)
  x2 <- pmin(pmax(x, lo), hi)
  (x2 - lo) / (hi - lo)
}

# Column-wise quantiles consistent with scale01_clip()
col_lo <- apply(mat_raw, 2, function(z) as.numeric(quantile(z, 0.05, na.rm = TRUE)))
col_hi <- apply(mat_raw, 2, function(z) as.numeric(quantile(z, 0.95, na.rm = TRUE)))

baseline_scaled <- matrix(NA_real_, nrow = nrow(mat_raw), ncol = ncol(mat_raw))
colnames(baseline_scaled) <- colnames(mat_raw)

cos_cols <- which(col_labels_raw == "Pairwise_Cosine")
if (length(cos_cols) > 0 && nrow(null_cos_lut) > 0) {
  for (j in cos_cols) {
    vj <- view_vec[j]
    lo <- col_lo[j]; hi <- col_hi[j]
    for (i in seq_len(nrow(mat_raw))) {
      tr <- as.character(df_matrix$Track[i])
      b0 <- get_baseline(tr, vj)
      baseline_scaled[i, j] <- scale01_point(b0, lo, hi)
    }
  }
}

col_std <- circlize::colorRamp2(c(0, 1), c("#E8F1FF", PAL$view[["Standard"]]))
col_sys <- circlize::colorRamp2(c(0, 1), c("#FFF3C4", PAL$view[["Systema"]]))

cell_bar_fun <- function(j, i, x, y, width, height, fill) {
  v <- mat_scaled[i, j]
  grid.rect(
    x, y,
    width = width * 0.90, height = height * 0.78,
    gp = gpar(fill = "#F5F5F5", col = NA)
  )
  if (is.na(v)) return()

  v_clamp <- max(v, 0.04)
  bar_w <- (width * 0.90) * v_clamp
  vname <- view_vec[j]
  bar_fill <- if (vname == "Standard") col_std(v) else col_sys(v)

  grid.rect(
    x - (width * 0.90) / 2 + bar_w / 2, y,
    width = bar_w, height = height * 0.78,
    gp = gpar(fill = bar_fill, col = NA, alpha = 1.0)
  )

  # Null baseline marker (cosine only; baseline_scaled is NA for other metrics)
  b <- baseline_scaled[i, j]
  if (is.finite(b)) {
    x0 <- x - (width * 0.90) / 2 + (width * 0.90) * b
    y0 <- y - (height * 0.78) / 2
    y1 <- y + (height * 0.78) / 2
    grid.lines(
      x = unit.c(x0, x0),
      y = unit.c(y0, y1),
      gp = gpar(col = "grey35", lwd = 0.6, lty = 2)
    )
  }
}

ht1c <- Heatmap(
  mat_scaled,
  name = "Scaled Score",
  cell_fun = cell_bar_fun,
  rect_gp = gpar(type = "none"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split_vec,
  bottom_annotation = col_ha,
  left_annotation = row_ha,
  row_names_gp = gpar(fontsize = 7),
  column_labels = col_labels_vec,
  column_names_gp = gpar(fontsize = 7),
  column_names_side = "top",
  column_names_rot = 45,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  column_title = "Fig 1C - Context Scoreboard (bar-based; Std and Sys)",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  border = FALSE
)

out1c <- file.path(OUT_DIR, "Fig1C_Scoreboard_BarBased_StdSys.pdf")
if (has_cairo) Cairo::cairo_pdf(out1c, width = 10.2, height = 10.8) else pdf(out1c, width = 10.2, height = 10.8)
draw(ht1c, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
message("Saved: ", out1c)

# ==============================================================================
# Fig 1D - Fingerprints (requires Step1_L1_Instance_Tidy.csv)
# ==============================================================================
message(">>> Fig 1D - Fingerprints (purple heatmap + overall bubble)...")

if (!file.exists(FILE_STEP1_INST)) {
  message("Skipping Fig 1D: Step1_L1_Instance_Tidy.csv not found.")
} else {

  inst <- readr::read_csv(FILE_STEP1_INST, show_col_types = FALSE)

  must_have_cols(inst, c("Track","Scenario","View","DomainType","Cell","Target","True_Rank","N_Gallery","LogDose","TimeBin","Success_Score","MRR"),
                 "Step1_L1_Instance_Tidy")

  VIEW_USE_FIG1D <- "Standard"


# ---- Retrieval null baselines (used to build Delta_Success and Delta_MRR) ----
null_retr <- read_csv_safe(FILE_TASK1_RETR_NULL)
base_success <- NA_real_
base_mrr <- NA_real_
if (nrow(null_retr) > 0 && all(c("Track","Scenario","View","Direction","LabelType","success_mu","mrr_mu") %in% names(null_retr))) {
  base_row <- null_retr %>%
    filter(
      Track == "Gene",
      Scenario == "A_L1Base",
      View == VIEW_USE_FIG1D,
      Direction == "Chem->Gene",
      LabelType == "Target"
    ) %>%
    slice_head(n = 1)
  if (nrow(base_row) == 1) {
    base_success <- as.numeric(base_row$success_mu[[1]])
    base_mrr <- as.numeric(base_row$mrr_mu[[1]])
  }
}
if (!is.finite(base_success)) base_success <- 0.5
if (!is.finite(base_mrr)) base_mrr <- 0.0

  inst1 <- inst %>%
    filter(Scenario == "A_L1Base", Track == "Gene", View == VIEW_USE_FIG1D) %>%
    mutate(
      Cell = as.character(Cell),
      Target = as.character(Target),
      CellK = std_key(Cell),
      TargetK = std_key(Target),
      True_Rank = as.numeric(True_Rank),
      N_Gallery = as.numeric(N_Gallery),
      Success_Score = as.numeric(Success_Score),
      MRR = as.numeric(MRR),
      LogDose = as.numeric(LogDose),
      TimeBin = as.character(TimeBin)
    ) %>%
    filter(is.finite(True_Rank), is.finite(N_Gallery), N_Gallery >= 10) %>%
    mutate(
      Top1  = (True_Rank <= 1),
      Top5  = (True_Rank <= 5),
      Top10 = (True_Rank <= 10)
    ) %>%
    group_by(DomainType, CellK, TargetK) %>%
    mutate(
      DoseBin = ifelse(is.finite(LogDose),
                       ifelse(LogDose >= median(LogDose, na.rm = TRUE), "High", "Low"),
                       NA_character_)
    ) %>%
    ungroup()

  proto_rates <- inst1 %>%
    filter(!is.na(TimeBin), !is.na(DoseBin)) %>%
    group_by(DomainType, CellK, TargetK, TimeBin, DoseBin) %>%
    summarise(Top10_rate = mean(Top10, na.rm = TRUE), n = n(), .groups = "drop")

  proto_ctx <- proto_rates %>%
    group_by(DomainType, CellK, TargetK) %>%
    summarise(
      Protocol_Range = ifelse(n() >= 2, max(Top10_rate, na.rm = TRUE) - min(Top10_rate, na.rm = TRUE), NA_real_),
      Protocol_SD    = ifelse(n() >= 2, sd(Top10_rate, na.rm = TRUE), NA_real_),
      .groups = "drop"
    )

  fp <- inst1 %>%
    group_by(DomainType, CellK, TargetK) %>%
    summarise(
      Cell = dplyr::first(Cell),
      Target = dplyr::first(Target),
      n_instances = n(),
      Mean_Success = mean(Success_Score, na.rm = TRUE),
      Mean_MRR     = mean(MRR, na.rm = TRUE),
      Top1_rate  = mean(Top1, na.rm = TRUE),
      Top5_rate  = mean(Top5, na.rm = TRUE),
      Top10_rate = mean(Top10, na.rm = TRUE),
      Top10_Early = mean(Top10[TimeBin == "Early"], na.rm = TRUE),
      Top10_Mid   = mean(Top10[TimeBin == "Mid"],   na.rm = TRUE),
      Top10_Late  = mean(Top10[TimeBin == "Late"],  na.rm = TRUE),
      Top10_LowDose  = mean(Top10[DoseBin == "Low"],  na.rm = TRUE),
      Top10_HighDose = mean(Top10[DoseBin == "High"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(proto_ctx, by = c("DomainType","CellK","TargetK")) %>%
    mutate(
      Dose_Sensitivity = Top10_HighDose - Top10_LowDose,
      Delta_Success = Mean_Success - base_success,
      Delta_MRR     = Mean_MRR - base_mrr,
      Time_Range = pmap_dbl(list(Top10_Early, Top10_Mid, Top10_Late), function(a,b,c){
        z <- c(a,b,c); z <- z[is.finite(z)]
        if (length(z) < 2) return(NA_real_)
        max(z) - min(z)
      }),
      ContextID = paste0(Cell, " | ", Target)
    ) %>%
    left_join(labels %>% select(CellK, TargetK, Class4, Mean_Success_Lab, Peak_Success_Lab),
              by = c("CellK","TargetK")) %>%
    mutate(Class4 = ifelse(is.na(Class4) | Class4 == "", "Intermediate", Class4))

  classes_keep <- c("Robust_high","Protocol_sensitive","Robust_low")

  fp2 <- fp %>%
    filter(n_instances >= 8, Class4 %in% classes_keep) %>%
    group_by(DomainType, Class4) %>%
    arrange(desc(Top10_rate), desc(Protocol_Range)) %>%
    slice_head(n = 10) %>%
    ungroup()

  if (nrow(fp2) == 0) stop("Fig 1D: empty fp2 after filtering. Check label join and thresholds.", call. = FALSE)

  fp2 <- fp2 %>%
    mutate(
      Overall_Score = dplyr::coalesce(0.65 * Mean_Success_Lab + 0.35 * Peak_Success_Lab, Top10_rate)
    )

  hm_cols <- c(
  "Delta_Success","Delta_MRR",
  "Top1_rate","Top5_rate","Top10_rate",
  "Top10_Early","Top10_Mid","Top10_Late",
  "Top10_LowDose","Top10_HighDose",
  "Dose_Sensitivity","Time_Range","Protocol_Range","Protocol_SD"
)
hm_cols <- hm_cols[hm_cols %in% names(fp2)]

  mat_raw <- as.matrix(fp2[, hm_cols, drop = FALSE])

scale_signed <- function(x, q = 0.95) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  m <- as.numeric(quantile(abs(x), q, na.rm = TRUE))
  if (!is.finite(m) || m == 0) return(rep(0.5, length(x)))
  x2 <- pmin(pmax(x, -m), m)
  (x2 / (2 * m)) + 0.5
}

mat_scaled <- mat_raw
for (cc in colnames(mat_scaled)) {
  if (cc %in% c("Delta_Success","Delta_MRR")) {
    mat_scaled[, cc] <- scale_signed(mat_scaled[, cc])
  } else {
    mat_scaled[, cc] <- scale01_clip(mat_scaled[, cc])
  }
}
mat_scaled <- as.matrix(mat_scaled)

  row_labs <- make.unique(paste0(fp2$Target, " @ ", fp2$Cell))
  rownames(mat_scaled) <- row_labs

  col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#FCFBFD", "#9E9AC8", "#3F007D"))

  class_col <- c(
    Robust_low = PAL$cls[["Robust_low"]],
    Protocol_sensitive = PAL$cls[["Protocol_sensitive"]],
    Robust_high = PAL$cls[["Robust_high"]]
  )

  row_ha <- rowAnnotation(
    Class = fp2$Class4,
    Overall = anno_points(
      fp2$Overall_Score,
      pch = 16,
      size = unit(1.4, "mm"),
      gp = gpar(col = class_col[fp2$Class4]),
      axis = TRUE
    ),
    N = anno_barplot(fp2$n_instances,
                     gp = gpar(fill = "grey70", col = NA),
                     border = FALSE,
                     bar_width = 0.9),
    col = list(Class = class_col),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 8),
    annotation_width = unit(c(3.5, 16, 12), "mm")
  )

  row_split_df <- data.frame(
    Domain = factor(fp2$DomainType, levels = c("SameDomain","CrossDomain")),
    Class  = factor(fp2$Class4, levels = c("Robust_low","Protocol_sensitive","Robust_high"))
  )

  ht1d <- Heatmap(
    mat_scaled,
    name = "Scaled",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = row_ha,
    row_split = row_split_df,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 9, fontface = "bold"),
    row_names_gp = gpar(fontsize = 6.2, col = "grey15"),
    column_names_gp = gpar(fontsize = 7.2, col = "grey15"),
    rect_gp = gpar(col = "grey92", lwd = 0.3),
    column_title = paste0("Fig 1D - Fingerprints (", VIEW_USE_FIG1D, ")"),
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Scaled",
      direction = "horizontal",
      legend_width = unit(40, "mm"),
      at = c(0, 0.5, 1),
      labels = c("low","mid","high")
    )
  )

  out1d <- file.path(OUT_DIR, paste0("Fig1D_Fingerprint_Purple_WithOverall_", VIEW_USE_FIG1D, ".pdf"))
  if (has_cairo) Cairo::cairo_pdf(out1d, width = 11.0, height = 7.2) else pdf(out1d, width = 11.0, height = 7.2)
  draw(ht1d, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()
  message("Saved: ", out1d)
}

# ==============================================================================
# Fig 1E - Scatter panels: columns=View (Standard/Systema), rows=Track (Gene/Pathway)
# ==============================================================================
message(">>> Fig 1E - Scatter panels (View x Track)...")

scat <- ctx_context %>%
  filter(is.finite(Pairwise_Cosine), is.finite(Mean_Success)) %>%
  mutate(
    View = factor(View, levels = c("Standard","Systema")),
    Track = factor(Track, levels = c("Gene","Pathway")),
    DomainType = factor(DomainType, levels = c("SameDomain","CrossDomain")),
    Class4 = factor(Class4, levels = c("Robust_low","Intermediate","Protocol_sensitive","Robust_high"))
  )

facet_stat <- scat %>%
  group_by(Track, View) %>%
  summarise(
    n = n(),
    r = suppressWarnings(cor(Pairwise_Cosine, Mean_Success, method = "spearman", use = "complete.obs")),
    p = suppressWarnings(cor.test(Pairwise_Cosine, Mean_Success, method = "spearman")$p.value),
    x = quantile(Pairwise_Cosine, 0.05, na.rm = TRUE),
    y = quantile(Mean_Success, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lab = paste0("Spearman r = ", sprintf("%.2f", r),
                 "\n", "p-value = ", format(p, digits = 2, scientific = TRUE),
                 "\n", "n = ", n)
  )

p1e <- ggplot(scat, aes(x = Pairwise_Cosine, y = Mean_Success)) +
  geom_point(aes(color = Class4, shape = DomainType), size = 1.15, alpha = 0.65) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = "grey35") +
  facet_grid(Track ~ View) +
  geom_text(
    data = facet_stat,
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 3.0, color = "grey15"
  ) +
  scale_color_manual(values = c(
    Robust_high = PAL$cls[["Robust_high"]],
    Protocol_sensitive = PAL$cls[["Protocol_sensitive"]],
    Intermediate = PAL$cls[["Intermediate"]],
    Robust_low = PAL$cls[["Robust_low"]]
  )) +
  scale_shape_manual(values = c(SameDomain = 16, CrossDomain = 17)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Fig 1E - Pairwise alignment vs retrieval success (panel by View and Track)",
    subtitle = "Points are contexts (Cell, Target) aggregated across source pairs; color=class; shape=domain.",
    x = "Pairwise cosine (chem vs gene)",
    y = "Mean success"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

save_gg(p1e, "Fig1E_Scatter_PanelByViewAndTrack.pdf", w = 10.2, h = 6.2, has_cairo = has_cairo)

message(">>> Done. Figure 1 outputs written to: ")
message(OUT_DIR)
