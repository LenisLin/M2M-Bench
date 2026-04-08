#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REVISION_ROOT="${M2M_REVISION_ROOT:-/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current}"
SOURCE_ANALYSIS_ROOT="${M2M_SOURCE_ANALYSIS_ROOT:-/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis}"
SOURCE_SUPPORT_ROOT="${M2M_SOURCE_SUPPORT_ROOT:-/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support}"
ANALYSIS_STAGED_ROOT="${M2M_ANALYSIS_STAGED_ROOT:-$REVISION_ROOT/analysis}"
PLOT_READY_ROOT="${M2M_PLOT_READY_ROOT:-$REVISION_ROOT/plot_ready}"
FIGURES_ROOT="${M2M_FIGURES_ROOT:-$REVISION_ROOT/figures}"
STAGED_MANIFEST_PATH="$ANALYSIS_STAGED_ROOT/framework_analysis_manifest.json"
COMPARISON_STATS_PATH="$ANALYSIS_STAGED_ROOT/manuscript_comparison_statistics.csv"
TASK1_BRIDGE_SUMMARY_PATH="${M2M_TASK1_BRIDGE_SUMMARY_PATH:-$SOURCE_SUPPORT_ROOT/group_bridge/task1_internal_vs_cross_group_bridge_summary.csv}"
TASK1_BRIDGE_DETAIL_PATH="${M2M_TASK1_BRIDGE_DETAIL_PATH:-$SOURCE_SUPPORT_ROOT/group_bridge/task1_internal_vs_cross_group_bridge.csv}"
SPATIAL_ENV_PREFIX="${M2M_SPATIAL_ENV_PREFIX:-/home/lenislin/miniconda3/envs/Spatial}"
SPATIAL_LD_LIBRARY_PATH="${M2M_SPATIAL_LD_LIBRARY_PATH:-$SPATIAL_ENV_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"

if [[ "$ANALYSIS_STAGED_ROOT" == "$SOURCE_ANALYSIS_ROOT" ]]; then
  echo "Analysis staging root must differ from the authoritative analysis root." >&2
  exit 1
fi
if [[ "$PLOT_READY_ROOT" == "$SOURCE_SUPPORT_ROOT/plot_ready" ]]; then
  echo "Plot-ready staging root must differ from the authoritative support root." >&2
  exit 1
fi
if [[ "$FIGURES_ROOT" == "$SOURCE_SUPPORT_ROOT/figures" ]]; then
  echo "Figure staging root must differ from the authoritative support root." >&2
  exit 1
fi

rm -rf "$ANALYSIS_STAGED_ROOT" "$PLOT_READY_ROOT" "$FIGURES_ROOT"
mkdir -p "$ANALYSIS_STAGED_ROOT"
mkdir -p "$PLOT_READY_ROOT"
mkdir -p "$FIGURES_ROOT"
rsync -a --no-o --no-g --delete "$SOURCE_ANALYSIS_ROOT/" "$ANALYSIS_STAGED_ROOT/"

python "$ROOT_DIR/scripts/manuscript_comparison_statistics.py" \
  --analysis-root "$SOURCE_ANALYSIS_ROOT" \
  --task1-bridge-summary-path "$TASK1_BRIDGE_SUMMARY_PATH" \
  --task1-bridge-detail-path "$TASK1_BRIDGE_DETAIL_PATH" \
  --output-path "$COMPARISON_STATS_PATH" \
  --manifest-path "$STAGED_MANIFEST_PATH"

conda run -n Spatial env LD_LIBRARY_PATH="$SPATIAL_LD_LIBRARY_PATH" \
  python "$ROOT_DIR/scripts/manuscript_plot_ready_tables.py" \
  --analysis-root "$SOURCE_ANALYSIS_ROOT" \
  --framework-manifest-path "$STAGED_MANIFEST_PATH" \
  --comparison-stats-path "$COMPARISON_STATS_PATH" \
  --plot-ready-root "$PLOT_READY_ROOT" \
  --build-target all

conda run -n Spatial env \
  M2M_ANALYSIS_ROOT="$ANALYSIS_STAGED_ROOT" \
  LD_LIBRARY_PATH="$SPATIAL_LD_LIBRARY_PATH" \
  M2M_PLOT_READY_ROOT="$PLOT_READY_ROOT" \
  M2M_FIGURES_ROOT="$FIGURES_ROOT" \
  Rscript "$ROOT_DIR/plotting/R/render_spatial_revision_panels.R"

if [[ ! -f "$FIGURES_ROOT/panel_export_manifest.tsv" ]]; then
  echo "Missing panel export manifest: $FIGURES_ROOT/panel_export_manifest.tsv" >&2
  exit 1
fi

unexpected_composed_output="$(
  find "$FIGURES_ROOT" -type f -name 'figure*.pdf' ! -name '*_panel_*' -print -quit
)"
if [[ -n "$unexpected_composed_output" ]]; then
  echo "Unexpected composed figure output: $unexpected_composed_output" >&2
  exit 1
fi

conda run -n Spatial env LD_LIBRARY_PATH="$SPATIAL_LD_LIBRARY_PATH" \
  python "$ROOT_DIR/scripts/manuscript_spatial_revision_reports.py" \
  --plot-ready-root "$PLOT_READY_ROOT" \
  --figures-root "$FIGURES_ROOT"

printf '%s\n' "$PLOT_READY_ROOT"
printf '%s\n' "$FIGURES_ROOT"
