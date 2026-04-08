#!/usr/bin/env python3
"""
Write semantic checks for the staged spatial figure revision pass.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

try:
    from manuscript_plot_ready_tables import (
        DEFAULT_FIGURE2_TARGET_HIGH_CONCORDANCE_PATH,
    )
except ModuleNotFoundError:
    from scripts.manuscript_plot_ready_tables import (
        DEFAULT_FIGURE2_TARGET_HIGH_CONCORDANCE_PATH,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REVISION_ROOT = Path("/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current")
DEFAULT_PLOT_READY_ROOT = DEFAULT_REVISION_ROOT / "plot_ready"
DEFAULT_FIGURES_ROOT = DEFAULT_REVISION_ROOT / "figures"
EXPECTED_PANEL_IDS = [
    "1B_LINCS",
    "1B_SCPERTURB",
    "2A",
    "2B",
    "2C",
    "2D",
    "2E",
    "2F",
    "3B",
    "3C",
    "3D",
    "3E",
    "3F",
    "EF2C",
    "EF4A",
    "EF4B",
    "EF5A",
    "EF5B",
    "EF5C",
    "EF6A",
    "EF6B",
    "EF7A",
    "EF7B",
    "EF7C",
    "EF7D",
    "EF8A",
    "EF8B",
    "EF8C",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write semantic checks for the staged spatial revision outputs.")
    parser.add_argument("--plot-ready-root", type=Path, default=DEFAULT_PLOT_READY_ROOT)
    parser.add_argument("--figures-root", type=Path, default=DEFAULT_FIGURES_ROOT)
    parser.add_argument(
        "--target-high-concordance-path",
        type=Path,
        default=DEFAULT_FIGURE2_TARGET_HIGH_CONCORDANCE_PATH,
    )
    parser.add_argument(
        "--markdown-output",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--tsv-output",
        type=Path,
        default=None,
    )
    args = parser.parse_args()
    if args.markdown_output is None:
        args.markdown_output = args.figures_root / "spatial_revision_semantics.md"
    if args.tsv_output is None:
        args.tsv_output = args.figures_root / "spatial_revision_semantics.tsv"
    return args


def read_csv_required(path: Path, label: str, required_columns: set[str] | None = None) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    read_kwargs = {"sep": "\t"} if path.suffix.lower() == ".tsv" else {}
    frame = pd.read_csv(path, **read_kwargs)
    if required_columns is not None:
        missing = sorted(required_columns - set(frame.columns))
        if missing:
            raise ValueError(f"{label} is missing required columns: {missing}")
    return frame


def format_range(min_value: int, max_value: int) -> str:
    if min_value == max_value:
        return str(min_value)
    return f"{min_value}-{max_value}"


def summarize_2a(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure2" / "figure2_panel_2a_task1_scope.csv"
    frame = read_csv_required(path, "Figure 2 panel 2A staged table")
    denom = (
        frame[["slice_label", "coverage_denominator"]]
        .drop_duplicates()
        .sort_values(["coverage_denominator", "slice_label"], kind="mergesort")
    )
    denom_text = "; ".join(
        f"{row.slice_label} = {int(row.coverage_denominator)}"
        for row in denom.itertuples(index=False)
    )
    count_types = (
        frame[["analysis_family", "count_annotation_type"]]
        .drop_duplicates()
        .sort_values(["analysis_family", "count_annotation_type"], kind="mergesort")
    )
    count_text = "; ".join(
        f"{row.analysis_family} -> {row.count_annotation_type}"
        for row in count_types.itertuples(index=False)
    )
    return {
        "panel_id": "2A",
        "topic": "Denominator",
        "detail": (
            "coverage_denominator is the lawful representation-surface count per Task1 slice: "
            f"{denom_text}. Count annotations use {count_text}."
        ),
    }


def summarize_1b(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure1" / "figure1_panel_1b_dataset_context_coverage.csv"
    frame = read_csv_required(
        path,
        "Figure 1 panel 1B staged table",
        {"task", "perturbation_type", "context_label", "context_folded_bool", "target_label_rule", "cell_line_label_rule"},
    )
    task_levels = ", ".join(sorted(frame["task"].dropna().astype(str).unique().tolist()))
    perturbation_levels = ", ".join(sorted(frame["perturbation_type"].dropna().astype(str).unique().tolist()))
    folded_counts = (
        frame.loc[frame["context_folded_bool"].astype(bool)]
        .groupby("dataset", dropna=False)["context_label"]
        .nunique()
        .reset_index(name="n_folded")
    )
    folded_text = "; ".join(
        f"{row.dataset} folded contexts = {int(row.n_folded)}"
        for row in folded_counts.itertuples(index=False)
    ) if not folded_counts.empty else "no folded contexts"
    return {
        "panel_id": "1B",
        "topic": "Sunburst Contract",
        "detail": (
            "The staged 1B source still preserves task provenance, but the live rendered panel collapses task identity "
            "and visualizes a compact sunburst with fixed hierarchy Perturbation type -> Cell line -> Target gene. "
            f"Observed task levels: {task_levels}. Perturbation levels: {perturbation_levels}. "
            f"Observed folding status: {folded_text}. "
            f"Context-label rule: {frame['target_label_rule'].dropna().astype(str).iloc[0]} "
            f"On-arc label rule: {frame['cell_line_label_rule'].dropna().astype(str).iloc[0]} "
            "In the rendered revision, dataset-local top-10 cell lines keep explicit colors and the folded Others slice is de-emphasized."
        ),
    }


def summarize_2c(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure2" / "figure2_panel_2c_gene_vs_pathway_matched_units.csv"
    frame = read_csv_required(
        path,
        "Figure 2 panel 2C staged table",
        {"row_kind", "comparison_context", "analysis_family", "metric_name", "unit_id"},
    )
    raw = frame.loc[frame["row_kind"].eq("raw")].copy()
    counts = (
        raw.groupby(["comparison_context", "analysis_family", "metric_name"], dropna=False)["unit_id"]
        .nunique()
        .reset_index(name="n_units")
    )
    ranges = (
        counts.groupby("comparison_context", dropna=False)["n_units"]
        .agg(["min", "max"])
        .reset_index()
    )
    range_text = "; ".join(
        f"{row.comparison_context} = {format_range(int(row['min']), int(row['max']))} matched units per metric"
        for _, row in ranges.iterrows()
    )
    return {
        "panel_id": "2C",
        "topic": "Triplet Units",
        "detail": (
            "The live panel is a unified Gene-vs-Pathway triplet surface over shared dataset + cell_line + target units, "
            "with cross restricted to the genetic bridge contract rather than mixed chemical/genetic slices. "
            f"In the staged table that yields {range_text}. Summary rows expose q10/q25/q50/q75/q90 so the panel can render "
            "10th-90th whiskers, IQR crossbars, and a median point without raw scatter."
        ),
    }


def summarize_2b(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure2" / "figure2_panel_2b_internal_performance_overview.csv"
    frame = read_csv_required(
        path,
        "Figure 2 panel 2B staged table",
        {"dataset_or_direction", "representation", "value"},
    )
    datasets = ", ".join(sorted(frame["dataset_or_direction"].dropna().astype(str).unique().tolist()))
    representations = (
        frame.groupby("dataset_or_direction", dropna=False)["representation"]
        .nunique()
        .reset_index(name="n_representations")
    )
    rep_text = "; ".join(
        f"{row.dataset_or_direction} = {int(row.n_representations)} lawful rows"
        for row in representations.itertuples(index=False)
    )
    return {
        "panel_id": "2B",
        "topic": "Raw-Scale Scoreboard",
        "detail": (
            "The live panel uses raw metric values as the scoreboard bar length rather than percentile fill. "
            f"Observed datasets: {datasets}. Lawful representation inventory: {rep_text}."
        ),
    }


def summarize_2d(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure2" / "figure2_panel_2d_internal_to_cross_degradation.csv"
    frame = read_csv_required(
        path,
        "Figure 2 panel 2D staged table",
        {"degradation_scope_note", "n_test_units"},
    )
    note = frame["degradation_scope_note"].dropna().astype(str).unique().tolist()
    detail = note[0] if note else "Internal-to-cross delta is computed over shared Task1 bridge units."
    min_n = int(frame["n_test_units"].min()) if not frame.empty else 0
    max_n = int(frame["n_test_units"].max()) if not frame.empty else 0
    return {
        "panel_id": "2D",
        "topic": "Intended Meaning",
        "detail": (
            f"{detail} The staged table now carries family-separated Group vs Retrieval summaries, "
            "scope-wise median/quantile intervals, and paired Gene/Pathway test metadata. "
            f"Observed n_test_units range = {format_range(min_n, max_n)}. "
            "Any eDist-family row is rendered from the logged value transform rather than the raw scale."
        ),
    }


def summarize_2e(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure2" / "figure2_panel_2e_cell_line_high_concordance_summary.csv"
    frame = read_csv_required(
        path,
        "Figure 2 panel 2E staged table",
        {"col_block", "row_block", "selection_direction", "shared_flag"},
    )
    counts = (
        frame[["col_block", "row_block", "selection_direction", "cell_line"]]
        .drop_duplicates()
        .groupby(["col_block", "row_block", "selection_direction"], dropna=False)["cell_line"]
        .nunique()
        .reset_index(name="n_examples")
    )
    count_text = "; ".join(
        f"{row.col_block}/{row.row_block}/{row.selection_direction} = {int(row.n_examples)}"
        for row in counts.itertuples(index=False)
    )
    shared_counts = (
        frame[["col_block", "row_block", "cell_line", "shared_flag"]]
        .drop_duplicates()
        .groupby(["col_block", "row_block"], dropna=False)["shared_flag"]
        .sum()
        .reset_index(name="n_shared")
    )
    shared_text = "; ".join(
        f"{row.col_block}/{row.row_block} shared exemplars = {int(row.n_shared)}"
        for row in shared_counts.itertuples(index=False)
    )
    return {
        "panel_id": "2E",
        "topic": "Enrichment Interpretation",
        "detail": (
            "Each row is an internal-only cell-line enrichment exemplar on log2(observed success rate / background success rate), "
            "with Gene and Pathway shown on the same axis and support moved to a right-side text column. "
            f"Observed dual-tail selection counts: {count_text}. Shared LINCS/scPerturb highlight counts: {shared_text}."
        ),
    }


def summarize_2f(plot_ready_root: Path, target_source_path: Path) -> dict[str, str]:
    path = plot_ready_root / "figure2" / "figure2_panel_2f_target_high_concordance_summary.csv"
    frame = read_csv_required(
        path,
        "Figure 2 panel 2F staged table",
        {"col_block", "row_block", "selection_direction", "shared_flag"},
    )
    raw = read_csv_required(target_source_path, "Figure 2 target high-concordance source", {"target_token"})
    raw_tokens = raw["target_token"].fillna("").astype(str)
    raw_unique = raw_tokens[raw_tokens.ne("")].nunique()
    underscore_rows = int(raw_tokens.str.contains("_", regex=False).sum())
    pipe_rows = int(raw_tokens.str.contains("|", regex=False).sum())
    both_rows = int(
        raw_tokens.str.contains("_", regex=False).mul(raw_tokens.str.contains("|", regex=False)).sum()
    )
    exemplar_counts = (
        frame[["col_block", "row_block", "selection_direction", "target_token"]]
        .drop_duplicates()
        .groupby(["col_block", "row_block", "selection_direction"], dropna=False)["target_token"]
        .nunique()
        .reset_index(name="n_examples")
    )
    exemplar_text = "; ".join(
        f"{row.col_block}/{row.row_block}/{row.selection_direction} = {int(row.n_examples)}"
        for row in exemplar_counts.itertuples(index=False)
    )
    return {
        "panel_id": "2F",
        "topic": "Atomic Target Enrichment",
        "detail": (
            f"The analysis-layer target summary now arrives at plot-ready staging already atomic: {raw_unique} unique target labels, "
            f"{underscore_rows} `_`-delimited rows, {pipe_rows} `|`-delimited rows, and {both_rows} rows with both delimiters. "
            f"The paired shared-axis staged table now uses dual-tail exemplar selection ({exemplar_text}); any family aliases are display-only labels and support is rendered as text rather than point size."
        ),
    }


def summarize_3b(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure3" / "figure3_panel_3b_c2g_performance_overview.csv"
    frame = read_csv_required(
        path,
        "Figure 3 panel 3B staged table",
        {"dataset", "cell_line", "analysis_family", "metric_name"},
    )
    row_counts = (
        frame[["dataset", "cell_line"]]
        .drop_duplicates()
        .groupby("dataset", dropna=False)["cell_line"]
        .nunique()
        .reset_index(name="n_rows")
    )
    row_text = "; ".join(
        f"{row.dataset} = {int(row.n_rows)} dataset|cell-line rows"
        for row in row_counts.itertuples(index=False)
    )
    metric_counts = (
        frame.groupby(["dataset", "analysis_family"], dropna=False)["metric_name"]
        .nunique()
        .reset_index(name="n_metrics")
    )
    metric_text = "; ".join(
        f"{row.dataset}/{row.analysis_family} = {int(row.n_metrics)} metrics"
        for row in metric_counts.itertuples(index=False)
    )
    return {
        "panel_id": "3B",
        "topic": "Benchmark Grammar",
        "detail": (
            "The live Figure 3B contract is a compact scoreboard bar surface: rows are dataset | cell line, "
            "columns are family-separated metric units, and each cell is rendered as paired Gene/Pathway lanes rather than a lollipop or fill-only tile. "
            f"Observed row inventory: {row_text}. Metric coverage: {metric_text}."
        ),
    }


def summarize_3c(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure3" / "figure3_panel_3c_cell_line_pattern.csv"
    frame = read_csv_required(
        path,
        "Figure 3 panel 3C staged table",
        {"dataset", "cell_line", "row_rank"},
    )
    max_rank = int(frame["row_rank"].max()) if not frame.empty else 0
    return {
        "panel_id": "3C",
        "topic": "Cell-Line Ranking",
        "detail": (
            "The live panel is now a LINCS-only ranked suitability surface rather than a dual-tail exemplar panel. "
            f"Observed maximum displayed row rank = {max_rank}."
        ),
    }


def summarize_3d(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure3" / "figure3_panel_3d_target_pattern_summary.csv"
    frame = read_csv_required(
        path,
        "Figure 3 panel 3D staged table",
        {"dataset", "target", "row_rank", "shared_across_modalities_bool"},
    )
    counts = (
        frame[["dataset", "row_rank", "target"]]
        .drop_duplicates()
        .groupby("dataset", dropna=False)["target"]
        .nunique()
        .reset_index(name="n_targets")
    )
    count_text = "; ".join(
        f"{row.dataset} = {int(row.n_targets)} ranked targets"
        for row in counts.itertuples(index=False)
    )
    shared_counts = (
        frame[["dataset", "target", "shared_across_modalities_bool"]]
        .drop_duplicates()
        .groupby("dataset", dropna=False)["shared_across_modalities_bool"]
        .sum()
        .reset_index(name="n_shared")
    )
    shared_text = "; ".join(
        f"{row.dataset} shared targets = {int(row.n_shared)}"
        for row in shared_counts.itertuples(index=False)
    )
    return {
        "panel_id": "3D",
        "topic": "Target Ranking",
        "detail": (
            "The displayed target set is now a ranked suitability surface rather than a dual-tail exemplar selection. "
            f"Observed ranked target counts: {count_text}. Shared-target highlight counts: {shared_text}."
        ),
    }


def summarize_3e(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure3" / "figure3_panel_3e_gene_vs_pathway_paired.csv"
    frame = read_csv_required(
        path,
        "Figure 3 panel 3E staged table",
        {"row_kind", "dataset", "analysis_family", "metric_name", "unit_id", "unit_source"},
    )
    raw = frame.loc[frame["row_kind"].eq("raw")].copy()
    dataset_counts = (
        raw.groupby(["dataset", "analysis_family", "metric_name"], dropna=False)["unit_id"]
        .nunique()
        .reset_index(name="n_units")
    )
    range_frame = (
        dataset_counts.groupby("dataset", dropna=False)["n_units"]
        .agg(["min", "max"])
        .reset_index()
    )
    range_text = "; ".join(
        f"{row.dataset} = {format_range(int(row['min']), int(row['max']))} paired units per metric"
        for _, row in range_frame.iterrows()
    )
    hidden_sources = ", ".join(sorted(raw["unit_source"].dropna().astype(str).unique().tolist()))
    return {
        "panel_id": "3E",
        "topic": "Target-Paired Contract",
        "detail": (
            f"The paired Gene-vs-Pathway panel now uses a single visible target-anchored surface ({hidden_sources}) "
            f"so each significance label aligns with the actual plotted paired units. The raw staged table contributes {range_text}. "
            "Both LINCS and scPerturb therefore use matched (cell_line, target) units rather than cell-line-only aggregation."
        ),
    }


def summarize_3f(plot_ready_root: Path) -> dict[str, str]:
    path = plot_ready_root / "figure3" / "figure3_panel_3f_fm_local_tradeoff.csv"
    frame = read_csv_required(
        path,
        "Figure 3 panel 3F staged table",
        {"representation", "test_status_vs_gene", "gene_reference_metric_value", "dataset", "cell_line"},
    )
    tested_rows = int(frame["test_status_vs_gene"].fillna("").eq("tested").sum())
    local_rows = int(frame["representation"].nunique())
    scope_text = ", ".join(sorted({f"{row.dataset}/{row.cell_line}" for row in frame.itertuples(index=False)}))
    return {
        "panel_id": "3F",
        "topic": "Local Trade-Off Framing",
        "detail": (
            "The live panel is plotted on absolute local performance for the supported scPerturb/K562 surface, not as a delta-vs-Gene effect axis. "
            f"The staged table carries {local_rows} representation rows, {tested_rows} Gene-reference significance rows in metadata, per-metric Gene reference values for the dashed vertical guide, and an explicit staged scope of {scope_text}."
        ),
    }


def supplementary_audit_rows(figures_root: Path) -> list[dict[str, str]]:
    manifest_path = figures_root / "panel_export_manifest.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Panel export manifest does not exist: {manifest_path}")
    manifest = read_csv_required(
        manifest_path,
        "Panel export manifest",
        {"panel_id", "qc_status", "legend_output_path", "rerender_count"},
    )
    observed_panel_ids = set(manifest["panel_id"].astype(str).tolist())
    missing_panel_ids = [panel_id for panel_id in EXPECTED_PANEL_IDS if panel_id not in observed_panel_ids]
    if missing_panel_ids:
        raise ValueError(
            "Panel export manifest is incomplete; missing panels: "
            + ", ".join(missing_panel_ids)
        )
    if manifest["qc_status"].astype(str).eq("manual_manifest").any():
        bad_panels = manifest.loc[manifest["qc_status"].astype(str).eq("manual_manifest"), "panel_id"].astype(str).tolist()
        raise ValueError(
            "Panel export manifest contains manual_manifest fallback rows: "
            + ", ".join(bad_panels)
        )
    missing_legend_panels: list[str] = []
    for row in manifest.itertuples(index=False):
        legend_path = getattr(row, "legend_output_path", "")
        if pd.isna(legend_path) or not str(legend_path):
            continue
        if not Path(str(legend_path)).exists():
            missing_legend_panels.append(str(row.panel_id))
    if missing_legend_panels:
        raise FileNotFoundError(
            "Legend PDFs are missing for: " + ", ".join(sorted(missing_legend_panels))
        )

    status_map = {
        panel_id: "auto_manifest_exported"
        for panel_id in sorted(observed_panel_ids)
        if panel_id.startswith("EF")
    }
    rows: list[dict[str, str]] = []
    for panel_id, audit_status in status_map.items():
        row = manifest.loc[manifest["panel_id"].eq(panel_id)].iloc[0]
        legend_state = "legend_exported" if pd.notna(row.get("legend_output_path", pd.NA)) and str(row.get("legend_output_path", "")) else "legend_not_required"
        rows.append(
            {
                "panel_id": panel_id,
                "topic": "Supplementary Audit",
                "detail": (
                    f"audit_status={audit_status}; manifest_qc_status={row['qc_status']}; "
                    f"rerender_count={int(row['rerender_count']) if pd.notna(row['rerender_count']) else 0}; "
                    f"{legend_state}"
                ),
            }
        )
    return rows


def build_rows(plot_ready_root: Path, target_source_path: Path, figures_root: Path) -> list[dict[str, str]]:
    return [
        summarize_1b(plot_ready_root),
        summarize_2a(plot_ready_root),
        summarize_2b(plot_ready_root),
        summarize_2c(plot_ready_root),
        summarize_2d(plot_ready_root),
        summarize_2e(plot_ready_root),
        summarize_2f(plot_ready_root, target_source_path),
        summarize_3b(plot_ready_root),
        summarize_3c(plot_ready_root),
        summarize_3d(plot_ready_root),
        summarize_3e(plot_ready_root),
        summarize_3f(plot_ready_root),
        *supplementary_audit_rows(figures_root),
    ]


def write_outputs(
    rows: list[dict[str, str]],
    markdown_output: Path,
    tsv_output: Path,
    figures_root: Path,
) -> None:
    markdown_output.parent.mkdir(parents=True, exist_ok=True)
    tsv_output.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = figures_root / "panel_export_manifest.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Panel export manifest does not exist: {manifest_path}")

    frame = pd.DataFrame(rows, columns=["panel_id", "topic", "detail"])
    frame.to_csv(tsv_output, sep="\t", index=False)

    lines = [
        "# Spatial Revision Semantics",
        "",
        f"- Panel manifest: `{manifest_path}`",
        f"- Legend source: `{ROOT / 'docs' / 'plotting' / 'manuscript_figure_legends.md'}`",
        "",
    ]
    for row in rows:
        lines.extend(
            [
                f"## Panel {row['panel_id']}",
                f"**{row['topic']}**",
                row["detail"],
                "",
            ]
        )
    markdown_output.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    rows = build_rows(
        args.plot_ready_root.resolve(),
        args.target_high_concordance_path.resolve(),
        args.figures_root.resolve(),
    )
    write_outputs(
        rows=rows,
        markdown_output=args.markdown_output.resolve(),
        tsv_output=args.tsv_output.resolve(),
        figures_root=args.figures_root.resolve(),
    )
    print(args.markdown_output.resolve())
    print(args.tsv_output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
