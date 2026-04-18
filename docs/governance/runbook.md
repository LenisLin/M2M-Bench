# Runbook

## Active Stage Map

| Stage | Purpose | Primary outputs |
| --- | --- | --- |
| `S0` | data inventory | `task1_data_inventory_long.csv`, `data_source_manifest.csv` |
| `S1` | Task1 internal metrics | `task1_group_concordance_long.csv`, `task1_retrieval_per_query.parquet`, `task1_retrieval_summary.csv`, `task1_leaderboard_long.csv` |
| `S2` | Task1 cross metrics | `task1_group_concordance_long.csv`, `task1_retrieval_per_query.parquet`, `task1_retrieval_summary.csv`, `task1_cross_alignment_proof.csv` |
| `S3` | Task2 multisource cohort build | `task2_pairs_coverage.csv`, `task2_row_membership.parquet`, `delta_meta.csv`, `representation_availability_registry.csv`, `snapshot_manifest.json` |
| `S4` | Task2 group concordance | `task2_group_concordance_long.csv`, `task2_group_leaderboard.csv` |
| `S5` | Task2 retrieval | `task2_retrieval_per_query.parquet`, `task2_retrieval_summary_long.csv`, `task2_retrieval_leaderboard.csv` |
| `S6` | Task2 synthesis | `task2_benchmark_summary_long.csv` |
| `S7` | project synthesis | `project_input_registry.csv`, `project_benchmark_summary_long.csv`, `project_axis_score_inputs_long.csv`, `project_representation_scorecard.csv` |

## Active Roots

- Task1 data: `/mnt/NAS_21T/ProjectData/M2M/data/task1`
- Task2 data: `/mnt/NAS_21T/ProjectData/M2M/data/task2`
- Stage runs: `/mnt/NAS_21T/ProjectData/M2M/runs`
- Manuscript analysis: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
- Plot review export: `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current`

## Stage Bundle Rule

Every stage writes:

- `run_manifest.json`
- `audit_assertions.json`
- `manifest.json`
- the stage tables defined in `docs/contracts/output-schemas.md`

## Stage Execution Rules

- One stage writes to one run directory.
- Task2 core metrics stay within each `dataset` and `cell_line`.
- `C2G` and `G2C` stay separate in Task2 retrieval outputs.
- The local checkout is source-only.
