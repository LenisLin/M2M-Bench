# M2M-Bench

M2M-Bench is an audit-oriented benchmark for perturbation-response concordance
across two locked tasks:

- `Task1`: modality concordance under frozen `internal` and `cross` scopes
- `Task2`: mechanism concordance under corrected multisource `dataset` and
  `cell_line` stratification

FM analyses remain local complements rather than a third primary task.

## Start here

- Project state: `docs/governance/state.md`
- Execution rules: `docs/governance/runbook.md`
- Repo conventions: `docs/governance/repo_conventions.md`
- Data-contract index: `docs/data_contracts.md`
- Task contracts: `docs/contracts/task1_spec.md`,
  `docs/contracts/task2_spec.md`, `docs/contracts/output-schemas.md`
- Manuscript logic: `docs/manuscript_master.md`
- Plotting handoff: `docs/plotting/plotting_preparation_freeze.md`,
  `docs/plotting/manuscript_figure_legends.md`

## Active pipeline entrypoints

```bash
python -m pip install -e ".[dev]"
python scripts/s0_build_data_inventory.py --run-id <run_id> --seed 619
python scripts/s1_task1_internal_metrics.py --run-id <run_id> --seed 619
python scripts/s2_task1_cross_metrics.py --run-id <run_id> --seed 619
python scripts/s3_build_task2_multisource_snapshot.py --run-id <run_id> --seed 619
python scripts/s4_task2_group_concordance_multisource.py --run-id <run_id> --seed 619
python scripts/s5_task2_retrieval_multisource.py --run-id <run_id> --seed 619
python scripts/s6_task2_result_synthesis_multisource.py --run-id <run_id> --seed 619
python scripts/s7_project_benchmark_synthesis.py --run-id <run_id> --seed 619
```

## Storage layout

- Active authoritative benchmark runs:
  `/mnt/NAS_21T/ProjectData/M2M/runs`
- Archived benchmark runs:
  `/mnt/NAS_21T/ProjectData/M2M/archive/runs`
- Canonical manuscript analysis root:
  `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
- Active manuscript support root:
  `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support`
- Historical manuscript archive:
  `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history`
- Retained reference/result namespaces:
  - `/mnt/NAS_21T/ProjectData/M2M/R_Vis_Ready`
  - `/mnt/NAS_21T/ProjectData/M2M/Model_Evaluation_Results`

Only `runs/manuscript_active/analysis` should be treated as the manuscript
canonical evidence root. `manuscript_support` is the active plotting/support
layer; `archive/manuscript_history` contains retained historical and retired
surfaces.

## Local checkout notes

- Local `runs/` remains a compatibility directory for staging and path-stable
  inspection.
- Root `reports` is a symlink to `/mnt/NAS_21T/ProjectData/M2M/archive/reports`.
- `project.yaml` is retained as compact repo metadata.

## Current repo role split

- `scripts/s0`–`scripts/s7`: active benchmark pipeline
- `scripts/manuscript_framework_analysis_objects.py` and
  `scripts/manuscript_comparison_statistics.py`: active manuscript builders
- `scripts/manuscript_plot_ready_tables.py` and R plotting scripts: active
  manuscript support layer
- legacy bridge/backfill/sensitivity helpers: historical retained only

## Retained legacy note

Older AVCP-template scaffolding has been removed from the active entry surface.
If an older note or path conflicts with the current docs above, the current docs
and audited manifests win.
