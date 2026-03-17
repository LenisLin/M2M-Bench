# M2M-Bench

M2M-Bench is an audit-oriented benchmark for perturbation-response concordance across two locked tasks:

- `Task1` evaluates modality concordance with frozen `internal` and `cross` scopes.
- `Task2` evaluates mechanism concordance with corrected multisource core metrics stratified by `dataset` and `cell_line`.

Local contracts in `docs/contracts/` and governance documents in `docs/governance/` are the source of truth for scope. The older public AVCP-template metadata is superseded by the local files in this repository.

## Current Project State

- Project objective:
  deliver an audit-grade benchmark for Task1 modality concordance and Task2 mechanism concordance across LINCS and scPerturb, with FM analyses treated as representation-space complements rather than a third primary task.
- Locked Task1 scope:
  `docs/governance/state.md`, `docs/governance/constraints.md`, and `docs/contracts/task1_spec.md` lock Task1 to `data/task1_snapshot_v1/`, with internal analyses on frozen eligible LINCS and scPerturb cohorts and cross analyses restricted to the frozen LINCS↔scPerturb matched contract.
- Locked Task2 scope:
  `docs/contracts/task2_spec.md` locks corrected Task2 to `mech_key=(dataset, cell_line, target_token)`, requires multisource core metrics from both `LINCS` and `scPerturb`, keeps `data/task2_snapshot_v1/` as legacy/interim scPerturb-K562 evidence, and limits FM representations to the scPerturb K562 subset.
- Current active pipelines:
  `scripts/s0_build_data_inventory.py`, `scripts/s1_task1_internal_metrics.py`, `scripts/s2_task1_cross_metrics.py`, `scripts/s3_build_task2_multisource_snapshot.py`, `scripts/s4_task2_group_concordance_multisource.py`, `scripts/s5_task2_retrieval_multisource.py`, `scripts/s6_task2_result_synthesis_multisource.py`, and `scripts/s7_project_benchmark_synthesis.py`.
- Completed stages with local audit evidence:
  `S0` (`runs/s0_build_data_inventory_0303/s0_build_data_inventory/`, 7/7 assertions passed), `S1` (`runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/`, 7/8 passed with one non-blocking soft FM-policy note), `S2` (`runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/`, 9/9 passed), legacy/interim `Task2 S3` (`runs/s3_build_task2_data_0305/s3_build_task2_snapshot/`, 7/7 passed), Task2 K562 FM extractors (`runs/fm_*` plus `runs/0303/extract_state/`, all audited), corrected multisource `S3` (`runs/0310_fix_1hae/s3_build_task2_multisource_snapshot/`, 14/14 passed), corrected multisource `S4` (`runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/`, 14/14 passed), corrected multisource `S5` (`runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/`, 15/15 passed), and corrected multisource `S6` (`runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/`, 11/11 passed).
- Implemented project-level synthesis:
  `S7` is now implemented as `scripts/s7_project_benchmark_synthesis.py`, with authoritative implcheck evidence at `runs/s7_project_benchmark_synthesis_implcheck_20260311_a/s7_project_benchmark_synthesis/`.
- Storage state:
  authoritative active results are now NAS-backed under `/mnt/NAS_21T/ProjectData/M2M/runs`, historical results are archived under `/mnt/NAS_21T/ProjectData/M2M/archive/runs`, local `runs/` intentionally remains a real directory, and local project-path compatibility is preserved via run-id-level symlinks under `runs/` plus a root `reports` symlink. See `docs/governance/local_storage_policy.md`.
- Pending or not yet evidenced:
  no `runs/**/*_key_summary.csv` files were found, so current status summaries come from manifests, audit assertions, and stage tables rather than key-summary CSVs.
- Authoritative output artifacts:
  `runs/<run_id>/s0_build_data_inventory/task1_data_inventory_long.csv`,
  `runs/<run_id>/s0_build_data_inventory/data_source_manifest.csv`,
  `runs/<run_id>/s1_task1_internal_metrics/task1_retrieval_per_query.parquet`,
  `runs/<run_id>/s1_task1_internal_metrics/task1_retrieval_summary.csv`,
  `runs/<run_id>/s1_task1_internal_metrics/task1_chance_identity_check.csv`,
  `runs/<run_id>/s1_task1_internal_metrics/task1_leaderboard_long.csv`,
  `runs/<run_id>/s1_task1_internal_metrics/task1_attrition.csv`,
  `runs/<run_id>/s2_task1_cross_metrics/task1_group_cross.parquet`,
  `runs/<run_id>/s2_task1_cross_metrics/task1_cross_retrieval_summary.csv`,
  `runs/<run_id>/s2_task1_cross_metrics/task1_cross_alignment_proof.csv`,
  `data/task2_snapshot_v2/snapshot_manifest.json`,
  `data/task2_snapshot_v2/task2_pairs_coverage.csv`,
  `data/task2_snapshot_v2/representation_availability_registry.csv`,
  `runs/<run_id>/s4_task2_group_concordance_multisource/task2_group_concordance.csv`,
  `runs/<run_id>/s4_task2_group_concordance_multisource/task2_group_attrition.csv`,
  `runs/<run_id>/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet`,
  `runs/<run_id>/s5_task2_retrieval_multisource/task2_retrieval_summary.csv`,
  `runs/<run_id>/s5_task2_retrieval_multisource/task2_retrieval_summary_long.csv`,
  `runs/<run_id>/s5_task2_retrieval_multisource/task2_retrieval_attrition.csv`,
  `runs/<run_id>/s5_task2_retrieval_multisource/task2_chance_identity_check.csv`,
  `runs/<run_id>/s6_task2_result_synthesis_multisource/task2_group_concordance_long.csv`,
  `runs/<run_id>/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv`,
  `runs/<run_id>/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv`,
  `runs/<run_id>/s6_task2_result_synthesis_multisource/task2_benchmark_summary_long.csv`,
  `runs/<run_id>/s7_project_benchmark_synthesis/project_input_registry.csv`,
  `runs/<run_id>/s7_project_benchmark_synthesis/project_benchmark_summary_long.csv`,
  `runs/<run_id>/s7_project_benchmark_synthesis/project_axis_score_inputs_long.csv`,
  `runs/<run_id>/s7_project_benchmark_synthesis/project_representation_scorecard.csv`.

## Script Inventory

### Current pipeline scripts

- `scripts/s0_build_data_inventory.py`: Task1 data inventory build from `data/task1_snapshot_v1/`.
- `scripts/s1_task1_internal_metrics.py`: Task1 internal group-level concordance and retrieval.
- `scripts/s2_task1_cross_metrics.py`: Task1 frozen cross-contract metrics and eligibility-gated retrieval.
- `scripts/s3_build_task2_multisource_snapshot.py`: corrected multisource Task2 snapshot materialization under `data/task2_snapshot_v2/`.
- `scripts/s4_task2_group_concordance_multisource.py`: corrected Task2 group concordance.
- `scripts/s5_task2_retrieval_multisource.py`: corrected Task2 target-level retrieval with multi-positive chance correction.
- `scripts/s6_task2_result_synthesis_multisource.py`: corrected Task2 leaderboard and benchmark summary synthesis.
- `scripts/s7_project_benchmark_synthesis.py`: project-level benchmark synthesis above Task1 S0-S2 and corrected Task2 S3-S6.

### Preserved historical/supporting scripts

- `scripts/s3_build_task2_snapshot.py`
- `scripts/s4_task2_group_concordance.py`
- `scripts/s5_task2_retrieval.py`
- `scripts/s6_task2_result_synthesis.py`
- `scripts/fm_extractors/extract_scgpt.py`
- `scripts/fm_extractors/extract_geneformer.py`
- `scripts/fm_extractors/extract_scbert.py`
- `scripts/fm_extractors/extract_scfoundation.py`
- `scripts/fm_extractors/extract_uce.py`
- `scripts/fm_extractors/extract_state.py`
- `scripts/fm_extractors/extract_tahoex1.py`

### Package surface

The installable package is currently minimal: `src/m2mbench/metrics/`, `src/m2mbench/utils/`, and `src/m2mbench/utils/logging.py`. Most benchmark workflow logic currently lives in `scripts/`.

## Representative Entrypoints

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

## Known Public Repo Drift

- The previous public `README.md` and `project.yaml` described this repository as `avcp-template`; the local source of truth identifies it as `m2m-bench`.
- No local authoritative owner field or license file was found, so public metadata is currently marked as `unknown` rather than guessed.
- No `runs/**/*_key_summary.csv` files were found in the local tree. Current status summaries are therefore derived from `run_manifest.json`, `audit_assertions.json`, and the stage tables they reference.

## Storage Layout

- Active authoritative run storage:
  `/mnt/NAS_21T/ProjectData/M2M/runs`
- Archived historical run storage:
  `/mnt/NAS_21T/ProjectData/M2M/archive/runs`
- Archived reports storage:
  `/mnt/NAS_21T/ProjectData/M2M/archive/reports`
- Local compatibility behavior:
  `runs/` intentionally remains a local directory, `runs/<run_id>` entry points are symlinks spanning the NAS active and NAS archive namespaces, and `reports` is a root symlink to the NAS archive path. Root-level `runs` normalization is not the current strategy.

## Version

- Package version: `0.4.0` (`pyproject.toml`)
- Repository version note: `v0.4.0` (`VERSION.md`)

The previous `VERSION.md` reboot note said Task2 was not yet frozen. That is superseded by the current corrected Task2 contracts and audited multisource S3-S6 runs.
