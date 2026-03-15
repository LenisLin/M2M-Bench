# Decisions

## D001 — Repository Identity Is M2M-Bench
- Context:
  older public metadata drifted from the local benchmark contracts and retained AVCP-template wording.
- Decision:
  public-facing repository identity tracks `M2M-Bench`; AVCP files remain governance/process tooling, not the project title.
- Evidence:
  `docs/contracts/project-positioning.md`, `docs/governance/state.md`, `README.md`, `project.yaml`
- Review Trigger:
  any approved rename or project repositioning contract

## D002 — Corrected Task2 Successor Root Is Frozen At `data/task2_snapshot_v2/`
- Context:
  earlier docs left the corrected Task2 successor root as a placeholder pending approval.
- Decision:
  corrected multisource Task2 S3-S6 consume `data/task2_snapshot_v2/`; `data/task2_snapshot_v1/` remains legacy/interim evidence only.
- Evidence:
  `data/task2_snapshot_v2/snapshot_manifest.json`, `docs/governance/state.md`, `runs/0310_fix_1hae/s3_build_task2_multisource_snapshot/run_manifest.json`
- Review Trigger:
  any approved snapshot version bump or successor-root replacement

## D003 — Corrected scPerturb Task2 Subtree Reuses Audited Legacy K562 Evidence
- Context:
  older migration text described rebuilding corrected scPerturb K562 from raw and referred generically to a `k562/` subtree.
- Decision:
  the corrected successor snapshot materializes `data/task2_snapshot_v2/scperturb_k562/` by reusing the audited legacy K562 bundle while preserving pairing and FM-policy evidence.
- Evidence:
  `runs/0310_fix_1hae/s3_build_task2_multisource_snapshot/run_manifest.json`, `data/task2_snapshot_v2/snapshot_manifest.json`
- Review Trigger:
  any approved change from reuse-based materialization to a new rebuild contract

## D004 — S7 Is Implemented With Local Implcheck Evidence
- Context:
  metadata drift left some public files describing S7 as pending or pointing to retired script names.
- Decision:
  repository state records S7 as implemented with implcheck evidence at `runs/s7_project_benchmark_synthesis_implcheck_20260311_a/s7_project_benchmark_synthesis/`.
- Evidence:
  `docs/governance/state.md`, `runs/s7_project_benchmark_synthesis_implcheck_20260311_a/s7_project_benchmark_synthesis/run_manifest.json`, `runs/s7_project_benchmark_synthesis_implcheck_20260311_a/s7_project_benchmark_synthesis/audit_assertions.json`
- Review Trigger:
  any approved change to the S7 contract, outputs, or versioned release state
