# M2M-Bench Project State (Onboarding SoT)

Last updated: 2026-02-19 (Tier-1 sanity + cleanup executed)

## Project objective
Deliver an audit-grade benchmark for:
- **Task1 (Modality):** internal reliability + cross-modality comparability
- **Task2 (Mechanism):** Chemical↔Genetic consistency within source

FM analyses are integrated as representation-space complements (not a separate primary task).

## Current analysis scope (locked)
- Task1 internal: LINCS + scPerturb full eligible cohorts
- Task1 cross: matched LINCS↔scPerturb subset only
- Task2 delta: full eligible mechanism cohorts
- Task2 FM: K562 contract subset only

## Primary active paths
- Intent lock: `docs/contracts/analysis-intent-lock.md`
- Execution runbook: `docs/governance/runbook.md`
- Core scripts:
  - `src/m2m_bench/m2m_v2/apps/task1_internal.py`
  - `src/m2m_bench/m2m_v2/apps/task1_cross_modality.py`
  - `src/m2m_bench/m2m_v2/apps/task2_mechanism_delta.py`
  - `src/m2m_bench/m2m_v2/apps/task1_fm_internal.py`
  - `src/m2m_bench/m2m_v2/apps/task2_fm_k562.py`
  - `src/m2m_bench/m2m_v2/apps/question_pack.py`

## Current engineering status
- Retrieval mainline unified to `instance_to_centroid`.
- `two_stage` removed from v2 mainline.
- `--retrieval-num-workers` wired into Task1/Task2 delta/FM scripts.
- Progress logging and memory cleanup hooks added in heavy pipelines.
- AVCP stage defaults centralized in `config/config.yaml` via `src/m2m_bench/m2m_v2/core/config.py`.
- Core analysis entrypoints (`task1_*`, `task2_*`, covariates, question_pack) now also read centralized defaults from `config/config.yaml`.
- AVCP tier-gate documents added and locked:
  - `docs/m2m_v2_tier_gate.md`
  - `docs/m2m_v2_schema_contract.md`
- AVCP 8-script orchestration layer added:
  - `docs/m2m_v2_avcp_8_scripts.md`
  - `scripts/m2m_v2_s1_*` ... `scripts/m2m_v2_s8_*`
- AVCP script-header contract coverage completed for active `scripts/*.py`.
- Legacy code/result cleanup executed inside `M2MBench/`:
  - legacy `scripts/task*` entrypoints removed from active tree
  - legacy result directories removed from `outputs/` except protected root-owned remnants
- Legacy inventory and audit assets:
  - `docs/m2m_v2_legacy_inventory.md`
  - `src/m2m_bench/m2m_v2/apps/legacy_reference_audit.py`

## Current stage
- Data Curation / Implementation
- 
## Next 3 tasks (approved queue)
1. Task1/Task2 script-flow normalization (no schema break).
2. Task1/Task2 covariate depth upgrade (enrichment + deconfounding).
3. Tier-1 sanity + dry-run pass, then Tier-2 full execution.

## Non-goals (do not drift)
- No new biological claim layers before finishing current locked pipeline.
- No silent metric redefinition outside documented SoT files.
- No full-run numbers reported without corresponding output tables/manifests.
