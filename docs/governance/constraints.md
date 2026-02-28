# M2M-Bench Constraints (Execution + Change Control)

## 1) Data and environment boundaries
- Do not modify NAS raw datasets.
- Do not modify legacy/original code under `Chem2Gene/` unless explicitly requested.
- Keep destructive cleanup limited to `outputs/m2m_v2` unless explicitly requested.

## 2) Analysis contract constraints
- `docs/contracts/analysis-intent-lock.md` is authoritative.
- Task scopes are fixed:
  - Task1 internal: full internal cohorts
  - Task1 cross: matched subset only
  - Task2 delta: full eligible mechanism cohorts
  - Task2 FM: K562 contract only
- Task1 scPerturb contract is fixed:
  - treated/control pairing is built from raw scRNA observations,
  - row selection/alignment must follow processed metadata (including manual metadata corrections),
  - gene/pathway recomputed delta and FM delta must share the same `control_pool_id` contract (`dataset+cell+modality`).
- Retrieval mainline mode is `instance_to_centroid`.

## 3) Engineering quality constraints
- Long-running scripts must print progress status.
- Heavy operations should expose parallel/splitting knobs.
- Memory-heavy loops should release intermediates (`del` + `gc.collect()`).
- Output schema changes require runbook update in same PR.
- Fail-fast only: no silent fallback that changes cohort, metric, or retrieval mode.
- If fallback is unavoidable, script must:
  - log explicit warning,
  - emit fallback tag in key summary/manifest,
  - keep primary endpoint tables separated from fallback outputs.

## 3.1 AVCP operational constraints
- Tiered execution is mandatory (`docs/m2m_v2_tier_gate.md`):
  - Tier-0: docs/spec lock,
  - Tier-1: code sanity and small dry-run checks,
  - Tier-2: full execution.
- Code change is not allowed to skip Tier-1 sanity unless explicitly approved.
- Every accepted run must be traceable to:
  - `run_manifest_*.json`,
  - `*_key_summary.csv`,
  - schema contract in `docs/m2m_v2_schema_contract.md`.
- Path/runtime defaults must be sourced from `config/config.yaml` (or explicit CLI args), not hidden hardcoded constants in stage wrappers or core analysis entry scripts.
- Before any legacy-script deletion/move, run legacy reference audit and confirm no mainline references:
  - `PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.legacy_reference_audit --output-dir ./outputs/m2m_v2/legacy_reference_audit`

## 4) Validation policy
- Before full execution, run sanity checks:
  - parser/help checks
  - `py_compile`
  - dry-run path checks if available
- Full execution is user-confirmed and user-driven.

## 5) Required artifacts for accepted analysis run
- `run_manifest_*.json`
- `*_key_summary.csv`
- pairwise detail + summary
- retrieval per-query + summary
- balanced/random baseline tables (retrieval-enabled runs)
- covariate/enrichment/LOO tables (covariate runs)

## 6) Reporting discipline
- Every manuscript-facing number must map to an output table path.
- If metric definitions are changed, document:
  - rationale,
  - compatibility impact,
  - migration note.
- Cross-task comparisons must state:
  - cohort scope (`internal_full`, `matched_subset`, `k562_contract`),
  - representation space (`gene_delta`, `pathway_delta`, `fm_delta`),
  - endpoint family (`pairwise`, `retrieval`, `covariate`).

## scPerturb Preprocessing Contract
- **Input Object**: 传入的 `.h5ad` 或 `adata.X` 矩阵，必须在上游（跨模态基准测试之前）完成 Size Factor Normalization (如 CPM/TPM) 和对数转换 (`log1p`)。
- **Delta Calculation**: 基准测试引擎在计算 `delta = x_treated - mean(control)` 时，不再进行二次归一化或插补操作。此契约保证了 Delta 空间的合法性，排除了测序深度带来的混杂偏差。