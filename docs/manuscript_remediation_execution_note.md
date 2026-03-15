# Manuscript Remediation Execution Note

## Attempted Remediations

This pass attempted the four highest-priority manuscript-support remediations from the current local plan:

1. reconcile manuscript-facing wording for Task1 cross chemical exclusion
2. repair Task1 scPerturb internal common-scope indexing
3. create a manuscript panel-support ledger for retained main-text candidate slices
4. export two downstream manuscript-support sensitivity tables from existing Task2 per-query evidence

## Outcomes

### 1. Cross-chemical wording

- Status: succeeded for manuscript-facing planning files
- Current safe wording: cross chemical is excluded by frozen policy/support gating
- Implementation detail: manuscript-facing planning files no longer rely on a numeric matched-context rationale
- Active caution: the frozen contract policy file still contains numeric prose that is not supported by the current materialized S2 support tables, so that number should not be quoted in manuscript prose

### 2. Task1 scPerturb internal common-scope indexing

- Status: succeeded
- Local check result: scPerturb internal `Gene` and `Pathway` rows are present in both `task1_leaderboard_long.csv` and `task1_retrieval_summary.csv`
- Interpretation: the previous issue was under-extraction, not absence of local materialization
- Implementation detail: manuscript-facing indexing now distinguishes the local common-scope scPerturb slice from the additional FM-local scPerturb slice

### 3. Manuscript panel-support ledger

- Status: succeeded
- Output: `docs/manuscript_panel_support_ledger.csv`
- Result: 10 rows written, one for each currently retained main-text candidate slice
- Support check: all listed denominator fields were found in the inspected source files for those retained rows

### 4. Downstream sensitivity tables

- Status: succeeded
- Output: `docs/task2_c2g_query_n_targets_sensitivity.csv`
  - 779 rows written
  - scope: Task2 retrieval, `direction=C2G` only, dataset/cell_line explicit
- Output: `docs/task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv`
  - 27 rows written
  - scope: Task2 retrieval, `dataset=scPerturb`, `cell_line=K562`, `direction=C2G` only

## Files Updated Or Created

Updated:

- `docs/manuscript_source_register.md`
- `docs/manuscript_claim_to_evidence_map.csv`
- `docs/manuscript_blocker_memo.md`

Created:

- `docs/manuscript_panel_support_ledger.csv`
- `docs/task2_c2g_query_n_targets_sensitivity.csv`
- `docs/task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv`
- `docs/manuscript_remediation_execution_note.md`
- `scripts/manuscript_support_remediation.py`

## Unresolved Cautions That Remain Active

- `specificity_tier` remains unsupported as a benchmark-wide or group-level main-text axis; the new table is local downstream sensitivity evidence only.
- `n_targets` remains unsupported as a frozen core bidirectional/group-level axis; the new table is a C2G-only downstream sensitivity output.
- FM interpretation remains local to the supported scPerturb K562 scope and must not be generalized benchmark-wide.
- Cross-chemical exclusion should remain non-numeric in manuscript prose unless the numeric rationale is separately reconciled against frozen local evidence.
