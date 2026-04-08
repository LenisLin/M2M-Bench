# Exclusions & Policies (Frozen)

## Cross (Task1 S5)
- Cross chemical is EXCLUDED.
  - Reason: frozen policy/support gating does not admit a current-phase evaluable canonical cross-chemical slice.
  - Policy: do not compute or report any chemical cross metrics.
  - Evidence: manuscript-facing scope/support tables should report the exclusion through eligibility fields and frozen exclusion labels (for example `matched_keys_lt5`), not numeric matched-context rhetoric.

- Cross genetic only.
- Cross alignment MUST use a single contract:
  - LINCS side index: `global_idx_lincs`
  - scPerturb side index: `sc_delta_row_idx` (derived from `sc_m2m_row_idx -> delta_row_idx` via CRISPR meta)
  - No alternative alignment keys are allowed.

## LINCS representation policy
- LINCS is bulk: no FM representations.
- LINCS pathway space is derived by project_on_load:
  - pathway_delta = gene_delta @ hallmark_W_2477x50
  - W is loaded from a frozen artifact file (npy), not re-derived.

## Systema
- Systema view is removed in the reboot phase.
- Only Standard view is computed and reported.

"Systema view is strictly excluded from all analytical pipelines and outputs."