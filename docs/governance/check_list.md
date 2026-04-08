# Preprocessing Review Checklist

Last updated: 2026-04-07

This checklist is the working review record for the preprocessing boundary of
M2M-Bench. It is a governance artifact, not a second contract source. When this
file conflicts with active contracts or audited manifests, the contract files
and manifests win.

## 1. Review Scope

This checklist covers the preprocessing boundary through `S3` only:

- `S0`: data inventory
- `S1`: Task1 internal metrics inputs and cohort construction
- `S2`: Task1 cross matched-input validation
- `S3`: corrected Task2 multisource snapshot build

This checklist does not review `S4-S7` benchmark-evaluation stages except where
`S3` must explicitly hand off a downstream interface.

Primary sources of truth:

- [runbook.md](/home/lenislin/Experiment/projects/M2M/docs/governance/runbook.md)
- [state.md](/home/lenislin/Experiment/projects/M2M/docs/governance/state.md)
- [task1_spec.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/task1_spec.md)
- [task2_spec.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md)
- [output-schemas.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/output-schemas.md)
- [ARCHITECTURE.md](/home/lenislin/Experiment/projects/M2M/scripts/ARCHITECTURE.md)

Status vocabulary:

- `PASS`: reviewed and aligned with current contract
- `FAIL`: reviewed and not aligned with current contract
- `BLOCKED`: cannot be declared aligned until another issue is fixed
- `NOT-IN-SCOPE`: outside the current preprocessing review boundary

## 2. Checklist

| Review item | Status | Evidence | Conclusion | Downstream interface | Notes |
| --- | --- | --- | --- | --- | --- |
| Authoritative preprocessing boundary is frozen to `S0-S3` | `PASS` | [runbook.md](/home/lenislin/Experiment/projects/M2M/docs/governance/runbook.md), [ARCHITECTURE.md](/home/lenislin/Experiment/projects/M2M/scripts/ARCHITECTURE.md) | Active preprocessing boundary is stable and does not include `S4-S7`. | `S0-S3` stage map only | Use this file for preprocessing review only. |
| Data and result roots are NAS-backed rather than repo-local | `PASS` | [state.md](/home/lenislin/Experiment/projects/M2M/docs/governance/state.md), [local_storage_policy.md](/home/lenislin/Experiment/projects/M2M/docs/governance/local_storage_policy.md) | Source-only local checkout and NAS-backed snapshots/runs are the active operating model. | All stages consume NAS-backed snapshots and emit NAS-backed runs | Do not interpret missing repo-local `data/` or `runs/` as contract failure. |
| `S0` inventory contract exists and is stage-bounded | `PASS` | [runbook.md](/home/lenislin/Experiment/projects/M2M/docs/governance/runbook.md), [output-schemas.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/output-schemas.md) | `S0` is a pure upstream inventory stage with explicit outputs. | `task1_data_inventory_long.csv`, `data_source_manifest.csv` | Inventory is upstream context, not a Task2 metric surface. |
| LINCS enters corrected Task2 as precomputed delta-space rather than recomputed pairing | `PASS` | [task2_spec.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md#L244), [s3_build_task2_multisource_snapshot.py](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py) | Corrected `S3` reuses Task1 LINCS delta-space assets and does not rebuild treated-control deltas for LINCS. | `gene_delta.npy`, `pathway_delta.npy`, `delta_meta.csv` under corrected Task2 snapshot | LINCS metadata is still retained for audit even when delta is reused. |
| scPerturb contribution to corrected Task2 remains explicit delta-space input with legacy provenance | `PASS` | [task2_spec.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md), [s3_build_task2_multisource_snapshot.py](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py), [s3_build_task2_snapshot.py](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_snapshot.py) | Corrected `S3` ingests the audited legacy scPerturb K562 delta-space subtree rather than silently rebuilding it. | `scperturb_k562/derived/*`, `scperturb_k562/fm/*` | Legacy is upstream evidence here, not live canonical stage logic. |
| Chemical target tokenization is source-aware | `PASS` | [task2_spec.md#L206](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md#L206), [s3_build_task2_multisource_snapshot.py#L221](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py#L221), [s3_build_task2_snapshot.py#L271](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_snapshot.py#L271), [test_task2_target_tokenization.py](/home/lenislin/Experiment/projects/M2M/tests/test_task2_target_tokenization.py) | LINCS chemical target parsing accepts `|`; scPerturb chemical parsing accepts raw `_` and canonical `;`. | `delta_meta.target_tokens`, `task2_row_membership.parquet` | Source-native delimiters are normalized before snapshot writing. |
| Corrected LINCS `delta_meta.target_tokens` preserves full instance-level target identity | `PASS` | [task2_spec.md#L256](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md#L256), [s3_build_task2_multisource_snapshot.py#L678](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py#L678), [test_task2_target_tokenization.py#L15](/home/lenislin/Experiment/projects/M2M/tests/test_task2_target_tokenization.py#L15) | Multi-target LINCS chemicals now keep the full token set in `delta_meta.target_tokens`. | `delta_meta.target_tokens` handed to `S4/S5` as instance identity | Eligibility filtering must stay in membership tables only. |
| Corrected membership tables restrict cohort eligibility without mutating row identity | `PASS` | [task2_spec.md#L228](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md#L228), [s3_build_task2_multisource_snapshot.py#L713](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py#L713) | `target_token` explode logic is downstream membership-only and no longer rewrites the stored snapshot identity. | `task2_row_membership.parquet`, `task2_pairs_coverage.csv` | This is the expected handoff for `S4/S5` cohort construction. |
| Corrected Task2 snapshot keeps required metadata fields for audit and downstream use | `PASS` | [task2_spec.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md), [output-schemas.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/output-schemas.md), [s3_build_task2_multisource_snapshot.py#L691](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py#L691) | Snapshot-level metadata includes `row_id`, `treated_cell_id`, `cell_line`, `target_raw`, `target_tokens`, `time`, `dose_value`, `n_controls_used`, and related audit fields. | `delta_meta.csv`, FM meta files, manifests | FM validity remains representation-specific and not universal. |
| Representation-scope policy is explicit at preprocessing boundary | `PASS` | [output-schemas.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/output-schemas.md), [task2_spec.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md), [ARCHITECTURE.md](/home/lenislin/Experiment/projects/M2M/scripts/ARCHITECTURE.md) | LINCS is limited to `Gene/Pathway`; scPerturb K562 keeps FM scope. Unsupported rows are scope policy, not attrition. | `representation_availability_registry.csv` and stage manifests | Downstream stages must preserve this distinction. |
| Pairing and manifest evidence are contractually present for preprocessing artifacts | `PASS` | [runbook.md](/home/lenislin/Experiment/projects/M2M/docs/governance/runbook.md), [output-schemas.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/output-schemas.md), [task2_spec.md#L398](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md#L398) | The preprocessing contract includes pairing evidence, manifests, and audit assertions rather than raw tensors alone. | `pair_list.parquet`, `snapshot_manifest.json`, `audit_assertions.json`, `manifest.json` | Review should always inspect manifests alongside data tables. |
| LINCS `cell_line` normalization is consistent across `S1`, `S2`, and corrected `S3` | `FAIL` | [s3_build_task2_multisource_snapshot.py#L116](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py#L116), [s3_build_task2_multisource_snapshot.py#L540](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py#L540), [s1_task1_internal_metrics.py#L701](/home/lenislin/Experiment/projects/M2M/scripts/s1_task1_internal_metrics.py#L701), [s2_task1_cross_metrics.py#L685](/home/lenislin/Experiment/projects/M2M/scripts/s2_task1_cross_metrics.py#L685) | Corrected `S3` normalizes `1HAE -> HA1E`, but current `S1/S2` still use raw LINCS `cell_line` text in cohort construction and cross checks. | `cell_line` enters Task1 and Task2 grouping keys | This must be harmonized before preprocessing can be called fully aligned. |
| Task1 and corrected Task2 preprocessing interfaces are fully harmonized | `BLOCKED` | Current checklist findings in this file | This cannot be marked aligned until the LINCS `cell_line` normalization mismatch is fixed or explicitly retired by contract. | `Task1 cohort keys`, `Task2 mech_key` handoff | Keep blocked even though the target-token contract is now fixed. |

## 3. Current Findings

Current preprocessing review outcome:

- `PASS`: delta-space routing, source-aware target parsing, full-target preservation, metadata retention, representation-scope policy, and manifest-based auditability are aligned with the current frozen contracts.
- `FAIL`: LINCS `cell_line` normalization is inconsistent between `S1/S2` and corrected `S3`.
- `BLOCKED`: full preprocessing-interface harmonization cannot be signed off until the `cell_line` alias issue is resolved.

Interpretation rules for this review:

- Do not mark preprocessing as globally complete while a `FAIL` item remains open.
- Do not open a second contract source here; contract changes belong in the contract documents, while this file records review state.
- When a failed item is repaired, update this file with the repair date, evidence path, and new status.

## 4. Six-Point Confirmation List

This section is the focused audit list for the six preprocessing guarantees that
must hold before downstream mechanism-concordance interpretation is trusted.

### 4.1 LINCS metadata is correctly built

- Confirm object:
  - LINCS metadata used by `Task1` and corrected `Task2`
- Core evidence:
  - [task1_spec.md](/home/lenislin/Experiment/projects/M2M/docs/contracts/task1_spec.md)
  - [s1_task1_internal_metrics.py#L697](/home/lenislin/Experiment/projects/M2M/scripts/s1_task1_internal_metrics.py#L697)
  - [s3_build_task2_multisource_snapshot.py#L517](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py#L517)
- Pass criteria:
  - `cell_line`, `target`, `pert_type`, `dose_val`, `time_val`, `treated_cell_id`, `paired_control_id`, and `sig_id` are present and normalized consistently
  - row identity is unique and stable
  - Task1 and Task2 do not apply conflicting metadata standardization rules
- Main risks:
  - LINCS `cell_line` alias drift across `S1/S2/S3`
  - target or perturbation type normalization mismatch

### 4.2 scPerturb pairing indices are correctly built

- Confirm object:
  - treated-control pairing evidence for scPerturb K562
- Core evidence:
  - [output-schemas.md#L363](/home/lenislin/Experiment/projects/M2M/docs/contracts/output-schemas.md#L363)
  - [s3_build_task2_snapshot.py#L295](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_snapshot.py#L295)
  - [s3_build_task2_snapshot.py#L374](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_snapshot.py#L374)
  - [s3_build_task2_snapshot.py#L474](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_snapshot.py#L474)
- Pass criteria:
  - each treated row resolves to a legal control pool
  - `pair_list.parquet` preserves `treated_cell_id`, `control_cell_id`, `control_rank`, `n_controls_used`, `dataset_side`, `time`, `dose_value`, and `seed`
  - CRISPR and Drug use the intended control-pool logic, and no-control cases are explicitly written to attrition
- Main risks:
  - cross-condition control mixing
  - unstable sampling or ranking
  - `n_controls_used` not matching actual pairs

### 4.3 scPerturb delta space is correctly built

- Confirm object:
  - `delta_meta.csv` + `gene_delta.npy` row universe for scPerturb K562
- Core evidence:
  - [task2_spec.md#L289](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md#L289)
  - [output-schemas.md#L369](/home/lenislin/Experiment/projects/M2M/docs/contracts/output-schemas.md#L369)
  - [s3_build_task2_snapshot.py#L645](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_snapshot.py#L645)
  - [s3_build_task2_snapshot.py#L680](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_snapshot.py#L680)
- Pass criteria:
  - `gene_delta.npy` row count equals `delta_meta.csv` row count
  - `row_id` is unique and contiguous
  - every delta row can be traced back to its pairing group
  - `target_raw`, `target_tokens`, `cell_line`, `perturbation_class`, and `n_controls_used` remain attached to the same row identity
- Main risks:
  - meta-array row drift
  - non-contiguous row ids
  - delta direction inconsistency

### 4.4 LINCS and scPerturb both enter pathway space correctly

- Confirm object:
  - pathway-space construction contract for both data sources
- Core evidence:
  - [exclusions-and-policies.md#L17](/home/lenislin/Experiment/projects/M2M/docs/contracts/exclusions-and-policies.md#L17)
  - [s2_task1_cross_metrics.py#L591](/home/lenislin/Experiment/projects/M2M/scripts/s2_task1_cross_metrics.py#L591)
  - [s3_build_task2_snapshot.py#L688](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_snapshot.py#L688)
  - [s3_build_task2_snapshot.py#L751](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_snapshot.py#L751)
- Pass criteria:
  - LINCS uses frozen `project_on_load` pathway projection with `hallmark_W_2477x50`
  - the pathway policy file matches the actual `W` artifact
  - scPerturb pathway rows stay perfectly aligned with scPerturb `delta_meta`
- Main risks:
  - wrong pathway policy mode
  - `W` version drift
  - pathway/meta row mismatch

### 4.5 scPerturb metadata is correctly built

- Confirm object:
  - scPerturb `delta_meta.csv` as the authoritative row-identity table
- Core evidence:
  - [task2_spec.md#L46](/home/lenislin/Experiment/projects/M2M/docs/contracts/task2_spec.md#L46)
  - [output-schemas.md#L369](/home/lenislin/Experiment/projects/M2M/docs/contracts/output-schemas.md#L369)
  - [s3_build_task2_multisource_snapshot.py#L321](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py#L321)
  - [s3_build_task2_multisource_snapshot.py#L425](/home/lenislin/Experiment/projects/M2M/scripts/s3_build_task2_multisource_snapshot.py#L425)
- Pass criteria:
  - required metadata columns are present: `row_id`, `treated_cell_id`, `perturbation_class`, `cell_line`, `target_raw`, `target_tokens`, `time`, `dose_value`, `specificity_tier`, `n_controls_used`
  - chemical `target_tokens` follow the source-aware parsing rule
  - membership explode logic does not create fake new row identities
- Main risks:
  - empty or malformed target tokens
  - metadata not matching pair/gene/pathway/FM layers
  - instance identity confused with target membership

### 4.6 scPerturb FM embeddings are correctly built

- Confirm object:
  - per-model FM deltas and FM metadata for scPerturb K562
- Core evidence:
  - [output-schemas.md#L373](/home/lenislin/Experiment/projects/M2M/docs/contracts/output-schemas.md#L373)
  - [extract_geneformer.py#L31](/home/lenislin/Experiment/projects/M2M/scripts/fm_extractors/extract_geneformer.py#L31)
  - [extract_geneformer.py#L1018](/home/lenislin/Experiment/projects/M2M/scripts/fm_extractors/extract_geneformer.py#L1018)
  - [extract_geneformer.py#L1338](/home/lenislin/Experiment/projects/M2M/scripts/fm_extractors/extract_geneformer.py#L1338)
- Pass criteria:
  - every FM model preserves the full `delta_meta` row universe with no row drops
  - `fm_delta.npy` row count equals `len(delta_meta)`
  - invalid rows are retained as `NaN` with `valid_mask=False`
  - `fm_delta_meta.csv` preserves `row_id`, `treated_cell_id`, `n_controls_used`, and `valid_mask`
- Main risks:
  - silent row dropping
  - invalid-mask and array-content mismatch
  - different FM models using different row universes

## 5. Stage Interfaces

### `S0 -> S1/S2`

- Interface object: inventory and source manifest
- Required outputs:
  - `task1_data_inventory_long.csv`
  - `data_source_manifest.csv`
- Required semantics:
  - enumerate available `dataset`, `perturbation_type`, `cell_line`, `target_token`
  - provide the audited inventory surface that justifies downstream Task1 coverage

### `Task1 snapshot -> S1`

- Interface object: Task1 internal input snapshot
- Required snapshot roots:
  - `/mnt/NAS_21T/ProjectData/M2M/snapshots/task1_snapshot_v1/`
- Required fields:
  - `cell_line`
  - `target`
  - `dose_val`
  - `time_val`
  - `pert_type`
  - `delta_row_index`
  - `delta_valid_bool`
- Required semantics:
  - build Task1 internal cohorts
  - load Gene/Pathway representations directly
  - keep FM only where allowed by Task1 scope

### `Task1 snapshot + cross contract -> S2`

- Interface object: Task1 cross matched subset
- Required fields:
  - `expected_cell_line`
  - `expected_target`
  - `global_idx_lincs`
  - `sc_delta_row_idx`
- Required semantics:
  - validate cross-source alignment
  - enforce cross eligibility gate
  - build Task1 cross retrieval and group inputs without leakage

### `Task1 snapshot + legacy scPerturb Task2 snapshot -> corrected S3`

- Interface object: corrected Task2 multisource build inputs
- Required input roots:
  - `/mnt/NAS_21T/ProjectData/M2M/snapshots/task1_snapshot_v1/`
  - `/mnt/NAS_21T/ProjectData/M2M/snapshots/task2_snapshot_v1/`
- Required output roots:
  - `/mnt/NAS_21T/ProjectData/M2M/snapshots/task2_snapshot_v2/`
- Required outputs:
  - `task2_pairs_coverage.csv`
  - `delta_meta.csv`
  - `task2_row_membership.parquet`
  - `representation_availability_registry.csv`
  - `snapshot_manifest.json`
- Required semantics:
  - reuse LINCS delta-space without recomputation
  - ingest audited scPerturb K562 delta-space subtree
  - normalize source-aware target tokens into canonical snapshot fields
  - preserve row identity separately from exploded membership

### `Corrected S3 -> S4`

- Interface object: Task2 group-concordance-ready snapshot
- Required downstream fields:
  - `dataset`
  - `cell_line`
  - `target_token`
  - `representation`
- Required semantics:
  - `S4` must operate on `dataset, cell_line, target_token` cohorts only
  - no raw pooling of LINCS and scPerturb

### `Corrected S3 -> S5`

- Interface object: Task2 retrieval-ready snapshot
- Required downstream fields:
  - `row_id`
  - `treated_cell_id`
  - `delta_meta.target_tokens`
  - `valid_mask`
  - `dataset`
  - `cell_line`
- Required semantics:
  - `S5` uses instance-level target identity from `delta_meta.target_tokens`
  - FM validity filtering is representation-specific
  - retrieval directions remain separate and in-scope only where the representation registry permits

## 6. Next Update Rules

- If `S1/S2` adopt the same LINCS `cell_line` alias map as corrected `S3`, update the failed checklist row to `PASS` and record the repair evidence.
- If contracts change, update the contract files first and then refresh this checklist.
- If preprocessing review expands beyond `S3`, add a new review section instead of silently widening this file's scope.
