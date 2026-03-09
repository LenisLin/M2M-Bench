## Completed Stages

- **S2 complete and validated**: `scripts/s2_task1_cross_metrics.py`
  - Validation evidence: `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/audit_assertions.json`
  - Status: **9/9 assertions passed**
  - Scope validated: Task1 snapshot isolation, cross-contract validation, pathway projection-on-load policy, eligibility gate, deterministic parallelism, denominator conservation, chance identity tolerance, output routing, input path isolation.

- **S3 complete and validated**: `scripts/s3_build_task2_snapshot.py`
  - Validation evidence: `runs/s3_build_task2_data_0305/s3_build_task2_snapshot/audit_assertions.json`
  - Status: **7/7 assertions passed**
  - Scope validated: seed lock, raw K562 file presence, raw alignment, pairing contract, pathway projection contract, output routing, input path isolation.

- **Decoupled FM extractors complete and validated**
  - `scripts/fm_extractors/extract_scgpt.py`
    - Validation: `runs/fm_scgpt_0305/extract_scgpt/audit_assertions.json`
    - Status: **10/10 assertions passed**
  - `scripts/fm_extractors/extract_geneformer.py`
    - Validation: `runs/fm_geneformer_0305/extract_geneformer/audit_assertions.json`
    - Status: **11/11 assertions passed**
  - `scripts/fm_extractors/extract_scbert.py`
    - Validation: `runs/fm_scbert_0305/extract_scbert/audit_assertions.json`
    - Status: **11/11 assertions passed**
  - `scripts/fm_extractors/extract_scfoundation.py`
    - Validation: `runs/fm_scfoundation_0305/extract_scfoundation/audit_assertions.json`
    - Status: **11/11 assertions passed**
  - `scripts/fm_extractors/extract_uce.py`
    - Validation: `runs/fm_uce_0305/extract_uce/audit_assertions.json`
    - Status: **10/10 assertions passed**
  - `scripts/fm_extractors/extract_state.py`
    - Validation: `runs/0303/extract_state/audit_assertions.json`
    - Status: **10/10 assertions passed**
  - `scripts/fm_extractors/extract_tahoex1.py`
    - Validation: `runs/fm_tahoex1_0305/extract_tahoex1/audit_assertions.json`
    - Status: **11/11 assertions passed**

## Data Assets

### `data/task1_snapshot_v1/` verified state

- Frozen Task1 input snapshot is present and populated with:
  - `cross_contract/` assets
  - `pathway/hallmark-w-2477x50.npy`
  - `scperturb_delta/` gene + pathway deltas
  - `fm_delta/` bundles for `scgpt`, `geneformer`, `scbert`, `scfoundation`, `uce`, `state`, `tahoex1_3b`, plus `pca50`, `pca100`, `pca200`

- Verified representation dimensions:
  - Pathway projection matrix: **(2477, 50)**
  - `scperturb-drug`: gene **(282133, 2477)**, pathway **(282133, 50)**, meta rows **282133**
  - `scperturb-crispr`: gene **(1334725, 2477)**, pathway **(1334725, 50)**, meta rows **1334725**

- Verified FM delta bundle dimensions:
  - Chemical rows are **282133** for all FM/PCA bundles
  - Genetic rows are **1334725** for `scgpt`, `scbert`, `scfoundation`, `uce`, `state`, `tahoex1_3b`, `pca50`, `pca100`, `pca200`
  - `geneformer/Genetic_delta_cell.npy` is present at **(775317, 768)** with matching meta rows **775317**
  - FM/PCA representation widths:
    - `scgpt`: **512**
    - `geneformer`: **768**
    - `scbert`: **200**
    - `scfoundation`: **3072**
    - `uce`: **1280**
    - `state`: **2058**
    - `tahoex1_3b`: **2560**
    - `pca50` / `pca100` / `pca200`: **50 / 100 / 200**

### `data/task2_snapshot_v1/` verified state

- Frozen Task2 K562 snapshot is present and populated with:
  - Raw copied assets: `CRISPR_counts.pt`, `CRISPR_meta.csv`, `Drug_counts.pt`, `Drug_meta.csv`, `Common_Targets_K562.csv`, `shared_var_names.csv`
  - Derived core assets: `derived/pair_list.parquet`, `derived/delta_meta.csv`, `derived/gene_delta.npy`, `derived/pathway_delta.npy`
  - FM outputs under `k562/fm/` for `scgpt`, `geneformer`, `scbert`, `scfoundation`, `uce`, `state`, `tahoe-x1`

- Verified raw/base dimensions:
  - `CRISPR_meta.csv` / `CRISPR_counts.pt`: **14315 rows**
  - `Drug_meta.csv` / `Drug_counts.pt`: **30829 rows**
  - Shared gene dimension: **8363**

- **Strict Task2 row alignment is frozen at 30518 rows**
  - `derived/delta_meta.csv`: **30518 rows**
  - `derived/gene_delta.npy`: **(30518, 8363)**
  - `derived/pathway_delta.npy`: **(30518, 50)**
  - Every Task2 FM `fm_delta.npy` and `fm_delta_meta.csv` pair is aligned to the same **30518** rows

- Verified Task2 FM representation widths:
  - `scgpt`: **(30518, 512)**
  - `geneformer`: **(30518, 1152)**
  - `scbert`: **(30518, 200)**
  - `scfoundation`: **(30518, 3072)**
  - `uce`: **(30518, 1280)**
  - `state`: **(30518, 2058)**
  - `tahoe-x1`: **(30518, 512)**

## Core Contracts Enforced

- `GLOBAL_SEED=619` is rigidly enforced; scripts fail fast if the provided seed differs.
- Task2 row preservation is absolute:
  - `delta_meta.row_id` is contiguous `0..N-1`
  - FM extraction never drops rows
  - invalid rows are preserved as all-`np.nan` vectors with `valid_mask=False`
- `input_path_isolation` is enforced using `Path.absolute()` rather than `resolve()` so symlinked model roots preserve logical isolation semantics.
- S3 pairing contract is rigid:
  - `n_controls_used = min(pool_size, 50)`
  - controls are unique per treated row
  - selection is without replacement
  - missing control pools soft-fail into attrition with explicit reason tracking
- OOM handling is hardened in the extractor family via retry loops with dynamic batch-size reduction rather than silent failure.
- Tahoe-x1 extraction is wrapped through Docker subprocess execution with explicit image and shared-memory parameters.
- Output routing is isolated per stage under `runs/<run_id>/<stage>/`, with `run_manifest.json`, `audit_assertions.json`, and `manifest.json` emitted as standard AVCP artifacts.

## Pending Tasks

- **S4 remains to be implemented**: `scripts/s4_task2_group_concordance.py`
  - Current state: header contract exists; no completed implementation checkpoint.
  - Known risks requiring a fresh architectural approach:
    - dynamic globbing across FM assets
    - file descriptor exhaustion from too many simultaneously-open arrays/artifacts
    - chemical target cardinality complexity / imbalance
- **S5 remains to be implemented**: `scripts/s5_task2_retrieval.py`
  - Current state: header contract exists; no completed implementation checkpoint.

## Working Tree Note

- Current local working tree also contains unrelated modifications outside the requested checkpoint commits:
  - modified: `scripts/s0_build_data_inventory.py`
  - modified: `scripts/s1_task1_internal_metrics.py`
  - untracked: `scripts/tmp_audit_task2_fm_snapshot.py`
  - untracked: `runs/`
- The commit script below intentionally stages only the S2, S3, and FM extractor script files.