# M2M v2 Runbook

Canonical analysis intent reference: `docs/contracts/analysis-intent-lock.md`.
All executions and refactors must satisfy that contract first.
Tier gate protocol: `docs/m2m_v2_tier_gate.md`.
Output schema contract: `docs/m2m_v2_schema_contract.md`.
AVCP 8-script map: `docs/m2m_v2_avcp_8_scripts.md`.
Function refactor map: `docs/m2m_v2_function_refactor_map.md`.
Legacy inventory map: `docs/m2m_v2_legacy_inventory.md`.

Centralized runtime/path config: `config/config.yaml`.
All AVCP stage scripts and core analysis entry scripts load defaults from this file
(override with env `M2M_BENCH_CONFIG`).

This runbook defines the restarted execution chain aligned to the locked scope:

- `Task1`: modality consistency
  - delta-space internal (LINCS + scPerturb)
  - FM-delta internal (scPerturb)
- `Task2`: mechanism consistency
  - delta-space Chemical↔Genetic within source
  - FM-delta Chemical↔Genetic on K562

---

## AVCP 8-Script Entrypoints (preferred)

Inspect active config first:

```bash
python -c "import sys; from pathlib import Path; sys.path.insert(0, str(Path('src').resolve())); from m2m_bench.m2m_v2.core.config import config_path; print(config_path())"
```

```bash
python scripts/m2m_v2_s1_task1_data_curation.py --dry-run
python scripts/m2m_v2_s2_task2_data_curation.py --dry-run
python scripts/m2m_v2_s3_fm_embedding_generation.py --dry-run
python scripts/m2m_v2_s4_task1_internal_analysis.py --dry-run
python scripts/m2m_v2_s5_task1_cross_analysis.py --dry-run
python scripts/m2m_v2_s6_task1_representation_interpretation.py --dry-run
python scripts/m2m_v2_s7_task2_analysis.py --dry-run
python scripts/m2m_v2_s8_task2_representation_interpretation.py --dry-run
```

Then rerun without `--dry-run` per approved tier.

---

## Task1 Curation Contract (internal-first)

- Stage S1 now executes in this order:
  1) `task1_data_curation` (build `internal_candidates_full` first, then cross candidates and matched pairs),
  2) `build_scperturb_fm_dataset`,
  3) `data_curation_report`.
- For scPerturb, treated/control pairing is built from raw scRNA observations, while row alignment/scope follows processed metadata.
- scPerturb `gene_delta`, `pathway_delta`, and FM embedding delta must share the same control contract:
  - `control_context = dataset_std + cell_std + modality`.
- Contract audit files from S1 FM dataset build:
  - `outputs/m2m_v2/fm_scperturb_dataset/analysis/treated_control_contract.csv`
  - `outputs/m2m_v2/fm_scperturb_dataset/analysis/control_pool_members.csv`

---

## Tier Gate Order (must pass in sequence)

1. **Tier-0 (docs/spec lock only):**
   - update intent/constraints/runbook/schema docs only,
   - no analysis-logic code edits.
2. **Tier-1 (code refactor + sanity):**
   - `py_compile`,
   - parser `--help`,
   - dry-run / tiny-sanity outputs.
3. **Tier-2 (full execute):**
   - full cohorts / full retrieval runs,
   - question pack and manuscript numbers.

Tier transition criteria are defined in `docs/m2m_v2_tier_gate.md`.

---

## Framework Lock (Code-Level, Canonical Reference)

This section is the authoritative implementation contract for all subsequent code changes.
Any refactor or new script must satisfy this section first.

### A) Scientific scope and endpoint hierarchy

- **Task1 (Modality)**:
  - evaluate internal reliability and cross-domain comparability under fixed perturbation semantics.
  - analysis engines: `pairwise` + `retrieval`.
- **Task2 (Mechanism)**:
  - evaluate Chemical↔Genetic consistency within source, then FM extension on K562.
  - analysis engines: `pairwise` + `retrieval`.
- **Endpoint hierarchy (locked)**:
  - primary endpoint in main claims: `Standard` view.
  - `Systema` is pre-defined **secondary sensitivity analysis**.
  - final presentation selection happens after full computation; no metric dropping at compute stage.

### B) Data objects and analysis units

- `row`: one perturbation profile record.
- `context`:
  - Task1 default key: source-specific reliability key (`source_db`, `modality`, `cell_std`, `target_std`).
  - Task2 mechanism key: (`source_db`, `context_dataset`, `cell_std`, `target_match`).
- `query`: one retrieval query instance.
- `gallery`: retrieval candidate pool for one query.
- `pairwise_context_result`: one row of context-level pairwise summary.
- `instance_result`: one row of per-query retrieval output.

### C) Representation spaces

- **Delta spaces**:
  - `gene_delta`
  - `pathway_delta`
- **FM-delta space**:
  - `delta_z = z_treated - mean(z_control_pool)`
  - control pool must be explicit in metadata (`control_pool_id`).

### D) Metric contract (Task1/Task2 unified)

- **Pairwise family**:
  - centroid-level cosine (primary for Task2; also reported in Task1).
  - sampled instance-level cross cosine (supplementary robustness).
  - HVG-Spearman (Top50 genes per context; see section E).
  - E-distance (legacy-compatible definition, with fixed downsampling protocol).
- **Retrieval family**:
  - per-query `true_rank`, `MRR`, `TopK`.
  - balanced retrieval (`fixed gallery`, `fixed true/query`) with repeated subsampling.
  - random baseline under balanced setting:
    - `Top1_random = 1/N`
    - `MRR_random = H_N/N`

### E) HVG and distance definitions (locked)

- **HVG-Top50** (context-specific):
  - for each context and representation, compute:
    - `mean_chem`, `mean_gen`
    - `score_g = abs(mean_chem[g] - mean_gen[g])`
  - select top 50 genes by `score_g`.
  - compute Spearman on the selected index set.
- **E-distance**:
  - must follow legacy sign and scaling convention.
  - output both raw and direction-harmonized form (`NegEDist = -EDist`) where needed.
  - use fixed `max_n` and `n_boot`/`n_repeat` with seed logging.

### F) Systema policy (secondary analysis)

- `Systema` must be implemented as a view transformation on the same base vectors:
  - `x_sys = x_std - b`
  - background `b` fallback policy: `cell -> domain -> global`.
- `Systema` vs `Standard` comparison:
  - compute paired deltas by identical units (`track × context × direction`).
  - statistical testing includes paired non-parametric + paired t-test (sensitivity).
  - control multiplicity using BH-FDR across the predefined comparison family.

### G) Covariate analysis policy (Task1/Task2 harmonized framework)

- Common 4-layer structure:
  - **Layer 1**: univariate association scan.
  - **Layer 2**: multivariable deconfounding model.
  - **Layer 3**: enrichment analysis on context classes.
  - **Layer 4**: robustness (LOO / stratified re-run / threshold perturbation).
- Task-specific covariates:
  - Task1 emphasis: `Δdose`, `Δtime`, match quality features.
  - Task2 emphasis: dose/time protocol features, target confidence/polypharmacology tiers.
- Output requirement:
  - every layer must produce machine-readable tables in `analysis/`.
  - no manuscript-only numbers without a corresponding table artifact.

### H) Parity requirement versus legacy

- m2m_v2 must maintain an auditable parity bridge to legacy definitions:
  - same key construction,
  - same tie/rank logic for retrieval,
  - same target-expansion rule for multi-target chemicals,
  - same Systema background fallback policy,
  - same metric formulas/sign conventions.
- Every major module upgrade must produce a parity report:
  - counts parity (exact),
  - deterministic metric parity (strict tolerance),
  - stochastic metric parity (distribution/CI tolerance).

### I) Script/module responsibility map

- Dataset contract build:
  - `src/m2m_bench/m2m_v2/apps/build_scperturb_fm_dataset.py`
  - `src/m2m_bench/m2m_v2/apps/build_evalset_k562_contract.py`
- FM embedding execution bridge:
  - `src/m2m_bench/m2m_v2/fm/run_embeddings.py`
  - `src/m2m_bench/m2m_v2/fm/embed_*.py`
- FM-delta build:
  - `src/m2m_bench/m2m_v2/apps/build_fm_delta_from_embeddings.py`
- Core analyses:
  - `src/m2m_bench/m2m_v2/apps/task1_internal.py`
  - `src/m2m_bench/m2m_v2/apps/task1_fm_internal.py`
  - `src/m2m_bench/m2m_v2/apps/task2_mechanism_delta.py`
  - `src/m2m_bench/m2m_v2/apps/task2_fm_k562.py`
- Shared metric/util layers:
  - `src/m2m_bench/m2m_v2/core/metrics.py`
  - `src/m2m_bench/m2m_v2/core/common.py`

### J) Execution discipline

- Compute first, curate figures later.
- Do not remove metrics at run time for presentation convenience.
- Persist seed, config, manifest, and key summary for every run.
- If a module violates this framework lock, refactor before generating new manuscript numbers.

---

## 0) Delta-space baselines (no FM dependency)

Data curation denominator report (run first):

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.data_curation_report \
  --unified-meta-path ./outputs/task0_curated/metadata/unified_meta.csv \
  --m1-candidates-path ./outputs/task1/data/m1_candidates.csv \
  --m1-matched-pairs-path ./outputs/task1/data/m1_matched_pairs.csv \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --fm-evalset-k562-contract-dir ./outputs/m2m_v2/fm_evalset_k562_contract \
  --output-dir ./outputs/m2m_v2/data_curation_report
```

Intent deep-check (must pass before full runs):

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.intent_deep_check \
  --output-dir ./outputs/m2m_v2/intent_deep_check
```

Schema contract check (recommended before Tier-2):

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.data_readiness_audit \
  --output-dir ./outputs/m2m_v2/data_readiness_audit
```

Legacy reference audit (recommended before any legacy cleanup/deletion):

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.legacy_reference_audit \
  --output-dir ./outputs/m2m_v2/legacy_reference_audit
```

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.task1_internal \
  --cohort-mode internal_full \
  --task0-unified-meta-path ./outputs/task0_curated/metadata/unified_meta.csv \
  --m1-candidates-path ./outputs/task1/data/m1_candidates.csv \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --output-dir ./outputs/m2m_v2/task1_internal \
  --retrieval-mode instance_to_centroid \
  --retrieval-num-workers 8 \
  --balanced-gallery-size 256 \
  --balanced-true-per-query 1
```

Task1 retrieval mode quick switch:

```bash
# exact instance retrieval (slowest)
--retrieval-mode instance

# query instance -> centroid gallery (recommended)
--retrieval-mode instance_to_centroid
```

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.task2_mechanism_delta \
  --unified-meta-path ./outputs/task0_curated/metadata/unified_meta.csv \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --output-dir ./outputs/m2m_v2/task2_mechanism_delta \
  --min-rows-per-modality 8 \
  --lincs-context-level cell_target \
  --scperturb-context-level target_only \
  --retrieval-num-workers 8 \
  --balanced-gallery-size 256 \
  --balanced-true-per-query 1
```

Task2 scPerturb relaxed coverage run (target-only context, for full-source sensitivity):

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.task2_mechanism_delta \
  --unified-meta-path ./outputs/task0_curated/metadata/unified_meta.csv \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --output-dir ./outputs/m2m_v2/task2_mechanism_delta_sc_target_only \
  --min-rows-per-modality 8 \
  --lincs-context-level cell_target \
  --scperturb-context-level target_only \
  --retrieval-num-workers 8 \
  --balanced-gallery-size 256 \
  --balanced-true-per-query 1
```

Task1 cross-modality (matched LINCS↔scPerturb subset):

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.task1_cross_modality \
  --task0-unified-meta-path ./outputs/task0_curated/metadata/unified_meta.csv \
  --m1-matched-pairs-path ./outputs/task1/data/m1_matched_pairs.csv \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --output-dir ./outputs/m2m_v2/task1_cross_modality \
  --min-rows-per-source 1 \
  --retrieval-num-workers 8 \
  --balanced-gallery-size 256 \
  --balanced-true-per-query 1
```

---

## 1) Build scPerturb FM dataset (treated + control)

`task1_full` mode (all scPerturb rows in Task1 scope):

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.build_scperturb_fm_dataset \
  --mode task1_full \
  --m1-candidates-path ./outputs/task1/data/m1_candidates.csv \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --raw-sc-dir /mnt/NAS_21T/ProjectData/OSMOSIS/raw/scPerturb_Processed \
  --output-dir ./outputs/m2m_v2/fm_scperturb_dataset \
  --min-treated-per-context 1 \
  --min-controls-per-pool 8
```

K562 contract mode (FM mechanism extension subset only):

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.build_evalset_k562_contract \
  --input-dir /mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562 \
  --output-dir ./outputs/m2m_v2/fm_evalset_k562_contract \
  --file-materialization symlink
```

---

## 2) Run FM embedding scripts (manual, model-specific env)

Use existing Task1 FM embedding wrappers, but point `--eval-dir` to:

- `./outputs/m2m_v2/fm_scperturb_dataset` (task1_full)
- or `./outputs/m2m_v2/fm_evalset_k562_contract` (task2_k562)

Embed outputs should be written under:

- `./outputs/m2m_v2/task1_fm_embeddings/<model>/` (task1_full)
- `./outputs/m2m_v2/task2_fm_embeddings_k562/<model>/` (task2_k562)

---

## 3) Build FM deltas from embeddings

Task1 FM delta:

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.build_fm_delta_from_embeddings \
  --fm-dataset-dir ./outputs/m2m_v2/fm_scperturb_dataset \
  --embedding-root ./outputs/m2m_v2/task1_fm_embeddings \
  --output-dir ./outputs/m2m_v2/fm_delta \
  --min-treated-per-context 8
```

Task2 K562 FM delta:

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.build_fm_delta_from_embeddings \
  --fm-dataset-dir ./outputs/m2m_v2/fm_evalset_k562_contract \
  --embedding-root ./outputs/m2m_v2/task2_fm_embeddings_k562 \
  --output-dir ./outputs/m2m_v2/fm_delta_k562 \
  --min-treated-per-context 8
```

---

## 4) FM-level analyses

Task1 FM internal:

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.task1_fm_internal \
  --fm-delta-root ./outputs/m2m_v2/fm_delta \
  --output-dir ./outputs/m2m_v2/task1_fm_internal \
  --retrieval-num-workers 8 \
  --balanced-gallery-size 256 \
  --balanced-true-per-query 1
```

Task2 FM K562:

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.task2_fm_k562 \
  --fm-delta-root ./outputs/m2m_v2/fm_delta_k562 \
  --output-dir ./outputs/m2m_v2/task2_fm_k562 \
  --min-rows-per-modality 8 \
  --retrieval-num-workers 8 \
  --balanced-gallery-size 256 \
  --balanced-true-per-query 1
```

---

## 5) Data readiness audit (before analysis claims)

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.data_readiness_audit \
  --output-dir ./outputs/m2m_v2/data_readiness_audit
```

Key artifact:
- `outputs/m2m_v2/data_readiness_audit/analysis/m2m_data_readiness_matrix.csv`

---

## 6) Covariate package (Task1 / Task2)

Task1 covariates (cross-modality context):

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.task1_covariates \
  --pairwise-context-path ./outputs/m2m_v2/task1_cross_modality/analysis/task1_cross_pairwise_context.csv \
  --matched-pairs-path ./outputs/task1/data/m1_matched_pairs.csv \
  --output-dir ./outputs/m2m_v2/task1_covariates \
  --n-perm 200
```

Task2 covariates (mechanism context):

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.task2_covariates \
  --pairwise-context-path ./outputs/m2m_v2/task2_mechanism_delta/analysis/task2_mechanism_pairwise_context.csv \
  --rows-expanded-path ./outputs/m2m_v2/task2_mechanism_delta/analysis/task2_mechanism_rows_expanded_for_matching.csv \
  --output-dir ./outputs/m2m_v2/task2_covariates \
  --n-perm 200
```

---

## 7) Question-oriented evidence pack

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.question_pack \
  --task1-internal-dir ./outputs/m2m_v2/task1_internal \
  --task1-cross-dir ./outputs/m2m_v2/task1_cross_modality \
  --task2-delta-dir ./outputs/m2m_v2/task2_mechanism_delta \
  --task1-fm-dir ./outputs/m2m_v2/task1_fm_internal \
  --task2-fm-dir ./outputs/m2m_v2/task2_fm_k562 \
  --task1-covariate-dir ./outputs/m2m_v2/task1_covariates \
  --task2-covariate-dir ./outputs/m2m_v2/task2_covariates \
  --data-curation-dir ./outputs/m2m_v2/data_curation_report \
  --output-dir ./outputs/m2m_v2/question_pack
```

---

## 8) One-command orchestrator (optional)

```bash
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.pipeline --dry-run
```

Then enable only needed steps by toggles.

---

## 9) Minimal Tier-1 sanity command set

```bash
python -m py_compile src/m2m_bench/m2m_v2/core/metrics.py
python -m py_compile src/m2m_bench/m2m_v2/apps/task1_internal.py src/m2m_bench/m2m_v2/apps/task1_cross_modality.py src/m2m_bench/m2m_v2/apps/task2_mechanism_delta.py src/m2m_bench/m2m_v2/apps/task1_fm_internal.py src/m2m_bench/m2m_v2/apps/task2_fm_k562.py src/m2m_bench/m2m_v2/apps/task1_covariates.py src/m2m_bench/m2m_v2/apps/task2_covariates.py src/m2m_bench/m2m_v2/apps/question_pack.py
PYTHONPATH=src python -m m2m_bench.m2m_v2.apps.pipeline --dry-run
```
