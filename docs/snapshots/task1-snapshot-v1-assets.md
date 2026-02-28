# Data Assets to Keep (Task1)

## Must-keep (for rebuild)
### scPerturb deltas
- data/task1_snapshot_v1/scperturb_delta/scperturb-crispr-gene-delta.npy
- data/task1_snapshot_v1/scperturb_delta/scperturb-crispr-pathway-delta.npy
- data/task1_snapshot_v1/scperturb_delta/scperturb-crispr-delta-meta.csv
- data/task1_snapshot_v1/scperturb_delta/scperturb-drug-gene-delta.npy
- data/task1_snapshot_v1/scperturb_delta/scperturb-drug-pathway-delta.npy
- data/task1_snapshot_v1/scperturb_delta/scperturb-drug-delta-meta.csv

### FM deltas (aligned)
- data/task1_snapshot_v1/fm_delta/<fm_model>/{Genetic,Chemical}_delta_cell.npy
- data/task1_snapshot_v1/fm_delta/<fm_model>/{Genetic,Chemical}_delta_meta.csv

### LINCS bulk gene delta + meta (NAS)
- /home/lenislin/Experiment/projects/M2M/data/task1_snapshot_v1/lincs/lincs-engine1-gene-delta.pt
- /home/lenislin/Experiment/projects/M2M/data/task1_snapshot_v1/lincs/lincs-engine1-meta.csv

### Cross pairing contract + evidence
- outputs/m2m_v2/_review_extract/review_20260228_152654_task1_freeze_finalize/cross_pairs_genetic_contract.*
- outputs/m2m_v2/_review_extract/review_20260228_152654_task1_freeze_finalize/data/task1_snapshot_v1/cross_contract/cross-pairs-alignment-summary.json

### LINCS pathway projection artifacts
- outputs/m2m_v2/_review_extract/review_20260228_152654_task1_freeze_finalize/data/task1_snapshot_v1/pathway/hallmark-w-2477x50.npy
- outputs/m2m_v2/_review_extract/review_20260228_152654_task1_freeze_finalize/data/task1_snapshot_v1/pathway/lincs-pathway-policy.json
- /home/lenislin/Experiment/projects/M2M/data/task1_snapshot_v1/lincs/lincs-gene-alignment.csv

## Evidence packs to keep (recommended)
Keep the entire review folders for traceability:
- outputs/m2m_v2/_review_extract/review_20260228_114538_task1_freeze/
- outputs/m2m_v2/_review_extract/review_20260228_121425_lincs_locate_nas/
- outputs/m2m_v2/_review_extract/review_20260228_152654_task1_freeze_finalize/

## Task2 status
Task2 inputs are not frozen; do not rebuild Task2 until a dedicated Task2 snapshot exists.
