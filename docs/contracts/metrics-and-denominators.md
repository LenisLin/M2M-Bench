# Metrics & Denominators (Frozen)

## Similarity / Distance
- Cosine similarity, Euclidean distance (and any other metric_family must not be mixed in aggregation).

## Retrieval
- Query -> centroid gallery retrieval; report AMRR/MRR/TopK as configured.
- Multi-positive chance correction is required when m>1:
  - corrected_MRR = MRR - E[MRR(N,m)]  (exact expectation used; implementation must be documented).
- LOO policy:
  - intra-dataset: strict recompute (query removed from centroid members).
  - cross-dataset: disjoint_no_leakage (query not in gallery by contract).

## Denominator policies
- Headline policy (Task1 reporting):
  - N_min_release fixed (e.g., 50) with sensitivity grid recorded (10/20/50).
  - Any reported headline must reference the policy row explicitly.
