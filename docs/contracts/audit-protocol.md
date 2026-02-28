# Audit Protocol (Bounded)

We use a bounded audit approach:
- Gate 1: Snapshot integrity (hashes, schemas, alignment).
- Gate 2: Metric integrity (no metric_family mixing; denominator policies bound).
- Gate 3: Leakage integrity (LOO strictness vs disjoint).

Every gate must have:
- manifest.json (sha256 + size),
- a fact-only report,
- sample-based evidence (small CSV/JSONL) sufficient for third-party review.
