# Audit Protocol

Every audited stage passes three gates:

1. Input integrity: hashes, schemas, row alignment, and required files
2. Metric integrity: metric family, denominator fields, and aggregation rules
3. Leakage integrity: leave-one-out or disjoint-gallery behavior

Each stage bundle includes:

- `run_manifest.json`
- `audit_assertions.json`
- `manifest.json`
- the stage tables defined in `docs/contracts/output-schemas.md`

Each assertion must be falsifiable through concrete files and table rows.
