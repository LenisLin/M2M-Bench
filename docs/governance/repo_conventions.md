# Repository Conventions

This file collects repo-local engineering conventions. It is not a benchmark
contract and it is not a manuscript-spec replacement.

## Source hierarchy

When documents disagree, use this order:

1. Audited manifests and stage outputs
2. `docs/contracts/*.md`
3. `docs/governance/*.md`
4. `docs/manuscript_master.md` and plotting freeze docs
5. Support-only notes, retained historical files, and local staging outputs

## Active storage roles

- Canonical manuscript analysis:
  `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
- Active manuscript support:
  `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support`
- Historical manuscript archive:
  `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history`
- Retained reference namespaces:
  - `/mnt/NAS_21T/ProjectData/M2M/R_Vis_Ready`
  - `/mnt/NAS_21T/ProjectData/M2M/Model_Evaluation_Results`

`R_Vis_Ready` and `Model_Evaluation_Results` remain available for reference and
legacy comparison, but they are not the primary manuscript evidence root.

## Script header contract

Executable scripts under `scripts/` should keep a short file-level header that
states:

- current status (`active canonical`, `support-only`, or `historical retained`)
- manuscript or pipeline role
- active storage/output root
- any scope restrictions that are easy to violate

Keep those headers synchronized with actual behavior in the same patch.

## Evidence-first reporting

Non-trivial repo conclusions should cite concrete evidence such as:

- file paths and sections
- audited run manifests
- CSV or parquet outputs
- command results

If evidence is incomplete, say so directly and name the next verification step.

## README policy

`README.md` is maintained manually as the top-level navigation document.
`project.yaml` is retained as compact repository metadata, not as the source for
a generated README pipeline.
