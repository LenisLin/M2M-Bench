# M2M-Bench Agent Guide

Start with `docs/redesign_checkpoint.md`. It defines the current benchmark
system and the shared terminology for the rest of the repo.

## Always-Active Repo Rules

### Source Hierarchy

When repo documents disagree, use this order:

1. Audited manifests and stage outputs
2. `docs/contracts/*.md`
3. `docs/governance/*.md`
4. `docs/manuscript_master.md` and `docs/plotting/*.md`

### Minimum Grounding Before Changes

Read the smallest relevant subset of:

- `docs/redesign_checkpoint.md`
- `docs/governance/state.md`
- `docs/governance/runbook.md`
- `docs/governance/repo_conventions.md`
- `docs/data_contracts.md` when touching task data, manifests, or result tables
- `docs/contracts/task1_spec.md` for Task1 semantics
- `docs/contracts/task2_spec.md` for Task2 semantics
- `docs/contracts/output-schemas.md` when editing outputs or validation logic
- `docs/manuscript_master.md` and plotting docs for manuscript or figure work

### Evidence And Storage Discipline

- Use NAS-backed roots for evidence discovery:
  - `/mnt/NAS_21T/ProjectData/M2M/runs`
  - `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
  - `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current`
- Treat the local checkout as source-only.
- Every non-trivial claim should cite a file path, manifest, table, or command
  result.
- If evidence is incomplete, say what still needs checking.

### Scientific Boundaries

- M2M-Bench is a benchmark/evaluation paper.
- Keep `Task1` and `Task2` separate.
- Keep `FM` scoped to the `scPerturb/K562` `Figure 3F` local-only panel.
- Update the relevant contract docs before changing benchmark semantics or
  figure meaning.

### Change Discipline

- Keep names, paths, and figure roles synchronized across docs.
- Remove inactive wording instead of layering alternate names on top of the
  current system.
- Do not fall back from corrected multisource Task2 outputs to a scPerturb-only
  path.

## Repo-Scoped Skills

Use repo-local skills when the task needs deeper workflow guidance:

- `.agents/skills/m2m-contract-grounding/SKILL.md`
- `.agents/skills/m2m-evidence-trace/SKILL.md`
- `.agents/skills/m2m-migration-runbook/SKILL.md`
