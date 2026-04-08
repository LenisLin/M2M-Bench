# M2M-Bench Operating Prompt

Use the repository-local contracts and governance docs as the only source of truth for scope, storage, and manuscript interpretation.

## Required context before changes

Read the smallest relevant subset of:

- `docs/governance/state.md`
- `docs/governance/runbook.md`
- `docs/governance/repo_conventions.md`
- `docs/data_contracts.md`
- `docs/contracts/task1_spec.md`
- `docs/contracts/task2_spec.md`
- `docs/contracts/output-schemas.md`
- `docs/manuscript_master.md`
- `docs/plotting/plotting_preparation_freeze.md`

## Working rules

1. Treat benchmark contracts and audited manifests as higher priority than old notes or retained historical files.
2. Distinguish these roots explicitly:
   - canonical manuscript analysis: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
   - active manuscript support: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support`
   - historical manuscript surfaces: `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history`
   - retained reference namespaces: `/mnt/NAS_21T/ProjectData/M2M/R_Vis_Ready` and `/mnt/NAS_21T/ProjectData/M2M/Model_Evaluation_Results`
3. Do not promote support-only or historical files into canonical evidence by accident.
4. Keep `Task1` and corrected multisource `Task2` separate; do not reintroduce retired bridge language into current Figure 2 or Figure 3 claims.
5. When editing scripts, keep script headers, CLI defaults, and output paths consistent with the live storage layout.
6. When making conclusions, tie them to files, manifests, tables, or command output.

## Output discipline

- Keep changes scoped and explicit.
- Prefer evidence-backed summaries over template boilerplate.
- If a requested change would alter benchmark semantics, stop and update contracts before code.
