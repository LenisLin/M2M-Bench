---
name: m2m-contract-grounding
description: Ground benchmark semantics, claim boundaries, and task definitions before changing M2M code or docs.
---

# M2M Contract Grounding

Use this skill when the task touches benchmark semantics, task wording, result
meaning, or figure roles.

## Load Order

Read only the smallest relevant subset, in this order:

1. `docs/redesign_checkpoint.md`
2. `docs/governance/state.md`
3. `docs/governance/repo_conventions.md`
4. the relevant contract docs under `docs/contracts/`
5. `docs/data_contracts.md` when task data or manifests are involved
6. `docs/manuscript_master.md` and plotting docs for manuscript or figure work

## Core Workflow

1. State the exact benchmark question or figure role the task touches.
2. Confirm whether the task belongs to `Task1`, `Task2`, or the
   `scPerturb/K562` `FM` panel.
3. Identify the highest-priority document that defines the behavior today.
4. Lock the active field names, unit keys, and panel meanings before editing.
5. Update docs or code only after the semantic boundary is clear.

## Guardrails

- Keep `Task1` and `Task2` separate.
- Keep `anchor_gene`, `perturbation_gene`, `query_instance_id`, and
  `pair_mean_enrichment` consistent across docs.
- Keep `FM` scoped to the `Figure 3F` local-only panel.
- Update contract docs before changing benchmark semantics.
