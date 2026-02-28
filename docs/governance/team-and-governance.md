# Team & Governance (M2M-Bench)

## Roles
### Human Lead (You)
- Owns: scope, benchmark positioning ("report facts"), and final decisions on disputes.
- Owns: mandatory audit checkpoints and release sign-off.
- Approves: any change to task definitions, denominators, common support, exclusions, and reporting policy.
- Decides: what is "headline" vs "diagnostic-only".

### GPT Reviewer (ChatGPT; this assistant)
- Owns: mathematical/statistical audit and bounded gating checks.
- Produces: falsifiable audit rules, minimum evidence requirements, and contract-level risk diagnosis.
- Does NOT: run code; only audits based on CODEX evidence artifacts.

### Gemini (Co-author assistant)
- Owns: drafting CODEX prompts, drafting docs, and providing decision options to Human Lead.
- Supports: GPT Reviewer by turning audit requirements into executable extraction tasks.
- Does NOT: change scope unilaterally.

### CODEX (Engineer)
- Owns: implementation and execution.
- Must: follow AVCP contract (dry-run, run-isolation, manifests/hashes, fail-fast on contract violation).
- Must NOT: reinterpret tasks, modify definitions/policies, or invent conclusions.

## Dispute Resolution
1) CODEX emits evidence artifacts.
2) GPT Reviewer identifies violated/ambiguous contract clauses.
3) Gemini drafts revised prompt/doc patch options.
4) Human Lead decides.

## Mandatory Review Checkpoints
- A: Task1 snapshot freeze (input contract & alignment).
- B: Task1 metrics rebuild (denominators, LOO, chance correction).
- C: Any change to task definition/policy/exclusion.