# Section 1: Purpose

This document freezes the implementation specification for Phase 1 manuscript-support analyses before full coding. It is a script-spec and skeleton-review artifact only. It does not execute analyses, does not generate manuscript result tables, and does not change benchmark-stage contracts.

# Section 2: Global coding principles

- downstream manuscript-support only
- no upstream stage mutation by default
- explicit local-only labeling where relevant
- explicit gating logic
- partial-output-over-forced-inference policy
- denominator and key validation requirements
- no silent coercion
- no silent join loss
- no silent fillna behavior

Shared status vocabulary:

- `successful`: the intended comparison or audit row was constructed on lawful inputs and passed required validation checks.
- `partial`: only a lawful subset of the intended output could be produced, with limitations carried explicitly in notes or status fields.
- `unavailable`: the requested row or comparison could not be materialized from current valid inputs, but no rule violation or blocking error occurred.
- `blocked`: implementation must stop for that row or comparison because a required lawful key, schema, or review decision is missing.
- `no_overlap`: exact overlap required for the audited comparison was not found.
- `partial_overlap`: some exact overlap exists, but not enough to treat the slice as fully usable without further support review.
- `usable_overlap`: exact overlap exists at a level that permits the downstream gated analysis to proceed, subject to denominator checks.
- `blocked_missing_covariates`: overlap evaluation could not proceed cleanly because required covariate fields are absent or unusable after explicit accounting.

# Section 3: A3 spec

- Analysis goal: direction-specific robustness audit for Task2 retrieval support without collapsing `C2G` and `G2C`.
- Allowed input files: S5 retrieval summary, S5 retrieval per-query, S5 chance identity check, S6 retrieval leaderboard.
- Required columns: `dataset`, `cell_line`, `direction`, `representation`, corrected retrieval metrics, gallery fields, `m_pos` fields, and chance-check fields keyed to the same retrieval slices.
- Legal grouping keys: `dataset`, `cell_line`, `direction`, `representation`, plus reviewed robustness-slice labels only.
- Required output schema: one row per lawful robustness slice with corrected retrieval fields, denominator/support fields, gallery and `m_pos` context, status fields, and provenance notes.
- Robustness logic to preserve: corrected retrieval is primary; gallery-size and `m_pos` context must remain explicit; direction-specific logic must remain explicit; chance-check coverage must be auditable.
- Things explicitly NOT to do: do not collapse `C2G` and `G2C`; do not replace corrected metrics with raw metrics; do not pre-implement unreviewed pairing logic; do not silently pool across datasets or cell lines.
- Failure modes / partial-output logic: emit `partial`, `unavailable`, or `blocked` rows when lawful robustness slices or denominator checks cannot be completed from current inputs.
- Validation checks: required-column assertions, key consistency across S5 and S6, chance-check coverage on the intended slice keys, and denominator conservation against source summaries.

# Section 4: B6 spec

- Analysis goal: audit whether exact dose/time covariate overlap exists in evaluable Task2 query pools before any controlled dose/time comparison is considered.
- Allowed input files: S5 retrieval per-query, LINCS `delta_meta.csv`, scPerturb `delta_meta.csv`.
- Required columns: S5 `dataset`, `cell_line`, `direction`, `representation`, `query_row_id`; metadata `row_id`, `time`, `dose_value`; and any explicit local-only label fields needed by the final audit output.
- Exact overlap definition: exact equality on `dataset`, `cell_line`, `time`, and `dose_value` after explicit `query_row_id -> row_id` metadata joins and explicit missing-covariate accounting.
- Allowed fallback behavior: missing `time` or `dose_value` may be carried forward only as explicit missing-covariate accounting; no silent binning, rounding, inferred harmonization, or hidden fallback rules are allowed.
- Required output schema: rows keyed by lawful audit slices with explicit overlap status, support counts, blocked-reason fields where needed, and any local-only label fields required by scope.
- Failure modes / partial-output logic: `no_overlap`, `partial_overlap`, `usable_overlap`, and `blocked_missing_covariates` are the allowed overlap-specific statuses; limitation-only outputs are allowed when overlap is not supportable.
- Stop condition: B6 may stop at a limitation-only output if exact shared strata are insufficient after explicit metadata joins, missing-covariate accounting, and exact overlap validation on dataset, cell_line, time, and dose_value.
- Validation checks: metadata join uniqueness, explicit accounting for missing covariates, denominator checks on evaluable-query pools, and exact-strata validation on the audited keys.

# Section 5: A1 spec

- Analysis goal: compare Task1 internal versus Task1 cross performance on the lawful common-scope slice without invoking ceiling language.
- Allowed input files: S1 leaderboard, S1 retrieval summary, S1 retrieval per-query when locally available and schema-compatible, S2 cross leaderboard, S2 cross retrieval summary, S2 cross retrieval per-query, and `task1_group_cross.parquet`.
- Common-scope definition: `Gene` and `Pathway` only; genetic-only cross logic only; internal rows must be restricted to lawful common-scope keys rather than broader internal scope.
- Legal comparison units: exact lawful paired units are preferred when valid per-query inputs exist; summary-based comparisons are allowed only as an explicit fallback when exact pairing cannot be lawfully supported from current frozen outputs.
- Expected output schema: rows with `internal_value`, `cross_value`, `delta_value`, denominator/support fields, status fields, and explicit notes describing pairing scope or pairing limitations.
- Handling if perfect pairing is unavailable: do not invent exact pairings or pseudo-pairs; degrade only to the highest-validity comparison supported by current frozen outputs and keep limitations explicit in output status and note fields.
- Contingency note: "If Task1 retrieval per-query inputs are absent locally or are schema-incompatible with the frozen pairing fields, A1 must degrade to the highest-validity summary-based comparison allowed by current frozen outputs, must carry explicit pairing limitations in output status and note fields, and must not invent exact pairings or summary-derived pseudo-pairs."
- Existing-output rule: existing outputs must be exhausted before any upstream S1 internal-group export is considered.
- Explicit prohibition: no ceiling language and no summary-derived claims that imply exact matched-unit support where none exists.
- Validation checks: required-column assertions, lawful common-scope filtering, denominator checks against the frozen summary tables, and explicit output labeling whenever exact pairing is unavailable.

# Section 6: A2 spec

- Analysis goal: shared-slice external-validity bridge comparison between Task1 and Task2 on lawful exact shared keys.
- Allowed input files: current downstream Task1 outputs, S4 Task2 group concordance, and S6 Task2 leaderboard outputs; any last-resort upstream touch remains outside this spec-freeze pass.
- Legal bridge keys: exact shared keys only, centered on `dataset`, `cell_line`, `target_token`, and `representation`, with any metric-family bridge rules reviewed explicitly during implementation.
- Bridge table semantics: rows represent lawful shared-slice bridge attempts only; no bridge rows may be invented and no non-identical keys may be collapsed into a forced bridge.
- Successful / unavailable / blocked rows: `successful` when both sides exist on lawful reviewed keys; `unavailable` when current downstream files do not support a row without violating rules; `blocked` when a required bridge rule, schema field, or lawful mapping remains unresolved.
- Planning note: "A2 is expected to have non-trivial exact shared-key coverage on current files, but lawful metric-family bridging remains unresolved and must be reviewed explicitly during implementation."
- Expected output schema: rows with Task1-side and Task2-side values where lawful, bridge status, denominator/support fields, and explicit scope and blocked-reason notes.
- No forced bridge rule: do not invent bridge rows, do not infer shared slices from summary-only coincidences, and do not turn unresolved metric-family mapping into hidden implementation defaults.
- Validation checks: required-column assertions, exact shared-key validation, denominator checks on any retained bridge rows, and explicit blocking whenever lawful bridge semantics are not yet frozen.

# Section 7: Implementation order

1. A3 skeleton review
2. B6 skeleton review
3. A1 skeleton review
4. A2 skeleton review
5. only then full implementation in a later pass

# Section 8: Review checklist

- Are keys lawful?
- Are filters lawful?
- Are stop conditions clear?
- Are partial outputs acceptable?
- Are local-only labels preserved?
- Are denominator checks sufficient?
