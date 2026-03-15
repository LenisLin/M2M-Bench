# Phase 1 Skeleton Review Note

## Files created

- `docs/manuscript_phase1_script_specs.md`
- `docs/manuscript_phase1_skeleton_review_note.md`

## Whether a common helper module was created

No helper module was created in this corrective pass. This pass only materializes the two documentation files needed to carry the frozen patch content.

## Unresolved assumptions per script

- `A3`: lawful robustness-slice definitions remain implementation-reviewed; corrected retrieval remains primary and direction-specific.
- `B6`: overlap-status logic is frozen, but exact support thresholds separating `partial_overlap` from `usable_overlap` still require explicit implementation review.
- `A1`: lawful exact pairing remains preferred; group-level internal support and any fallback beyond current frozen outputs remain review-sensitive.
- `A2`: exact shared keys are the only lawful bridge base, but metric-family bridging remains unresolved and must be reviewed explicitly during implementation.

## Places where lawful keys or schemas are still ambiguous

- `A3`: the final reviewed robustness-slice labels remain to be frozen during implementation review.
- `B6`: the exact audit-row layout may vary by representation or scope fields, but overlap logic and stop conditions are fixed.
- `A1`: summary-fallback labeling is fixed, but the exact output-row granularity under fallback remains implementation-reviewed.
- `A2`: lawful metric-family mapping remains the main unresolved schema-level ambiguity.

## Whether any script currently appears too large or too implicit

This corrective pass did not create or modify code files. The main human review risk remains implicit bridge semantics in `A2` and fallback labeling in `A1`, not script size.

## Recommended human review points before implementation begins

- Confirm that `A1` fallback behavior never implies exact pairing when only summary-level comparison is lawful.
- Confirm that `A2` bridge semantics remain exact-key-only and do not infer metric-family mappings silently.
- Confirm that `B6` implementation does not introduce bins, rounding, or hidden overlap rescue logic.
- Confirm that `A3` preserves corrected retrieval, direction separation, and gallery / `m_pos` context throughout.

## Provisional local observation

Provisional local observation: current inspected files may have little or no exact shared common-scope C2G/G2C dose/time overlap on dataset, cell_line, time, and dose_value. This is not a frozen result and must be revalidated during implementation.
