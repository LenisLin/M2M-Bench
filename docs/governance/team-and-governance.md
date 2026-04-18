# Team And Governance

## Roles

### Human Lead

- owns benchmark scope, release sign-off, and final wording choices
- approves changes to task definitions, denominators, and figure claims

### Codex

- owns implementation, doc updates, verification commands, and manifest-aware
  delivery
- keeps benchmark semantics aligned with the contract docs

### Review Model

- reviews semantic consistency, metric logic, and wording drift against the
  contract docs and audited outputs

## Review Flow

1. Codex updates code or docs.
2. Verification commands confirm the claimed state.
3. Review checks focus on correctness, risk, and contract alignment.
4. Human Lead signs off on semantic changes.
