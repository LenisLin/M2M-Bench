# M2M-Bench Version

Version: v0.4.0

Reasoning:
- `pyproject.toml` declares `name = "m2m-bench"` and `version = "0.4.0"`.
- The older reboot note that said "Task2 not yet frozen" no longer matches the current local source of truth.
- `docs/contracts/task2_spec.md` now locks corrected Task2 to a multisource contract spanning `LINCS` and `scPerturb`.
- Local audited run evidence covers Task1 `S0-S2`, legacy/interim Task2 `S3`, the Task2 K562 FM extractor family, and corrected multisource Task2 `S3-S6`.
- `scripts/s7_benchmark_health_and_story_pack.py` exists, but no audited S7 run artifact was found, so `v0.4.0` does not claim S7 completion.
