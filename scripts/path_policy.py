from __future__ import annotations

import os
from pathlib import Path


DEFAULT_NAS_ROOT = Path(os.getenv("M2M_NAS_ROOT", "/mnt/NAS_21T/ProjectData/M2M"))
DEFAULT_SNAPSHOTS_ROOT = Path(
    os.getenv("M2M_SNAPSHOTS_ROOT", str(DEFAULT_NAS_ROOT / "snapshots"))
)

DEFAULT_TASK1_SNAPSHOT_ROOT = Path(
    os.getenv(
        "M2M_TASK1_SNAPSHOT_ROOT",
        str(DEFAULT_SNAPSHOTS_ROOT / "task1_snapshot_v1"),
    )
)
DEFAULT_TASK1_SNAPSHOT_V2_CANDIDATE_ROOT = Path(
    os.getenv(
        "M2M_TASK1_SNAPSHOT_V2_CANDIDATE_ROOT",
        str(DEFAULT_SNAPSHOTS_ROOT / "task1_snapshot_v2_candidate"),
    )
)
DEFAULT_TASK2_LEGACY_SNAPSHOT_ROOT = Path(
    os.getenv(
        "M2M_TASK2_SNAPSHOT_V1_ROOT",
        str(DEFAULT_SNAPSHOTS_ROOT / "task2_snapshot_v1"),
    )
)
DEFAULT_TASK2_CORRECTED_SNAPSHOT_ROOT = Path(
    os.getenv(
        "M2M_TASK2_SNAPSHOT_V2_ROOT",
        str(DEFAULT_SNAPSHOTS_ROOT / "task2_snapshot_v2"),
    )
)
DEFAULT_TASK2_SNAPSHOT_V3_CANDIDATE_ROOT = Path(
    os.getenv(
        "M2M_TASK2_SNAPSHOT_V3_CANDIDATE_ROOT",
        str(DEFAULT_SNAPSHOTS_ROOT / "task2_snapshot_v3_candidate"),
    )
)

DEFAULT_OSMOSIS_LINCS_ROOT = Path(
    os.getenv(
        "M2M_OSMOSIS_LINCS_ROOT",
        "/mnt/NAS_21T/ProjectData/OSMOSIS/raw/CMap_LINCS/LINCS_level5",
    )
)
DEFAULT_OSMOSIS_SCPERTURB_ROOT = Path(
    os.getenv(
        "M2M_OSMOSIS_SCPERTURB_ROOT",
        "/mnt/NAS_21T/ProjectData/OSMOSIS/raw/scPerturb_Processed",
    )
)

DEFAULT_RUNS_ROOT = Path(os.getenv("M2M_RUNS_ROOT", str(DEFAULT_NAS_ROOT / "runs")))
DEFAULT_ARCHIVE_ROOT = Path(
    os.getenv("M2M_ARCHIVE_ROOT", str(DEFAULT_NAS_ROOT / "archive"))
)

DEFAULT_MANUSCRIPT_ACTIVE_ROOT = DEFAULT_RUNS_ROOT / "manuscript_active"
DEFAULT_MANUSCRIPT_ACTIVE_ANALYSIS_ROOT = DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "analysis"
DEFAULT_MANUSCRIPT_SUPPORT_ROOT = DEFAULT_RUNS_ROOT / "manuscript_support"
DEFAULT_MANUSCRIPT_SUPPORT_ANALYSIS_ROOT = DEFAULT_MANUSCRIPT_SUPPORT_ROOT / "analysis"
DEFAULT_MANUSCRIPT_SUPPORT_PLOT_READY_ROOT = DEFAULT_MANUSCRIPT_SUPPORT_ROOT / "plot_ready"
DEFAULT_MANUSCRIPT_SUPPORT_FIGURES_ROOT = DEFAULT_MANUSCRIPT_SUPPORT_ROOT / "figures"
DEFAULT_MANUSCRIPT_HISTORY_ROOT = DEFAULT_ARCHIVE_ROOT / "manuscript_history"
DEFAULT_MANUSCRIPT_PHASE1_ROOT = DEFAULT_RUNS_ROOT / "manuscript_phase1"


def resolve_path(project_root: Path, raw_path: str | Path) -> Path:
    path = Path(str(raw_path))
    if path.is_absolute():
        return path.resolve()
    return (project_root / path).resolve()
