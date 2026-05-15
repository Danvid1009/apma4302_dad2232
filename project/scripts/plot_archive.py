"""Copy an existing plot PNG into ``visuals/archive/`` before regenerating."""

from __future__ import annotations

import shutil
from datetime import datetime, timezone
from pathlib import Path


def archive_prior_png(path: Path) -> Path | None:
    """
    If ``path`` exists, copy it to ``<parent>/archive/<stem>_YYYYMMDD_HHMMSSUTC<suffix>``.

    Returns the archive path if a copy was made, else None.
    """
    if not path.is_file():
        return None
    archive_dir = path.parent / "archive"
    archive_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%SUTC")
    dest = archive_dir / f"{path.stem}_{stamp}{path.suffix}"
    shutil.copy2(path, dest)
    print(f"previous plot preserved: {dest}")
    return dest
