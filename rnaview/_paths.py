from __future__ import annotations

from pathlib import Path


def data_dir() -> Path:
    # Package install directory; used as RNAVIEW root for locating BASEPARS.
    return Path(__file__).resolve().parent


def basepars_dir() -> Path:
    return data_dir() / "BASEPARS"


def embedded_hotcore_path() -> Path:
    return data_dir() / "_bin" / "rnaview-hotcore"

