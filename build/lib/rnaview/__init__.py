from __future__ import annotations

from importlib import metadata as _metadata

from ._paths import basepars_dir, data_dir
from .api import analyze, from_structure, render_2d

try:
    __version__ = _metadata.version("rnaview")
except Exception:  # pragma: no cover
    __version__ = "0.0.0"

__all__ = [
    "__version__",
    "analyze",
    "basepars_dir",
    "data_dir",
    "from_structure",
    "render_2d",
]

