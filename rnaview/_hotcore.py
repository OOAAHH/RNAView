from __future__ import annotations

import os
import shutil
import stat
import subprocess
import sys
import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from ._paths import data_dir, embedded_hotcore_path


class HotcoreNotFoundError(FileNotFoundError):
    pass


@dataclass(frozen=True)
class HotcorePaths:
    exe: Path
    rnaview_root: Path


def _default_cache_dir() -> Path:
    xdg = os.environ.get("XDG_CACHE_HOME")
    if xdg:
        return Path(xdg)
    return Path.home() / ".cache"


def _ensure_executable(path: Path) -> None:
    try:
        st = path.stat()
    except OSError:
        return
    if st.st_mode & stat.S_IXUSR:
        return
    try:
        path.chmod(st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    except OSError:
        # We may not be able to chmod in site-packages; caller can fall back to a cache copy.
        pass


def _prepare_embedded_hotcore(path: Path) -> Path:
    _ensure_executable(path)
    try:
        if os.access(path, os.X_OK):
            return path
    except OSError:
        pass

    # Fall back to a user-writable cache location.
    cache_dir = _default_cache_dir() / "rnaview" / "bin"
    cache_dir.mkdir(parents=True, exist_ok=True)
    dst = cache_dir / path.name
    try:
        src_mtime = path.stat().st_mtime
        dst_mtime = dst.stat().st_mtime if dst.exists() else -1
        if dst_mtime < src_mtime:
            shutil.copy2(path, dst)
    except OSError:
        shutil.copy2(path, dst)
    _ensure_executable(dst)
    return dst


def find_hotcore() -> HotcorePaths:
    override = os.environ.get("RNAVIEW_HOTCORE")
    if override:
        p = Path(override)
        if p.exists():
            return HotcorePaths(exe=p, rnaview_root=data_dir())

    embedded = embedded_hotcore_path()
    if embedded.exists():
        if not sys.platform.startswith("linux"):
            raise HotcoreNotFoundError(
                f"embedded rnaview-hotcore is only supported on Linux; sys.platform={sys.platform!r}"
            )
        mach = platform.machine().lower()
        if mach not in ("x86_64", "amd64"):
            raise HotcoreNotFoundError(f"embedded rnaview-hotcore is only supported on x86_64; machine={mach!r}")
        return HotcorePaths(exe=_prepare_embedded_hotcore(embedded), rnaview_root=data_dir())

    which = shutil.which("rnaview-hotcore")
    if which:
        return HotcorePaths(exe=Path(which), rnaview_root=data_dir())

    raise HotcoreNotFoundError(
        "rnaview-hotcore not found. Install a wheel that bundles the binary, "
        "or set RNAVIEW_HOTCORE, or add rnaview-hotcore to PATH."
    )


def run_hotcore(
    args: Sequence[str],
    *,
    cwd: str | os.PathLike[str] | None = None,
    env: dict[str, str] | None = None,
    check: bool = True,
    capture_output: bool = False,
    text: bool = True,
) -> subprocess.CompletedProcess[str]:
    hc = find_hotcore()
    merged_env = dict(os.environ)
    merged_env.setdefault("RNAVIEW", str(hc.rnaview_root))
    if env:
        merged_env.update(env)

    cmd = [str(hc.exe), *args]
    return subprocess.run(
        cmd,
        cwd=cwd,
        env=merged_env,
        check=check,
        capture_output=capture_output,
        text=text,
    )
