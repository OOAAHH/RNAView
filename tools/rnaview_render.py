#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


_BACKEND_VERSIONS = {
    "legacy": "legacy-shell-v1",
    "rustcore": "rustcore-shell-v1",
    "rustcore-release": "rustcore-release-shell-v1",
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _infer_format(path: Path) -> str | None:
    ext = path.suffix.lower()
    if ext == ".cif":
        return "cif"
    if ext in (".pdb", ".ent"):
        return "pdb"
    return None


def _json_dumps(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n"


def _run_cli_render(
    *,
    repo: Path,
    input_path: Path,
    out_dir: Path,
    renderer_bin: Path,
    want_ps_xml: bool,
    want_wrl: bool,
    log_path: Path,
) -> dict[str, Path]:
    fmt = _infer_format(input_path)
    if fmt not in ("pdb", "cif"):
        raise RuntimeError(f"unsupported input: {input_path}")

    try:
        rel = input_path.resolve().relative_to(repo.resolve())
        local_input = (out_dir / rel).resolve()
    except Exception:  # noqa: BLE001
        local_input = (out_dir / input_path.name).resolve()
    local_input.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_path, local_input)

    flags: list[str] = []
    if want_ps_xml:
        flags.append("-p")
    if want_wrl:
        flags.append("-v")
    flags.append("--pdb" if fmt == "pdb" else "--cif")

    env = dict(os.environ)
    env["RNAVIEW"] = str(repo)
    env.setdefault("USER", "unknown")
    env.setdefault("LOGNAME", env.get("USER", "unknown"))

    cmd = [str(renderer_bin), *flags, local_input.name]
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("wb") as log:
        proc = subprocess.run(
            cmd,
            cwd=str(local_input.parent),
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"render failed (code={proc.returncode}); see {log_path}")

    base = local_input.name
    outputs: dict[str, Path] = {}

    if want_ps_xml:
        ps = local_input.parent / f"{base}.ps"
        xml = local_input.parent / f"{base}.xml"
        if ps.exists() and ps.stat().st_size > 0:
            outputs["ps"] = ps
        if xml.exists() and xml.stat().st_size > 0:
            outputs["xml"] = xml

    if want_wrl:
        wrl_new = local_input.parent / f"{base}_new.wrl"
        wrl_old = local_input.parent / f"{base}.wrl"
        if wrl_new.exists() and wrl_new.stat().st_size > 0:
            outputs["wrl"] = wrl_new
        elif wrl_old.exists() and wrl_old.stat().st_size > 0:
            outputs["wrl"] = wrl_old

    return outputs


def _cmd_render(args: argparse.Namespace) -> int:
    repo = _repo_root()
    input_path = Path(args.input).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    requested_formats = [s.strip() for s in str(args.formats).split(",") if s.strip()]
    supported = {"ps", "xml", "wrl"}
    for f in requested_formats:
        if f not in supported:
            raise RuntimeError(f"unsupported format: {f} (supported: {sorted(supported)})")

    log_path = Path(args.log).resolve() if args.log else (out_dir / "render.log")

    backend = str(args.backend or "legacy")
    if backend not in _BACKEND_VERSIONS:
        raise RuntimeError(f"unsupported backend: {backend} (supported: {sorted(_BACKEND_VERSIONS)})")

    if args.rnaview_bin:
        renderer_bin = Path(args.rnaview_bin).resolve()
    else:
        if backend == "legacy":
            renderer_bin = (repo / "bin" / "rnaview").resolve()
        elif backend == "rustcore":
            renderer_bin = (repo / "bin" / "rnaview_rustcore").resolve()
        elif backend == "rustcore-release":
            renderer_bin = (repo / "bin" / "rnaview_rustcore_release").resolve()
        else:
            raise RuntimeError(f"unsupported backend: {backend}")

    if not (renderer_bin.exists() and os.access(renderer_bin, os.X_OK)):
        if backend == "legacy":
            build = repo / "tools" / "build_legacy_rnaview.sh"
        elif backend == "rustcore":
            build = repo / "tools" / "build_rnaview_rustcore.sh"
        else:
            build = repo / "tools" / "build_rnaview_rustcore_release.sh"
        raise RuntimeError(f"missing renderer binary: {renderer_bin} (build via {build})")

    want_ps_xml = any(f in requested_formats for f in ("ps", "xml"))
    want_wrl = "wrl" in requested_formats

    outputs = _run_cli_render(
        repo=repo,
        input_path=input_path,
        out_dir=out_dir,
        renderer_bin=renderer_bin,
        want_ps_xml=want_ps_xml,
        want_wrl=want_wrl,
        log_path=log_path,
    )

    missing = [f for f in requested_formats if f not in outputs]
    if missing:
        raise RuntimeError(f"render missing outputs: {missing} (see {log_path})")

    sys.stdout.write(
        _json_dumps(
            {
                "schema_version": 1,
                "status": "ok",
                "renderer_version": _BACKEND_VERSIONS[backend],
                "outputs": {k: str(v) for k, v in sorted(outputs.items())},
                "log": str(log_path),
            }
        )
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Phase 4 render runner (new pipeline entrypoint)")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("render", help="Render requested formats for one structure input")
    r.add_argument("--input", required=True, help="Input file (PDB/mmCIF)")
    r.add_argument("--out-dir", required=True, help="Output directory")
    r.add_argument("--formats", default="ps,xml,wrl", help="Comma-separated output formats (ps,xml,wrl)")
    r.add_argument(
        "--backend",
        default="legacy",
        help="Backend renderer (legacy|rustcore|rustcore-release). Use --rnaview-bin to override the binary path.",
    )
    r.add_argument("--rnaview-bin", default=None, help="Renderer binary path (overrides --backend default)")
    r.add_argument("--log", default=None, help="Write renderer log here (default: <out-dir>/render.log)")
    r.set_defaults(func=_cmd_render)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
