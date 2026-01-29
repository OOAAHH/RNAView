#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import difflib
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _json_dumps(obj: Any, *, indent: int | None = 2) -> str:
    if indent is None:
        return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n"
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, indent=indent) + "\n"


def _has_glob_chars(text: str) -> bool:
    return any(c in text for c in ("*", "?", "["))


def _iter_list_file(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        yield s


def _collect_inputs(items: list[str], list_files: list[str]) -> list[Path]:
    repo = _repo_root()
    out: list[Path] = []
    allowed_exts = {".pdb", ".ent", ".cif"}
    excluded_name_suffixes = ("_tmp.pdb",)

    expanded: list[str] = []
    expanded.extend(items)
    for lf in list_files:
        expanded.extend(_iter_list_file(Path(lf)))

    for item in expanded:
        if _has_glob_chars(item):
            import glob

            matches = sorted(glob.glob(item, recursive=True))
            for m in matches:
                cand = Path(m)
                if cand.is_file():
                    if cand.suffix.lower() not in allowed_exts:
                        continue
                    if cand.name.lower().endswith(excluded_name_suffixes):
                        continue
                out.append(cand)
            continue

        p = Path(item)
        if not p.is_absolute():
            p = (repo / p).resolve() if (repo / p).exists() else p.resolve()
        if p.is_dir():
            for cand in sorted(p.rglob("*")):
                if cand.is_file():
                    if cand.suffix.lower() not in allowed_exts:
                        continue
                    if cand.name.lower().endswith(excluded_name_suffixes):
                        continue
                    out.append(cand)
            continue
        if p.is_file():
            if p.suffix.lower() not in allowed_exts:
                continue
            if p.name.lower().endswith(excluded_name_suffixes):
                continue
        out.append(p)

    seen: set[Path] = set()
    uniq: list[Path] = []
    for p in out:
        try:
            rp = p.resolve()
        except OSError:
            continue
        if rp in seen:
            continue
        seen.add(rp)
        uniq.append(rp)
    return uniq


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _sanitize_job_id(text: str) -> str:
    out_chars: list[str] = []
    for ch in text:
        if ch.isalnum() or ch in ("-", "_", "."):
            out_chars.append(ch)
        else:
            out_chars.append("_")
    out = "".join(out_chars).strip("._")
    return out or "job"


def _hash8(text: str) -> str:
    return hashlib.sha1(text.encode("utf-8", errors="replace")).hexdigest()[:8]


def _job_id_for_input(path: Path) -> str:
    base = _sanitize_job_id(path.stem)
    return f"{base}__{_hash8(str(path))}"


def _infer_format(path: Path) -> str | None:
    ext = path.suffix.lower()
    if ext == ".cif":
        return "cif"
    if ext in (".pdb", ".ent"):
        return "pdb"
    return None


def _relpath(repo: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(repo.resolve()).as_posix()
    except Exception:  # noqa: BLE001
        return str(path)


def _rel_to_test(repo: Path, path: Path) -> Path:
    rel = Path(_relpath(repo, path))
    if rel.parts and rel.parts[0] == "test":
        return Path(*rel.parts[1:])
    return Path("external") / _sanitize_job_id(path.name) / rel.name


def _gzip_write_deterministic(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename=path.name, mode="wb", compresslevel=9, fileobj=raw, mtime=0) as g:
            g.write(data)


def _gzip_read(path: Path) -> bytes:
    with gzip.open(path, "rb") as f:
        return f.read()


def _canonicalize_text(*, raw: bytes, fmt: str) -> bytes:
    # Canonicalize to LF text with no trailing whitespace, and strip unstable headers.
    if fmt == "gltf":
        obj = json.loads(raw.decode("utf-8", errors="replace"))
        return (json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n").encode("utf-8")

    text = raw.decode("utf-8", errors="replace")
    lines = text.splitlines()
    out_lines: list[str] = []

    for line in lines:
        s = line.rstrip()
        if fmt == "ps":
            if s.startswith("%%CreationDate:"):
                continue
        elif fmt == "wrl":
            if s.startswith("# Creation Date:") or s.startswith("# UserName:"):
                continue
        out_lines.append(s)

    return ("\n".join(out_lines) + "\n").encode("utf-8")


def _unified_diff(a: str, b: str, *, fromfile: str, tofile: str, max_lines: int) -> list[str]:
    diff = list(
        difflib.unified_diff(
            a.splitlines(),
            b.splitlines(),
            fromfile=fromfile,
            tofile=tofile,
            lineterm="",
        )
    )
    if len(diff) <= max_lines:
        return diff
    return diff[:max_lines] + [f"... (diff truncated; total lines={len(diff)})"]


def _load_allowlist(path: Path | None) -> dict[str, Any]:
    if path is None or not path.exists():
        return {"schema_version": 1, "approved": []}
    text = path.read_text(encoding="utf-8")
    stripped = text.strip()
    if not stripped:
        return {"schema_version": 1, "approved": []}

    try:
        data = json.loads(stripped)
    except json.JSONDecodeError:
        approved: list[Any] = []
        for line in text.splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith("-"):
                s = s[1:].strip()
            if s:
                approved.append(s)
        data = {"schema_version": 1, "approved": approved}

    if not isinstance(data, dict):
        return {"schema_version": 1, "approved": []}
    if "approved" not in data or not isinstance(data["approved"], list):
        data["approved"] = []
    if "schema_version" not in data:
        data["schema_version"] = 1
    return data


def _approved_ids(allowlist: dict[str, Any]) -> set[str]:
    out: set[str] = set()
    for item in allowlist.get("approved", []):
        if isinstance(item, str):
            out.add(item)
            continue
        if isinstance(item, dict):
            ident = item.get("id")
            if isinstance(ident, str) and ident:
                out.add(ident)
    return out


def _stable_id(*, input_id: str, fmt: str, semantics: str, renderer_version: str, payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    h = hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
    return f"{input_id}|render:{fmt}|{semantics}|{renderer_version}|{h}"


@dataclass(frozen=True)
class GoldenEntry:
    input: str
    input_sha256: str
    outputs: dict[str, str]  # fmt -> gz path (repo-relative)


def _default_rnaview_bin(repo: Path) -> Path:
    return repo / "bin" / "rnaview"


def _ensure_legacy_rnaview(repo: Path, rnaview_bin: Path) -> None:
    if rnaview_bin.exists() and os.access(rnaview_bin, os.X_OK):
        return
    build = repo / "tools" / "build_legacy_rnaview.sh"
    raise RuntimeError(f"missing legacy binary: {rnaview_bin} (build via {build})")


def _run_legacy_render(
    *,
    repo: Path,
    input_path: Path,
    work_dir: Path,
    rnaview_bin: Path,
    want_ps_xml: bool,
    want_wrl: bool,
) -> tuple[Path | None, Path | None, Path | None, Path]:
    fmt = _infer_format(input_path)
    if fmt not in ("pdb", "cif"):
        raise RuntimeError(f"unsupported input: {input_path}")

    # Copy the file into a per-input directory under work_dir, preserving test/ layout when possible,
    # so legacy outputs are local and stable (avoid absolute path effects).
    rel = Path(_relpath(repo, input_path))
    local_input = (work_dir / rel).resolve()
    local_input.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_path, local_input)

    flags: list[str] = []
    if want_ps_xml:
        flags.append("-p")
    if want_wrl:
        flags.append("-v")
    flags.append("--pdb" if fmt == "pdb" else "--cif")

    log_path = local_input.parent / "legacy_render.log"
    env = dict(os.environ)
    env["RNAVIEW"] = str(repo)
    # Ensure USER is always defined to avoid legacy VRML header issues.
    env.setdefault("USER", "unknown")
    env.setdefault("LOGNAME", env.get("USER", "unknown"))

    cmd = [str(rnaview_bin), *flags, local_input.name]
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
        raise RuntimeError(f"legacy render failed (code={proc.returncode}); see {log_path}")

    base = local_input.name
    ps_path = (local_input.parent / f"{base}.ps") if want_ps_xml else None
    xml_path = (local_input.parent / f"{base}.xml") if want_ps_xml else None
    wrl_candidates: list[Path] = []
    if want_wrl:
        wrl_candidates.append(local_input.parent / f"{base}_new.wrl")
        wrl_candidates.append(local_input.parent / f"{base}.wrl")
    wrl_path: Path | None = None
    for cand in wrl_candidates:
        if cand.exists() and cand.stat().st_size > 0:
            wrl_path = cand
            break

    return (
        ps_path if ps_path and ps_path.exists() else None,
        xml_path if xml_path and xml_path.exists() else None,
        wrl_path,
        log_path,
    )


def _run_command_render(
    *,
    repo: Path,
    input_path: Path,
    work_dir: Path,
    command: list[str],
    want_ps: bool,
    want_xml: bool,
    want_wrl: bool,
) -> tuple[Path | None, Path | None, Path | None, Path]:
    formats: list[str] = []
    if want_ps:
        formats.append("ps")
    if want_xml:
        formats.append("xml")
    if want_wrl:
        formats.append("wrl")

    render_dir = (work_dir / "candidate_render").resolve()
    render_dir.mkdir(parents=True, exist_ok=True)
    log_path = render_dir / "candidate_render.log"

    cmd = [
        *command,
        "--input",
        str(input_path),
        "--out-dir",
        str(render_dir),
        "--formats",
        ",".join(formats),
        "--log",
        str(log_path),
    ]

    env = dict(os.environ)
    env["RNAVIEW"] = str(repo)
    env.setdefault("USER", "unknown")
    env.setdefault("LOGNAME", env.get("USER", "unknown"))

    proc = subprocess.run(
        cmd,
        cwd=str(repo),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        try:
            log_path.write_text(proc.stdout, encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
        raise RuntimeError(f"candidate render command failed (code={proc.returncode}); see {log_path}")

    try:
        data = json.loads(proc.stdout or "{}")
    except json.JSONDecodeError as exc:
        try:
            log_path.write_text(proc.stdout, encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
        raise RuntimeError(f"candidate render returned non-JSON stdout: {exc}; see {log_path}") from exc

    if not isinstance(data, dict) or data.get("status") != "ok":
        try:
            log_path.write_text(proc.stdout, encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
        raise RuntimeError(f"candidate render returned error JSON; see {log_path}")

    outputs = data.get("outputs") or {}
    if not isinstance(outputs, dict):
        raise RuntimeError(f"candidate render JSON missing outputs map; see {log_path}")

    def get_path(key: str) -> Path | None:
        v = outputs.get(key)
        if not isinstance(v, str) or not v:
            return None
        p = Path(v)
        if not p.is_absolute():
            p = (render_dir / p).resolve()
        return p

    ps_path = get_path("ps") if want_ps else None
    xml_path = get_path("xml") if want_xml else None
    wrl_path = get_path("wrl") if want_wrl else None

    return ps_path, xml_path, wrl_path, log_path


def _cmd_freeze(args: argparse.Namespace) -> int:
    repo = _repo_root()
    golden_dir = Path(args.golden_dir) if args.golden_dir else (repo / "test" / "golden_render")
    rnaview_bin = Path(args.rnaview_bin) if args.rnaview_bin else _default_rnaview_bin(repo)

    _ensure_legacy_rnaview(repo, rnaview_bin)

    requested_formats = [s.strip() for s in str(args.formats).split(",") if s.strip()]
    supported = {"ps", "xml", "wrl", "svg", "gltf"}
    for f in requested_formats:
        if f not in supported:
            raise RuntimeError(f"unsupported format in --formats: {f} (supported: {sorted(supported)})")

    svg_converter_version: str | None = None
    if "svg" in requested_formats:
        import rnaview_ps_svg

        svg_converter_version = rnaview_ps_svg.RENDERER_VERSION

    gltf_converter_version: str | None = None
    if "gltf" in requested_formats:
        import rnaview_vrml_gltf

        gltf_converter_version = rnaview_vrml_gltf.RENDERER_VERSION

    inputs = _collect_inputs(list(args.inputs), list(args.list_file or []))
    if not inputs:
        raise RuntimeError("no inputs selected")

    started = time.time()
    entries: list[dict[str, Any]] = []
    counts = {"ok": 0, "failed": 0, "written": 0}

    def one(input_path: Path) -> dict[str, Any]:
        input_rel = _relpath(repo, input_path)
        input_sha = _sha256_path(input_path)
        job_id = _job_id_for_input(input_path)

        with tempfile.TemporaryDirectory(prefix="rnaview_gate_d_freeze.") as tmp:
            work = Path(tmp)
            want_ps_xml = any(f in requested_formats for f in ("ps", "xml", "svg"))
            want_wrl = any(f in requested_formats for f in ("wrl", "gltf"))

            ps_path, xml_path, wrl_path, log_rel = _run_legacy_render(
                repo=repo,
                input_path=input_path,
                work_dir=work,
                rnaview_bin=rnaview_bin,
                want_ps_xml=want_ps_xml,
                want_wrl=want_wrl,
            )

            out_map: dict[str, str] = {}
            produced: dict[str, str] = {}

            rel_to_test = _rel_to_test(repo, input_path)
            base_out = (golden_dir / rel_to_test).resolve()

            def write_one(fmt: str, *, raw: bytes, produced_name: str) -> None:
                canon = _canonicalize_text(raw=raw, fmt=fmt)
                out_path = base_out.parent / f"{base_out.name}.{fmt}.canon.gz"
                _gzip_write_deterministic(out_path, canon)
                out_map[fmt] = _relpath(repo, out_path)
                produced[fmt] = produced_name

            if "ps" in requested_formats:
                if ps_path is None:
                    raise RuntimeError("missing legacy output for ps")
                write_one("ps", raw=ps_path.read_bytes(), produced_name=ps_path.name)
            if "xml" in requested_formats:
                if xml_path is None:
                    raise RuntimeError("missing legacy output for xml")
                write_one("xml", raw=xml_path.read_bytes(), produced_name=xml_path.name)
            if "wrl" in requested_formats:
                if wrl_path is None:
                    raise RuntimeError("missing legacy output for wrl")
                write_one("wrl", raw=wrl_path.read_bytes(), produced_name=wrl_path.name)
            if "svg" in requested_formats:
                if ps_path is None:
                    raise RuntimeError("missing legacy output for svg (needs ps)")
                import rnaview_ps_svg

                svg_text = rnaview_ps_svg.svg_from_ps(ps_path.read_text(encoding="utf-8", errors="replace"))
                write_one(
                    "svg",
                    raw=svg_text.encode("utf-8"),
                    produced_name=f"{ps_path.name}->{rnaview_ps_svg.RENDERER_VERSION}",
                )
            if "gltf" in requested_formats:
                if wrl_path is None:
                    raise RuntimeError("missing legacy output for gltf (needs wrl)")
                import rnaview_vrml_gltf

                gltf_text = rnaview_vrml_gltf.gltf_from_vrml(
                    wrl_path.read_text(encoding="utf-8", errors="replace")
                )
                write_one(
                    "gltf",
                    raw=gltf_text.encode("utf-8"),
                    produced_name=f"{wrl_path.name}->{rnaview_vrml_gltf.RENDERER_VERSION}",
                )

            return {
                "status": "ok",
                "job_id": job_id,
                "input": input_rel,
                "input_sha256": input_sha,
                "outputs": out_map,
                "produced_names": produced,
            }

    workers = int(args.workers or 1)
    if workers <= 1:
        for p in inputs:
            try:
                e = one(p)
                entries.append(e)
                counts["ok"] += 1
                counts["written"] += 1
            except Exception as exc:  # noqa: BLE001
                counts["failed"] += 1
                entries.append(
                    {
                        "status": "failed",
                        "job_id": _job_id_for_input(p),
                        "input": _relpath(repo, p),
                        "input_sha256": None,
                        "outputs": {},
                        "error": str(exc),
                    }
                )
                if not args.keep_going:
                    break
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(one, p): p for p in inputs}
            for fut in concurrent.futures.as_completed(futs):
                p = futs[fut]
                try:
                    e = fut.result()
                    counts["ok"] += 1
                    counts["written"] += 1
                    entries.append(e)
                except Exception as exc:  # noqa: BLE001
                    counts["failed"] += 1
                    entries.append(
                        {
                            "status": "failed",
                            "job_id": _job_id_for_input(p),
                            "input": _relpath(repo, p),
                            "input_sha256": None,
                            "outputs": {},
                            "error": str(exc),
                        }
                    )
                    if not args.keep_going:
                        break

    manifest = {
        "schema_version": 1,
        "canonical_version": 1,
        "formats": requested_formats,
        "baseline": {"engine": "legacy", "rnaview_bin": _relpath(repo, rnaview_bin)},
        "converters": {
            k: v
            for k, v in {
                "svg": svg_converter_version,
                "gltf": gltf_converter_version,
            }.items()
            if v
        },
        "elapsed_ms": int((time.time() - started) * 1000),
        "counts": counts,
        "entries": sorted(entries, key=lambda e: str(e.get("input") or "")),
    }

    golden_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = golden_dir / "manifest.json"
    manifest_path.write_text(_json_dumps(manifest, indent=2), encoding="utf-8")

    if counts["failed"]:
        sys.stderr.write(f"freeze: ok={counts['ok']} failed={counts['failed']} (see {manifest_path})\n")
        return 1
    sys.stderr.write(f"freeze: ok={counts['ok']} (manifest={manifest_path})\n")
    return 0


def _load_golden_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _cmd_compare(args: argparse.Namespace) -> int:
    repo = _repo_root()
    out_dir = Path(args.out_dir) if args.out_dir else (repo / "out_phase4_gate_d")
    golden_dir = Path(args.golden_dir) if args.golden_dir else (repo / "test" / "golden_render")
    manifest_path = Path(args.manifest) if args.manifest else (golden_dir / "manifest.json")
    allowlist_path = Path(args.allowlist) if args.allowlist else (repo / "test" / "gate_d_allowlist.yaml")
    rnaview_bin = Path(args.rnaview_bin) if args.rnaview_bin else _default_rnaview_bin(repo)
    candidate_engine = str(args.candidate_engine or "command")

    semantics = str(args.semantics or "legacy-v1")
    renderer_version = str(args.renderer_version or "dev")

    allowlist = _load_allowlist(allowlist_path)
    approved = _approved_ids(allowlist)

    manifest = _load_golden_manifest(manifest_path)
    formats = list(manifest.get("formats") or [])
    if not formats:
        raise RuntimeError(f"invalid render golden manifest (no formats): {manifest_path}")

    svg_converter_version = None
    if "svg" in formats:
        import rnaview_ps_svg

        svg_converter_version = rnaview_ps_svg.RENDERER_VERSION

    gltf_converter_version = None
    if "gltf" in formats:
        import rnaview_vrml_gltf

        gltf_converter_version = rnaview_vrml_gltf.RENDERER_VERSION

    # Optional filtering: if inputs are passed, only compare those. Otherwise compare all manifest entries.
    wanted: set[str] = set()
    if args.inputs or args.list_file:
        selected = _collect_inputs(list(args.inputs), list(args.list_file or []))
        wanted = {_relpath(repo, p) for p in selected}

    entries = [e for e in (manifest.get("entries") or []) if isinstance(e, dict)]
    if wanted:
        entries = [e for e in entries if str(e.get("input") or "") in wanted]

    if not entries:
        raise RuntimeError("no entries selected for compare")

    if candidate_engine == "legacy":
        _ensure_legacy_rnaview(repo, rnaview_bin)
    elif candidate_engine == "command":
        if args.candidate_cmd:
            candidate_cmd = list(args.candidate_cmd)
            if candidate_cmd and candidate_cmd[0] == "--":
                candidate_cmd = candidate_cmd[1:]
            if not candidate_cmd:
                raise RuntimeError("--candidate-cmd is empty (it must be last and include a command)")
        else:
            candidate_cmd = [sys.executable, str(repo / "tools" / "rnaview_render.py"), "render"]
    else:
        raise RuntimeError(f"unsupported --candidate-engine: {candidate_engine}")

    if out_dir.exists():
        shutil.rmtree(out_dir)
    (out_dir / "cases").mkdir(parents=True, exist_ok=True)

    started = time.time()
    results: list[dict[str, Any]] = []
    counts = {"ok": 0, "changed": 0, "unapproved": 0, "failed": 0}

    for entry in entries:
        input_rel = str(entry.get("input") or "")
        input_path = (repo / input_rel).resolve() if input_rel else None
        job_id = str(entry.get("job_id") or _job_id_for_input(Path(input_rel) if input_rel else Path("job")))
        case_dir = out_dir / "cases" / job_id
        case_dir.mkdir(parents=True, exist_ok=True)

        if input_path is None or not input_path.exists():
            counts["failed"] += 1
            results.append({"input": input_rel, "job_id": job_id, "status": "failed", "error": "missing input"})
            continue

        try:
            with tempfile.TemporaryDirectory(prefix="rnaview_gate_d_compare.") as tmp:
                work = Path(tmp)
                want_ps = "ps" in formats
                want_xml = "xml" in formats
                want_wrl = any(f in formats for f in ("wrl", "gltf"))
                if "svg" in formats:
                    want_ps = True
                if candidate_engine == "command":
                    ps_path, xml_path, wrl_path, log_rel = _run_command_render(
                        repo=repo,
                        input_path=input_path,
                        work_dir=work,
                        command=candidate_cmd,
                        want_ps=want_ps,
                        want_xml=want_xml,
                        want_wrl=want_wrl,
                    )
                elif candidate_engine == "legacy":
                    ps_path, xml_path, wrl_path, log_rel = _run_legacy_render(
                        repo=repo,
                        input_path=input_path,
                        work_dir=work,
                        rnaview_bin=rnaview_bin,
                        want_ps_xml=want_ps or want_xml,
                        want_wrl=want_wrl,
                    )

                produced_map: dict[str, Path | None] = {
                    "ps": ps_path,
                    "xml": xml_path,
                    "wrl": wrl_path,
                    "svg": None,
                    "gltf": None,
                }
                if "svg" in formats:
                    if ps_path is None:
                        raise RuntimeError("missing candidate ps needed for svg")
                    import rnaview_ps_svg

                    svg_text = rnaview_ps_svg.svg_from_ps(ps_path.read_text(encoding="utf-8", errors="replace"))
                    svg_path = case_dir / "candidate.svg"
                    svg_path.write_text(svg_text, encoding="utf-8")
                    produced_map["svg"] = svg_path
                if "gltf" in formats:
                    if wrl_path is None:
                        raise RuntimeError("missing candidate wrl needed for gltf")
                    import rnaview_vrml_gltf

                    gltf_text = rnaview_vrml_gltf.gltf_from_vrml(
                        wrl_path.read_text(encoding="utf-8", errors="replace")
                    )
                    gltf_path = case_dir / "candidate.gltf"
                    gltf_path.write_text(gltf_text, encoding="utf-8")
                    produced_map["gltf"] = gltf_path

                events: list[dict[str, Any]] = []
                diffs_by_fmt: dict[str, list[str]] = {}
                unapproved_ids: list[str] = []

                for fmt in formats:
                    renderer_version_fmt = renderer_version
                    if fmt == "svg" and svg_converter_version:
                        renderer_version_fmt = f"{renderer_version}+{svg_converter_version}"
                    if fmt == "gltf" and gltf_converter_version:
                        renderer_version_fmt = f"{renderer_version}+{gltf_converter_version}"

                    golden_gz = entry.get("outputs", {}).get(fmt)
                    if not isinstance(golden_gz, str) or not golden_gz:
                        continue
                    golden_path = (repo / golden_gz).resolve()
                    if not golden_path.exists():
                        raise RuntimeError(f"missing golden: {golden_gz}")
                    golden_canon = _gzip_read(golden_path)
                    golden_sha = _sha256_bytes(golden_canon)

                    src = produced_map.get(fmt)
                    if src is None:
                        payload = {"change": "missing", "golden_sha256": golden_sha, "candidate_sha256": None}
                        ident = _stable_id(
                            input_id=input_rel,
                            fmt=fmt,
                            semantics=semantics,
                            renderer_version=renderer_version_fmt,
                            payload=payload,
                        )
                        ev = {"id": ident, "format": fmt, "change": "missing", "payload": payload}
                        ev["approved"] = ident in approved
                        events.append(ev)
                        if not ev["approved"]:
                            unapproved_ids.append(ident)
                        continue

                    candidate_raw = src.read_bytes()
                    candidate_canon = _canonicalize_text(raw=candidate_raw, fmt=fmt)
                    cand_sha = _sha256_bytes(candidate_canon)
                    if candidate_canon == golden_canon:
                        continue

                    payload = {
                        "change": "changed",
                        "golden_sha256": golden_sha,
                        "candidate_sha256": cand_sha,
                    }
                    ident = _stable_id(
                        input_id=input_rel,
                        fmt=fmt,
                        semantics=semantics,
                        renderer_version=renderer_version_fmt,
                        payload=payload,
                    )

                    diff_lines = _unified_diff(
                        golden_canon.decode("utf-8", errors="replace"),
                        candidate_canon.decode("utf-8", errors="replace"),
                        fromfile=golden_gz,
                        tofile=_relpath(repo, src),
                        max_lines=int(args.max_diff_lines or 200),
                    )
                    diffs_by_fmt[fmt] = diff_lines

                    ev = {
                        "id": ident,
                        "format": fmt,
                        "change": "changed",
                        "payload": payload,
                        "diff": diff_lines,
                    }
                    ev["approved"] = ident in approved
                    events.append(ev)
                    if not ev["approved"]:
                        unapproved_ids.append(ident)

                diff_path = case_dir / "diff.json"
                diff_path.write_text(_json_dumps({"schema_version": 1, "events": events}, indent=2), encoding="utf-8")

                log_dst = case_dir / "candidate_render.log"
                try:
                    shutil.copy2(log_rel, log_dst)
                    log_out = _relpath(repo, log_dst)
                except Exception:  # noqa: BLE001
                    log_out = None

                if unapproved_ids:
                    counts["unapproved"] += 1
                    status = "unapproved"
                elif events:
                    counts["changed"] += 1
                    status = "changed"
                else:
                    counts["ok"] += 1
                    status = "ok"

                results.append(
                    {
                        "input": input_rel,
                        "job_id": job_id,
                        "status": status,
                        "events": len(events),
                        "unapproved_events": len(unapproved_ids),
                        "log": log_out,
                        "diff": _relpath(repo, diff_path),
                    }
                )
        except Exception as exc:  # noqa: BLE001
            counts["failed"] += 1
            results.append({"input": input_rel, "job_id": job_id, "status": "failed", "error": str(exc)})

    summary = {
        "schema_version": 1,
        "baseline_manifest": _relpath(repo, manifest_path),
        "allowlist": _relpath(repo, allowlist_path),
        "candidate_engine": candidate_engine,
        "candidate_cmd": candidate_cmd if candidate_engine == "command" else None,
        "semantics": semantics,
        "renderer_version": renderer_version,
        "elapsed_ms": int((time.time() - started) * 1000),
        "counts": counts,
        "results": results,
    }
    (out_dir / "summary.json").write_text(_json_dumps(summary, indent=2), encoding="utf-8")

    # Minimal human report.
    lines: list[str] = []
    lines.append("# Gate D report")
    lines.append(f"- semantics: `{semantics}`")
    lines.append(f"- renderer_version: `{renderer_version}`")
    lines.append(f"- candidate_engine: `{candidate_engine}`")
    if candidate_engine == "command":
        lines.append(f"- candidate_cmd: `{candidate_cmd}`")
    lines.append(f"- baseline manifest: `{_relpath(repo, manifest_path)}`")
    lines.append(f"- allowlist: `{_relpath(repo, allowlist_path)}`")
    lines.append(f"- counts: `{counts}`")
    lines.append("")
    for k in ("failed", "unapproved", "changed"):
        subset = [r for r in results if r.get("status") == k]
        if not subset:
            continue
        lines.append(f"## {k.title()}")
        lines.append("")
        for r in subset:
            lines.append(f"- {k}: `{r.get('input')}` job_id=`{r.get('job_id')}` events={r.get('events')}")
        lines.append("")
    (out_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    if counts["failed"] or counts["unapproved"]:
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Gate D: render regressions (canonicalize then diff-0)")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_freeze = sub.add_parser("freeze", help="Run legacy render and write canonical golden outputs")
    p_freeze.add_argument("inputs", nargs="*", default=["test/pdb", "test/mmcif"])
    p_freeze.add_argument("--list-file", action="append", default=[], help="Text file listing extra inputs")
    p_freeze.add_argument("--golden-dir", default=None, help="Where to write golden files (default: test/golden_render)")
    p_freeze.add_argument("--rnaview-bin", default=None, help="Path to legacy rnaview binary (default: bin/rnaview)")
    p_freeze.add_argument(
        "--formats",
        default="ps,xml,wrl",
        help="Comma-separated formats to freeze (supported: ps,xml,svg,wrl,gltf; default: ps,xml,wrl)",
    )
    p_freeze.add_argument("--workers", type=int, default=1, help="Parallel workers (default: 1)")
    p_freeze.add_argument("--keep-going", action="store_true", help="Continue after failures")
    p_freeze.set_defaults(func=_cmd_freeze)

    p_cmp = sub.add_parser("compare", help="Run a renderer and compare against render goldens")
    p_cmp.add_argument("inputs", nargs="*", default=[])
    p_cmp.add_argument("--list-file", action="append", default=[], help="Text file listing extra inputs")
    p_cmp.add_argument("--out-dir", default=None, help="Where to write the Gate D report (default: out_phase4_gate_d)")
    p_cmp.add_argument("--golden-dir", default=None, help="Golden dir (default: test/golden_render)")
    p_cmp.add_argument("--manifest", default=None, help="Golden manifest path (default: <golden_dir>/manifest.json)")
    p_cmp.add_argument("--allowlist", default=None, help="Allowlist path (default: test/gate_d_allowlist.yaml)")
    p_cmp.add_argument("--rnaview-bin", default=None, help="Legacy rnaview binary for candidate run (default: bin/rnaview)")
    p_cmp.add_argument(
        "--candidate-engine",
        choices=["command", "legacy"],
        default="command",
        help="How to run the candidate renderer (default: command)",
    )
    p_cmp.add_argument(
        "--candidate-cmd",
        nargs=argparse.REMAINDER,
        default=None,
        help="Candidate render command (must be last; default: python3 tools/rnaview_render.py render)",
    )
    p_cmp.add_argument(
        "--semantics",
        default="legacy-v1",
        help="Semantics label to include in diff event ids (default: legacy-v1)",
    )
    p_cmp.add_argument(
        "--renderer-version",
        default="dev",
        help="Renderer version label to include in diff event ids (default: dev)",
    )
    p_cmp.add_argument("--max-diff-lines", type=int, default=200, help="Max unified diff lines per format")
    p_cmp.set_defaults(func=_cmd_compare)

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
