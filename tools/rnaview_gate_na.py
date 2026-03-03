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


def _sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


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


def _canonical_core_for_regress(core: Any) -> Any:
    # Keep NA gate focused on scientific content; renderer/writer bookkeeping should be covered
    # by the `.out`/render artifacts. See doc/spec.md for `out_index`.
    if not isinstance(core, dict):
        return core
    base_pairs = core.get("base_pairs")
    if isinstance(base_pairs, list):
        for bp in base_pairs:
            if isinstance(bp, dict):
                bp.pop("out_index", None)
    return core


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


def _sysroot_env() -> dict[str, str]:
    sysroot = Path.home() / ".cache" / "rnaview-toolchain" / "root"
    if not sysroot.exists():
        return {}
    paths = [
        sysroot / "usr" / "lib" / "x86_64-linux-gnu",
        sysroot / "lib" / "x86_64-linux-gnu",
        sysroot / "usr" / "lib",
        sysroot / "lib",
    ]
    ld = ":".join(str(p) for p in paths if p.exists())
    if not ld:
        return {}
    env: dict[str, str] = {}
    prev = os.environ.get("LD_LIBRARY_PATH")
    env["LD_LIBRARY_PATH"] = f"{ld}:{prev}" if prev else ld
    return env


def _find_rust_hotcore_binary(repo: Path) -> Path | None:
    release = repo / "rust" / "target" / "release" / "rnaview-hotcore"
    debug = repo / "rust" / "target" / "debug" / "rnaview-hotcore"

    found: list[tuple[float, Path]] = []
    for p in (release, debug):
        if not (p.exists() and os.access(p, os.X_OK)):
            continue
        try:
            mtime = float(p.stat().st_mtime)
        except OSError:
            mtime = 0.0
        found.append((mtime, p))

    if not found:
        return None

    found.sort(key=lambda t: t[0], reverse=True)
    return found[0][1]


def _rust_hotcore_needs_rebuild(repo: Path, binary: Path) -> bool:
    try:
        bin_mtime = binary.stat().st_mtime
    except OSError:
        return True

    candidates: list[Path] = [
        repo / "rust" / "Cargo.toml",
        repo / "rust" / "Cargo.lock",
    ]
    src_dir = repo / "rust" / "src"
    if src_dir.exists():
        candidates.extend(sorted(src_dir.rglob("*.rs")))

    for p in candidates:
        try:
            if p.exists() and p.stat().st_mtime > bin_mtime:
                return True
        except OSError:
            continue

    return False


def _ensure_rust_hotcore_binary(repo: Path) -> Path:
    found = _find_rust_hotcore_binary(repo)
    if found is not None and not _rust_hotcore_needs_rebuild(repo, found):
        return found

    cargo = shutil.which("cargo")
    if cargo is None:
        raise RuntimeError("missing cargo; cannot build rnaview-hotcore")

    manifest = repo / "rust" / "Cargo.toml"
    use_sysroot = shutil.which("cc") is None
    if use_sysroot:
        cmd = ["bash", str(repo / "tools" / "cargo_sysroot.sh"), "build", "--manifest-path", str(manifest)]
    else:
        cmd = [cargo, "build", "--manifest-path", str(manifest)]

    proc = subprocess.run(cmd, cwd=str(repo), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"failed to build rnaview-hotcore (code={proc.returncode}):\n{proc.stdout}")

    found = _find_rust_hotcore_binary(repo)
    if found is None:
        raise RuntimeError("build succeeded but rnaview-hotcore not found under rust/target/(debug|release)")
    return found


@dataclass(frozen=True)
class Produced:
    core: Any
    out_bytes: bytes
    ps_bytes: bytes
    xml_bytes: bytes
    wrl_bytes: bytes


def _run_hotcore_pipeline(
    *,
    repo: Path,
    rust_bin: Path,
    input_path: Path,
    semantics: str,
) -> Produced:
    fmt = _infer_format(input_path)
    if fmt not in ("pdb", "cif"):
        raise RuntimeError(f"unsupported input: {input_path}")

    env = dict(os.environ)
    env["RNAVIEW"] = str(repo)
    env.update(_sysroot_env())
    env.setdefault("USER", "unknown")
    env.setdefault("LOGNAME", env.get("USER", "unknown"))

    with tempfile.TemporaryDirectory(prefix="rnaview_gate_na.") as td:
        work = Path(td)
        pairs_path = work / "pairs.json"
        out_path = work / "engine.out"
        ps_path = work / "render.ps"
        xml_path = work / "render.xml"
        wrl_path = work / "render.wrl"

        cmd_compute = [
            str(rust_bin),
            "from-structure",
            str(input_path),
            "--format",
            fmt,
            "--oracle",
            "compute",
            "--semantics",
            semantics,
            "-o",
            str(pairs_path),
            "--emit-out",
            str(out_path),
        ]
        proc = subprocess.run(
            cmd_compute,
            cwd=str(repo),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"from-structure failed (code={proc.returncode}):\n{proc.stdout}")

        cmd_2d = [
            str(rust_bin),
            "render-2d",
            str(pairs_path),
            "--source",
            str(input_path),
            "--semantics",
            semantics,
            "--out-xml",
            str(xml_path),
            "--out-ps",
            str(ps_path),
        ]
        proc = subprocess.run(
            cmd_2d,
            cwd=str(repo),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"render-2d failed (code={proc.returncode}):\n{proc.stdout}")

        cmd_wrl = [
            str(rust_bin),
            "render-wrl",
            str(pairs_path),
            "--source",
            str(input_path),
            "--semantics",
            semantics,
            "-o",
            str(wrl_path),
        ]
        proc = subprocess.run(
            cmd_wrl,
            cwd=str(repo),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"render-wrl failed (code={proc.returncode}):\n{proc.stdout}")

        pairs = json.loads(pairs_path.read_text(encoding="utf-8", errors="replace"))
        core = _canonical_core_for_regress(pairs.get("core"))

        return Produced(
            core=core,
            out_bytes=out_path.read_bytes(),
            ps_bytes=ps_path.read_bytes(),
            xml_bytes=xml_path.read_bytes(),
            wrl_bytes=wrl_path.read_bytes(),
        )


@dataclass(frozen=True)
class FreezeResult:
    input: Path
    job_id: str
    status: str  # ok|failed
    input_sha256: str | None
    outputs: dict[str, Path] | None
    error: str | None = None


def _freeze_one(
    *,
    repo: Path,
    rust_bin: Path,
    input_path: Path,
    semantics: str,
    golden_dir: Path,
    formats: set[str],
) -> FreezeResult:
    job_id = _job_id_for_input(input_path)
    try:
        sha = _sha256_path(input_path)
    except Exception:  # noqa: BLE001
        sha = None

    try:
        produced = _run_hotcore_pipeline(repo=repo, rust_bin=rust_bin, input_path=input_path, semantics=semantics)
    except Exception as e:  # noqa: BLE001
        return FreezeResult(input=input_path, job_id=job_id, status="failed", input_sha256=sha, outputs=None, error=str(e))

    rel_to_test = _rel_to_test(repo, input_path)
    base_out = (golden_dir / rel_to_test).resolve()

    out_paths: dict[str, Path] = {}

    if "core" in formats:
        core_path = base_out.parent / f"{base_out.name}.core.json"
        core_path.parent.mkdir(parents=True, exist_ok=True)
        core_path.write_text(_json_dumps(produced.core, indent=None), encoding="utf-8")
        out_paths["core"] = core_path

    def write_gz(fmt: str, raw: bytes) -> None:
        canon = _canonicalize_text(raw=raw, fmt=fmt)
        out_path = base_out.parent / f"{base_out.name}.{fmt}.canon.gz"
        _gzip_write_deterministic(out_path, canon)
        out_paths[fmt] = out_path

    if "out" in formats:
        write_gz("out", produced.out_bytes)
    if "ps" in formats:
        write_gz("ps", produced.ps_bytes)
    if "xml" in formats:
        write_gz("xml", produced.xml_bytes)
    if "wrl" in formats:
        write_gz("wrl", produced.wrl_bytes)

    return FreezeResult(input=input_path, job_id=job_id, status="ok", input_sha256=sha, outputs=out_paths)


def _cmd_freeze(args: argparse.Namespace) -> int:
    repo = _repo_root()
    golden_dir = Path(args.golden_dir).resolve() if args.golden_dir else (repo / "test" / "golden_na")

    inputs = _collect_inputs(list(args.inputs), list(args.list_file or []))
    if not inputs:
        raise RuntimeError("no inputs selected")

    semantics = str(args.semantics).strip()
    if semantics not in ("science-v1",):
        raise RuntimeError("--semantics must be science-v1 for Gate NA")

    requested_formats = [s.strip() for s in str(args.formats).split(",") if s.strip()]
    supported = {"core", "out", "ps", "xml", "wrl"}
    for f in requested_formats:
        if f not in supported:
            raise RuntimeError(f"unsupported format in --formats: {f} (supported: {sorted(supported)})")
    formats = set(requested_formats)

    rust_bin = Path(args.rust_bin).resolve() if args.rust_bin else _ensure_rust_hotcore_binary(repo)

    started = time.time()
    results: list[FreezeResult] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(args.workers)) as ex:
        futs = [
            ex.submit(
                _freeze_one,
                repo=repo,
                rust_bin=rust_bin,
                input_path=inp,
                semantics=semantics,
                golden_dir=golden_dir,
                formats=formats,
            )
            for inp in inputs
        ]
        for fut in concurrent.futures.as_completed(futs):
            results.append(fut.result())

    results.sort(key=lambda r: (r.status, _relpath(repo, r.input)))

    ok = sum(1 for r in results if r.status == "ok")
    failed = sum(1 for r in results if r.status != "ok")

    manifest = {
        "schema_version": 1,
        "semantics": semantics,
        "formats": sorted(formats),
        "counts": {"ok": ok, "failed": failed},
        "elapsed_ms": int((time.time() - started) * 1000),
        "entries": [
            {
                "job_id": r.job_id,
                "input": _relpath(repo, r.input),
                "input_sha256": r.input_sha256,
                "status": r.status,
                "error": r.error,
                "outputs": (
                    {k: _relpath(repo, v) for k, v in (r.outputs or {}).items()} if r.outputs is not None else None
                ),
            }
            for r in results
        ],
    }
    golden_dir.mkdir(parents=True, exist_ok=True)
    (golden_dir / "manifest.json").write_text(_json_dumps(manifest, indent=2), encoding="utf-8")

    if failed:
        for r in results:
            if r.status != "ok":
                sys.stderr.write(f"FAILED: {r.input}: {r.error}\n")
        return 1

    sys.stderr.write(f"frozen: ok={ok} golden_dir={golden_dir}\n")
    return 0


@dataclass(frozen=True)
class CompareResult:
    input: str
    job_id: str
    status: str  # ok|changed|failed
    error: str | None
    diffs: dict[str, list[str]] | None


def _load_manifest(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise RuntimeError("invalid manifest JSON")
    return data


def _cmd_compare(args: argparse.Namespace) -> int:
    repo = _repo_root()
    golden_dir = Path(args.golden_dir).resolve() if args.golden_dir else (repo / "test" / "golden_na")
    manifest_path = golden_dir / "manifest.json"
    if not manifest_path.exists():
        raise RuntimeError(f"missing golden manifest: {manifest_path}")

    manifest = _load_manifest(manifest_path)
    semantics = str(args.semantics).strip() if args.semantics else str(manifest.get("semantics") or "science-v1")
    if semantics not in ("science-v1",):
        raise RuntimeError("--semantics must be science-v1 for Gate NA")

    formats = set(manifest.get("formats") or [])
    if not formats:
        raise RuntimeError("manifest is missing formats")

    entries = [e for e in (manifest.get("entries") or []) if isinstance(e, dict)]
    if not entries:
        raise RuntimeError("manifest has no entries")

    wanted: set[str] = set()
    if args.inputs or args.list_file:
        selected = _collect_inputs(list(args.inputs), list(args.list_file or []))
        wanted = {_relpath(repo, p) for p in selected}
    if wanted:
        entries = [e for e in entries if str(e.get("input") or "") in wanted]
    if not entries:
        raise RuntimeError("no entries selected for compare")

    out_dir = Path(args.out_dir).resolve() if args.out_dir else (repo / "out_gate_na")
    if out_dir.exists():
        shutil.rmtree(out_dir)
    (out_dir / "cases").mkdir(parents=True, exist_ok=True)

    rust_bin = Path(args.rust_bin).resolve() if args.rust_bin else _ensure_rust_hotcore_binary(repo)

    started = time.time()
    results: list[CompareResult] = []
    counts = {"ok": 0, "changed": 0, "failed": 0}

    for entry in entries:
        input_rel = str(entry.get("input") or "")
        job_id = str(entry.get("job_id") or _job_id_for_input(Path(input_rel) if input_rel else Path("job")))
        case_dir = out_dir / "cases" / job_id
        case_dir.mkdir(parents=True, exist_ok=True)

        input_path = (repo / input_rel).resolve() if input_rel else None
        if input_path is None or not input_path.exists():
            counts["failed"] += 1
            results.append(CompareResult(input=input_rel, job_id=job_id, status="failed", error="missing input", diffs=None))
            continue

        outputs = entry.get("outputs") or {}
        if not isinstance(outputs, dict):
            counts["failed"] += 1
            results.append(
                CompareResult(input=input_rel, job_id=job_id, status="failed", error="manifest missing outputs", diffs=None)
            )
            continue

        try:
            produced = _run_hotcore_pipeline(repo=repo, rust_bin=rust_bin, input_path=input_path, semantics=semantics)
        except Exception as e:  # noqa: BLE001
            counts["failed"] += 1
            results.append(CompareResult(input=input_rel, job_id=job_id, status="failed", error=str(e), diffs=None))
            continue

        diffs_by_fmt: dict[str, list[str]] = {}
        changed = False

        if "core" in formats:
            golden_core_path = outputs.get("core")
            if not isinstance(golden_core_path, str) or not golden_core_path:
                diffs_by_fmt["core"] = ["missing golden core path in manifest"]
                changed = True
            else:
                gp = Path(golden_core_path)
                gp = (repo / gp).resolve() if not gp.is_absolute() else gp
                golden_core = json.loads(gp.read_text(encoding="utf-8", errors="replace"))
                golden_core = _canonical_core_for_regress(golden_core)
                if produced.core != golden_core:
                    changed = True
                    a = _json_dumps(golden_core, indent=2)
                    b = _json_dumps(produced.core, indent=2)
                    diffs_by_fmt["core"] = _unified_diff(a, b, fromfile="golden.core.json", tofile="candidate.core.json", max_lines=int(args.max_diff_lines))

        def compare_gz(fmt: str, raw: bytes) -> None:
            nonlocal changed
            golden_path = outputs.get(fmt)
            if not isinstance(golden_path, str) or not golden_path:
                diffs_by_fmt[fmt] = [f"missing golden {fmt} path in manifest"]
                changed = True
                return
            gp = Path(golden_path)
            gp = (repo / gp).resolve() if not gp.is_absolute() else gp
            golden_canon = _gzip_read(gp)
            cand_canon = _canonicalize_text(raw=raw, fmt=fmt)
            if cand_canon != golden_canon:
                changed = True
                diffs_by_fmt[fmt] = _unified_diff(
                    golden_canon.decode("utf-8", errors="replace"),
                    cand_canon.decode("utf-8", errors="replace"),
                    fromfile=f"golden.{fmt}",
                    tofile=f"candidate.{fmt}",
                    max_lines=int(args.max_diff_lines),
                )
                (case_dir / f"candidate.{fmt}").write_bytes(raw)
                (case_dir / f"candidate.{fmt}.canon").write_bytes(cand_canon)
                (case_dir / f"golden.{fmt}.canon").write_bytes(golden_canon)

        if "out" in formats:
            compare_gz("out", produced.out_bytes)
        if "ps" in formats:
            compare_gz("ps", produced.ps_bytes)
        if "xml" in formats:
            compare_gz("xml", produced.xml_bytes)
        if "wrl" in formats:
            compare_gz("wrl", produced.wrl_bytes)

        if changed:
            counts["changed"] += 1
            results.append(CompareResult(input=input_rel, job_id=job_id, status="changed", error=None, diffs=diffs_by_fmt))
        else:
            counts["ok"] += 1
            results.append(CompareResult(input=input_rel, job_id=job_id, status="ok", error=None, diffs=None))

    summary = {
        "schema_version": 1,
        "semantics": semantics,
        "counts": counts,
        "elapsed_ms": int((time.time() - started) * 1000),
        "out_dir": str(out_dir),
        "results": [
            {
                "input": r.input,
                "job_id": r.job_id,
                "status": r.status,
                "error": r.error,
                "diffs": r.diffs,
            }
            for r in results
        ],
    }
    (out_dir / "summary.json").write_text(_json_dumps(summary, indent=2), encoding="utf-8")

    # Markdown report
    lines: list[str] = []
    lines.append("# Gate NA report (DNA/RNA hybrid, science-v1 self-oracle)")
    lines.append("")
    lines.append(f"- semantics: `{semantics}`")
    lines.append(f"- golden: `{_relpath(repo, golden_dir)}`")
    lines.append(f"- out_dir: `{_relpath(repo, out_dir)}`")
    lines.append(f"- counts: `{counts}`")
    lines.append("")

    for r in results:
        if r.status == "ok":
            continue
        lines.append(f"## {r.status}: {r.input}")
        lines.append("")
        lines.append(f"- job_id: `{r.job_id}`")
        if r.error:
            lines.append(f"- error: `{r.error}`")
        if r.diffs:
            for fmt, diff_lines in r.diffs.items():
                lines.append("")
                lines.append(f"### diff: {fmt}")
                lines.append("")
                lines.append("```diff")
                lines.extend(diff_lines)
                lines.append("```")
        lines.append("")

    (out_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    if counts["failed"] or counts["changed"]:
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Gate NA: DNA/RNA hybrid regression (science-v1 self-oracle)")
    sub = p.add_subparsers(dest="cmd", required=True)

    pf = sub.add_parser("freeze", help="Freeze science-v1 goldens for NA gate (core + out + render)")
    pf.add_argument("inputs", nargs="+", help="Input files/dirs/globs")
    pf.add_argument("--list-file", action="append", default=None, help="Read extra inputs from a list file")
    pf.add_argument("--golden-dir", default=None, help="Golden output dir (default: test/golden_na)")
    pf.add_argument("--semantics", default="science-v1", help="Semantics baseline (must be science-v1)")
    pf.add_argument("--formats", default="core,out,xml,ps,wrl", help="Comma-separated: core,out,xml,ps,wrl")
    pf.add_argument("--workers", type=int, default=2, help="Parallel workers (default: 2)")
    pf.add_argument("--rust-bin", default=None, help="Path to rnaview-hotcore (default: auto-build/find)")
    pf.set_defaults(func=_cmd_freeze)

    pc = sub.add_parser("compare", help="Compare current outputs against frozen goldens")
    pc.add_argument("--out-dir", default=None, help="Where to write the report (default: out_gate_na)")
    pc.add_argument("--golden-dir", default=None, help="Golden dir (default: test/golden_na)")
    pc.add_argument("--semantics", default=None, help="Semantics baseline override (default: from manifest)")
    pc.add_argument("--max-diff-lines", type=int, default=200, help="Max diff lines per format (default: 200)")
    pc.add_argument("--rust-bin", default=None, help="Path to rnaview-hotcore (default: auto-build/find)")
    pc.add_argument("inputs", nargs="*", help="Optional filter: only compare these inputs")
    pc.add_argument("--list-file", action="append", default=None, help="Optional filter: list file of inputs to compare")
    pc.set_defaults(func=_cmd_compare)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
