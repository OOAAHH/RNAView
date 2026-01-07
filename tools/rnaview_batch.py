#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


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


def _job_id_for_input(path: Path, mode: str) -> str:
    if mode == "stem":
        base = path.stem
        return _sanitize_job_id(base)
    if mode == "name":
        return _sanitize_job_id(path.name)
    if mode == "stem-hash":
        base = _sanitize_job_id(path.stem)
        return f"{base}__{_hash8(str(path))}"
    if mode == "name-hash":
        base = _sanitize_job_id(path.name)
        return f"{base}__{_hash8(str(path))}"
    raise ValueError(f"unknown job-id-mode: {mode}")


def _infer_format(path: Path) -> str | None:
    ext = path.suffix.lower()
    if ext == ".cif":
        return "cif"
    if ext in (".pdb", ".ent"):
        return "pdb"
    if ext == ".out":
        return "out"
    if ext in (".xml", ".rnaml"):
        return "xml"
    return None


def _legacy_cli_flags(fmt: str | None, ps: bool) -> list[str]:
    flags: list[str] = []
    if ps:
        flags.append("-p")
    if fmt == "pdb":
        flags.append("--pdb")
    elif fmt == "cif":
        flags.append("--cif")
    return flags


@dataclass(frozen=True)
class JobResult:
    input: str
    job_id: str
    engine: str  # "legacy" | "rust"
    status: str  # "ok" | "skipped" | "failed"
    job_dir: str
    pairs_json: str | None = None
    legacy_out: str | None = None
    error: str | None = None
    regress_ok: bool | None = None
    regress_diffs: list[str] | None = None
    elapsed_ms: int | None = None


@dataclass(frozen=True)
class GoldenEntry:
    out_path: Path
    core_path: Path


@dataclass(frozen=True)
class RegressIndex:
    by_exact_input: dict[Path, GoldenEntry]
    by_dir_canon: dict[tuple[Path, str], GoldenEntry]
    by_dir_single: dict[Path, GoldenEntry]


def _maybe_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        if dst.exists():
            dst.unlink()
        os.link(src, dst)
        return
    except OSError:
        pass
    shutil.copy2(src, dst)


def _canon_stem(text: str) -> str:
    return "".join(ch.lower() for ch in text if ch.isalnum())


def _lookup_golden_entry(input_path: Path, idx: RegressIndex) -> GoldenEntry | None:
    direct = idx.by_exact_input.get(input_path)
    if direct is not None:
        return direct

    directory = input_path.parent.resolve()
    canon = _canon_stem(input_path.name)
    by_canon = idx.by_dir_canon.get((directory, canon))
    if by_canon is not None:
        return by_canon

    return idx.by_dir_single.get(directory)


def _unified_diff(a: str, b: str, *, fromfile: str, tofile: str, max_lines: int) -> list[str]:
    import difflib

    diff = difflib.unified_diff(
        a.splitlines(),
        b.splitlines(),
        fromfile=fromfile,
        tofile=tofile,
        lineterm="",
    )
    out: list[str] = []
    for line in diff:
        out.append(line)
        if len(out) >= max_lines:
            break
    return out


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
    candidates = [
        repo / "rust" / "target" / "release" / "rnaview-hotcore",
        repo / "rust" / "target" / "debug" / "rnaview-hotcore",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


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
        candidates.extend(sorted(src_dir.glob("*.rs")))

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
        raise RuntimeError("missing cargo; cannot build rust engine")

    manifest = repo / "rust" / "Cargo.toml"
    use_sysroot = shutil.which("cc") is None

    if use_sysroot:
        cmd = ["bash", str(repo / "tools" / "cargo_sysroot.sh"), "build", "--manifest-path", str(manifest)]
    else:
        cmd = [cargo, "build", "--manifest-path", str(manifest)]

    proc = subprocess.run(cmd, cwd=str(repo), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"failed to build rust engine (code={proc.returncode}):\n{proc.stdout}")

    found = _find_rust_hotcore_binary(repo)
    if found is None:
        raise RuntimeError("rust engine build succeeded but binary not found under rust/target/(debug|release)")
    return found


def _run_one_legacy(
    *,
    input_path: Path,
    out_dir: Path,
    job_id_mode: str,
    ps: bool,
    overwrite: bool,
    rnaview_bin: Path,
    out_core_mod: Any,
    pairs_mod: Any,
    regress_index: RegressIndex | None,
    regress_mode: str,
    max_diffs: int,
    keep_going: bool,
) -> JobResult:
    repo = _repo_root()
    job_id = _job_id_for_input(input_path, job_id_mode)
    job_dir = (out_dir / job_id).resolve()
    job_dir.mkdir(parents=True, exist_ok=True)

    pairs_path = job_dir / "pairs.json"
    legacy_out_path = job_dir / "legacy.out"

    if not overwrite and pairs_path.exists() and legacy_out_path.exists():
        return JobResult(
            input=str(input_path),
            job_id=job_id,
            engine="legacy",
            status="skipped",
            job_dir=str(job_dir),
            pairs_json=str(pairs_path),
            legacy_out=str(legacy_out_path),
        )

    started = time.time()
    if not rnaview_bin.exists():
        return JobResult(
            input=str(input_path),
            job_id=job_id,
            engine="legacy",
            status="failed",
            job_dir=str(job_dir),
            error=f"missing legacy binary: {rnaview_bin} (build via tools/build_legacy_rnaview.sh)",
            elapsed_ms=int((time.time() - started) * 1000),
        )

    fmt = _infer_format(input_path)
    flags = _legacy_cli_flags(fmt, ps)

    input_arg: str
    local_input: Path
    try:
        rel = input_path.resolve().relative_to(repo.resolve())
    except ValueError:
        rel = None

    if rel is not None:
        input_arg = rel.as_posix()
        local_input = job_dir / rel
    else:
        input_arg = input_path.name
        local_input = job_dir / input_arg
    try:
        _maybe_copy(input_path, local_input)
    except Exception as e:  # noqa: BLE001
        return JobResult(
            input=str(input_path),
            job_id=job_id,
            engine="legacy",
            status="failed",
            job_dir=str(job_dir),
            error=f"failed to prepare input copy: {e}",
            elapsed_ms=int((time.time() - started) * 1000),
        )

    log_path = job_dir / "legacy.log"
    env = dict(os.environ)
    env["RNAVIEW"] = str(repo)
    env["PATH"] = f"{repo / 'bin'}:{env.get('PATH', '')}"

    cmd = [str(rnaview_bin), *flags, input_arg]
    try:
        with log_path.open("wb") as log:
            proc = subprocess.run(
                cmd,
                cwd=str(job_dir),
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if proc.returncode != 0:
            return JobResult(
                input=str(input_path),
                job_id=job_id,
                engine="legacy",
                status="failed",
                job_dir=str(job_dir),
                error=f"legacy rnaview failed (code={proc.returncode}); see {log_path}",
                elapsed_ms=int((time.time() - started) * 1000),
            )

        produced_out = job_dir / Path(f"{input_arg}.out")
        if not produced_out.exists():
            candidates = sorted(job_dir.rglob("*.out"))
            if candidates:
                produced_out = candidates[0]
            else:
                return JobResult(
                    input=str(input_path),
                    job_id=job_id,
                    engine="legacy",
                    status="failed",
                    job_dir=str(job_dir),
                    error=f"legacy rnaview produced no .out; see {log_path}",
                    elapsed_ms=int((time.time() - started) * 1000),
                )

        shutil.copy2(produced_out, legacy_out_path)

        core = out_core_mod.extract_core(legacy_out_path)
        pairs = pairs_mod.pairs_json_from_core(
            core,
            source_path=str(input_path),
            source_format=fmt or "unknown",
            options={"engine": "legacy"},
        )
        pairs_path.write_text(_json_dumps(pairs, indent=None), encoding="utf-8")

        regress_ok: bool | None = None
        regress_diffs: list[str] | None = None
        if regress_index is not None:
            golden = _lookup_golden_entry(input_path, regress_index)
            if golden is None:
                regress_ok = None
            elif regress_mode == "core":
                golden_core = json.loads(golden.core_path.read_text(encoding="utf-8"))
                if core == golden_core:
                    regress_ok = True
                else:
                    diffs = list(out_core_mod._iter_differences(golden_core, core, path=""))
                    regress_ok = False
                    regress_diffs = diffs[:max_diffs]
                    if not keep_going:
                        return JobResult(
                            input=str(input_path),
                            job_id=job_id,
                            engine="legacy",
                            status="failed",
                            job_dir=str(job_dir),
                            pairs_json=str(pairs_path),
                            legacy_out=str(legacy_out_path),
                            error="core regression mismatch",
                            regress_ok=False,
                            regress_diffs=regress_diffs,
                            elapsed_ms=int((time.time() - started) * 1000),
                        )
            elif regress_mode == "out":
                golden_text = golden.out_path.read_text(encoding="utf-8", errors="replace")
                candidate_text = legacy_out_path.read_text(encoding="utf-8", errors="replace")
                if candidate_text == golden_text:
                    regress_ok = True
                else:
                    regress_ok = False
                    regress_diffs = _unified_diff(
                        golden_text,
                        candidate_text,
                        fromfile=str(golden.out_path),
                        tofile=str(legacy_out_path),
                        max_lines=max_diffs,
                    )
                    if not keep_going:
                        return JobResult(
                            input=str(input_path),
                            job_id=job_id,
                            engine="legacy",
                            status="failed",
                            job_dir=str(job_dir),
                            pairs_json=str(pairs_path),
                            legacy_out=str(legacy_out_path),
                            error=".out regression mismatch",
                            regress_ok=False,
                            regress_diffs=regress_diffs,
                            elapsed_ms=int((time.time() - started) * 1000),
                        )
            else:
                raise ValueError(f"unknown regress_mode: {regress_mode}")

        return JobResult(
            input=str(input_path),
            job_id=job_id,
            engine="legacy",
            status="ok",
            job_dir=str(job_dir),
            pairs_json=str(pairs_path),
            legacy_out=str(legacy_out_path),
            regress_ok=regress_ok,
            regress_diffs=regress_diffs,
            elapsed_ms=int((time.time() - started) * 1000),
        )
    except Exception as e:  # noqa: BLE001
        return JobResult(
            input=str(input_path),
            job_id=job_id,
            engine="legacy",
            status="failed",
            job_dir=str(job_dir),
            error=str(e),
            elapsed_ms=int((time.time() - started) * 1000),
        )


def _run_one_rust(
    *,
    input_path: Path,
    out_dir: Path,
    job_id_mode: str,
    overwrite: bool,
    rust_bin: Path,
    legacy_bin: Path,
    mmcif_parser: str,
    out_core_mod: Any,
    regress_index: RegressIndex | None,
    regress_mode: str,
    max_diffs: int,
    keep_going: bool,
) -> JobResult:
    repo = _repo_root()
    job_id = _job_id_for_input(input_path, job_id_mode)
    job_dir = (out_dir / job_id).resolve()
    job_dir.mkdir(parents=True, exist_ok=True)

    pairs_path = job_dir / "pairs.json"
    if not overwrite and pairs_path.exists():
        return JobResult(
            input=str(input_path),
            job_id=job_id,
            engine="rust",
            status="skipped",
            job_dir=str(job_dir),
            pairs_json=str(pairs_path),
        )

    started = time.time()
    fmt = _infer_format(input_path)
    if fmt not in ("pdb", "cif"):
        return JobResult(
            input=str(input_path),
            job_id=job_id,
            engine="rust",
            status="failed",
            job_dir=str(job_dir),
            error=f"unsupported input format: {input_path}",
            elapsed_ms=int((time.time() - started) * 1000),
        )

    log_path = job_dir / "rust.log"
    env = dict(os.environ)
    env["RNAVIEW"] = str(repo)
    env.update(_sysroot_env())

    cmd = [
        str(rust_bin),
        "from-structure",
        str(input_path),
        "--format",
        fmt,
        "--mmcif-parser",
        str(mmcif_parser),
        "--legacy-bin",
        str(legacy_bin),
        "--rnaview-root",
        str(repo),
        "-o",
        str(pairs_path),
    ]

    try:
        with log_path.open("wb") as log:
            proc = subprocess.run(
                cmd,
                cwd=str(job_dir),
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if proc.returncode != 0:
            return JobResult(
                input=str(input_path),
                job_id=job_id,
                engine="rust",
                status="failed",
                job_dir=str(job_dir),
                error=f"rust engine failed (code={proc.returncode}); see {log_path}",
                elapsed_ms=int((time.time() - started) * 1000),
            )
        if not pairs_path.exists():
            return JobResult(
                input=str(input_path),
                job_id=job_id,
                engine="rust",
                status="failed",
                job_dir=str(job_dir),
                error=f"rust engine produced no pairs.json; see {log_path}",
                elapsed_ms=int((time.time() - started) * 1000),
            )

        pairs = json.loads(pairs_path.read_text(encoding="utf-8"))
        core = pairs.get("core", {})

        regress_ok: bool | None = None
        regress_diffs: list[str] | None = None
        if regress_index is not None:
            if regress_mode != "core":
                raise RuntimeError("--engine rust only supports --regress-mode core (for now)")

            golden = _lookup_golden_entry(input_path, regress_index)
            if golden is None:
                regress_ok = None
            else:
                golden_core = json.loads(golden.core_path.read_text(encoding="utf-8"))
                if core == golden_core:
                    regress_ok = True
                else:
                    diffs = list(out_core_mod._iter_differences(golden_core, core, path=""))
                    regress_ok = False
                    regress_diffs = diffs[:max_diffs]
                    if not keep_going:
                        return JobResult(
                            input=str(input_path),
                            job_id=job_id,
                            engine="rust",
                            status="failed",
                            job_dir=str(job_dir),
                            pairs_json=str(pairs_path),
                            error="core regression mismatch",
                            regress_ok=False,
                            regress_diffs=regress_diffs,
                            elapsed_ms=int((time.time() - started) * 1000),
                        )

        return JobResult(
            input=str(input_path),
            job_id=job_id,
            engine="rust",
            status="ok",
            job_dir=str(job_dir),
            pairs_json=str(pairs_path),
            regress_ok=regress_ok,
            regress_diffs=regress_diffs,
            elapsed_ms=int((time.time() - started) * 1000),
        )
    except Exception as e:  # noqa: BLE001
        return JobResult(
            input=str(input_path),
            job_id=job_id,
            engine="rust",
            status="failed",
            job_dir=str(job_dir),
            error=str(e),
            elapsed_ms=int((time.time() - started) * 1000),
        )


def _build_regress_index(manifest_path: Path) -> RegressIndex:
    repo = _repo_root()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    by_exact_input: dict[Path, GoldenEntry] = {}
    by_dir: dict[Path, list[tuple[str, GoldenEntry]]] = {}
    for entry in manifest.get("entries", []):
        out_rel = Path(entry["out"])
        input_rel = out_rel.with_suffix("")
        input_path = (repo / input_rel).resolve()
        out_path = (repo / out_rel).resolve()
        core_path = (repo / entry["core_json"]).resolve()
        golden = GoldenEntry(out_path=out_path, core_path=core_path)
        by_exact_input[input_path] = golden

        directory = (repo / out_rel).resolve().parent
        canon = _canon_stem(input_rel.name)
        by_dir.setdefault(directory, []).append((canon, golden))

    by_dir_canon: dict[tuple[Path, str], GoldenEntry] = {}
    by_dir_single: dict[Path, GoldenEntry] = {}
    for directory, entries in by_dir.items():
        if len(entries) == 1:
            by_dir_single[directory] = entries[0][1]
            continue
        for canon, core_path in entries:
            by_dir_canon[(directory, canon)] = core_path

    return RegressIndex(by_exact_input=by_exact_input, by_dir_canon=by_dir_canon, by_dir_single=by_dir_single)


def _cmd_run(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    inputs = _collect_inputs(args.inputs, args.list)
    if not inputs:
        sys.stderr.write("no inputs\n")
        return 2

    engine = str(getattr(args, "engine", "legacy"))
    if engine not in ("legacy", "rust"):
        sys.stderr.write(f"invalid engine: {engine}\n")
        return 2
    if getattr(args, "ps", False) and engine != "legacy":
        sys.stderr.write("--ps is only supported with --engine legacy\n")
        return 2

    repo = _repo_root()
    rnaview_bin = Path(args.rnaview_bin).resolve() if args.rnaview_bin else (repo / "bin" / "rnaview")
    if not rnaview_bin.exists():
        sys.stderr.write(f"missing rnaview binary: {rnaview_bin} (build via tools/build_legacy_rnaview.sh)\n")
        return 2

    regress_index: RegressIndex | None = None
    if args.regress:
        manifest = Path(args.manifest) if args.manifest else _repo_root() / "test" / "golden_core" / "manifest.json"
        if not manifest.exists():
            sys.stderr.write(f"missing manifest: {manifest}\n")
            return 2
        regress_index = _build_regress_index(manifest)
        if str(getattr(args, "regress_mode", "core")) != "core" and engine == "rust":
            sys.stderr.write("--engine rust only supports --regress-mode core (for now)\n")
            return 2

    out_core_mod = _load_module("rnaview_out_core", repo / "tools" / "rnaview_out_core.py")
    pairs_mod = _load_module("rnaview_pairs_json", repo / "tools" / "rnaview_pairs_json.py")

    rust_bin: Path | None = None
    if engine == "rust":
        try:
            rust_bin = _ensure_rust_hotcore_binary(repo)
        except Exception as e:  # noqa: BLE001
            sys.stderr.write(str(e) + "\n")
            return 3

    started = time.time()

    results: list[JobResult] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(args.workers)) as ex:
        if engine == "legacy":
            futs = [
                ex.submit(
                    _run_one_legacy,
                    input_path=inp,
                    out_dir=out_dir,
                    job_id_mode=args.job_id_mode,
                    ps=bool(args.ps),
                    overwrite=bool(args.overwrite),
                    rnaview_bin=rnaview_bin,
                    out_core_mod=out_core_mod,
                    pairs_mod=pairs_mod,
                    regress_index=regress_index,
                    regress_mode=str(getattr(args, "regress_mode", "core")),
                    max_diffs=int(args.max_diffs),
                    keep_going=bool(args.keep_going),
                )
                for inp in inputs
            ]
        else:
            assert rust_bin is not None
            futs = [
                ex.submit(
                    _run_one_rust,
                    input_path=inp,
                    out_dir=out_dir,
                    job_id_mode=args.job_id_mode,
                    overwrite=bool(args.overwrite),
                    rust_bin=rust_bin,
                    legacy_bin=rnaview_bin,
                    mmcif_parser=str(getattr(args, "mmcif_parser", "legacy")),
                    out_core_mod=out_core_mod,
                    regress_index=regress_index,
                    regress_mode=str(getattr(args, "regress_mode", "core")),
                    max_diffs=int(args.max_diffs),
                    keep_going=bool(args.keep_going),
                )
                for inp in inputs
            ]
        for fut in concurrent.futures.as_completed(futs):
            results.append(fut.result())

    results.sort(key=lambda r: (r.status, r.job_id, r.input))

    ok = sum(1 for r in results if r.status == "ok")
    skipped = sum(1 for r in results if r.status == "skipped")
    failed = sum(1 for r in results if r.status == "failed")
    regress_failed = sum(1 for r in results if r.regress_ok is False)

    summary = {
        "schema_version": 1,
        "counts": {"ok": ok, "skipped": skipped, "failed": failed, "regress_failed": regress_failed},
        "elapsed_ms": int((time.time() - started) * 1000),
        "results": [
            {
                "input": r.input,
                "job_id": r.job_id,
                "engine": r.engine,
                "status": r.status,
                "job_dir": r.job_dir,
                "pairs_json": r.pairs_json,
                "legacy_out": r.legacy_out,
                "error": r.error,
                "regress_ok": r.regress_ok,
                "regress_diffs": r.regress_diffs,
                "elapsed_ms": r.elapsed_ms,
            }
            for r in results
        ],
    }
    (out_dir / "summary.json").write_text(_json_dumps(summary, indent=2), encoding="utf-8")

    if failed or regress_failed:
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Batch runner (legacy/rust engines) for RNAVIEW modernization")
    sub = p.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Run an engine and write <out_dir>/<job_id>/pairs.json (and legacy.out for legacy)")
    run.add_argument("inputs", nargs="+", help="Input file/dir/glob (repeatable)")
    run.add_argument("--list", action="append", default=[], help="A file with one input path per line (repeatable)")
    run.add_argument("--out-dir", required=True, help="Output directory")
    run.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2), help="Parallel workers")
    run.add_argument("--engine", choices=["legacy", "rust"], default="legacy", help="Which engine to run")
    run.add_argument("--rnaview-bin", default=None, help="Path to rnaview binary (default: bin/rnaview)")
    run.add_argument(
        "--mmcif-parser",
        choices=["legacy", "pdbtbx"],
        default="legacy",
        help="For --engine rust: use legacy RNAVIEW mmCIF parsing, or convert via pdbtbx then run legacy PDB parser",
    )
    run.add_argument(
        "--job-id-mode",
        choices=["stem-hash", "name-hash", "stem", "name"],
        default="stem-hash",
        help="How to derive job_id from input path",
    )
    run.add_argument("--ps", action="store_true", help="Enable legacy -p (produce PS/XML); legacy engine only")
    run.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    run.add_argument("--regress", action="store_true", help="Compare against frozen golden_core manifest")
    run.add_argument(
        "--regress-mode",
        choices=["core", "out"],
        default="core",
        help="For --regress: compare extracted core JSON, or compare full .out text (legacy engine only)",
    )
    run.add_argument("--manifest", default=None, help="Path to golden_core/manifest.json (default: test/golden_core/manifest.json)")
    run.add_argument("--max-diffs", type=int, default=50, help="Max diff lines to record on mismatch")
    run.add_argument("--keep-going", action="store_true", help="Keep running even if a mismatch is found")
    run.set_defaults(func=_cmd_run)

    args = p.parse_args(argv)
    try:
        return int(args.func(args))
    except KeyboardInterrupt:
        return 130
    except Exception as e:  # noqa: BLE001
        sys.stderr.write(f"internal error: {e}\n")
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
