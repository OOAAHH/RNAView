#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import statistics
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


def _collect_inputs(items: list[str]) -> list[Path]:
    repo = _repo_root()
    out: list[Path] = []
    allowed_exts = {".pdb", ".ent", ".cif"}
    excluded_name_suffixes = ("_tmp.pdb",)

    for item in items:
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
                if not cand.is_file():
                    continue
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


def _infer_format(path: Path) -> str | None:
    ext = path.suffix.lower()
    if ext == ".cif":
        return "cif"
    if ext in (".pdb", ".ent"):
        return "pdb"
    return None


def _engine_flags(fmt: str | None, ps: bool, *, label: bool) -> list[str]:
    flags: list[str] = []
    if ps:
        flags.append("-p")
    if fmt == "pdb":
        flags.append("--pdb")
    elif fmt == "cif":
        flags.append("--cif")
        if label:
            flags.append("--label")
    return flags


_SUITE_PHASE2_FALLBACK = [
    "test/pdb/pdb1nvy/pdb1nvy.pdb",
    "test/pdb/test1/test1.pdb",
    "test/pdb/tr0001/tr0001.pdb",
    "test/pdb/url064/url064.pdb",
    "test/pdb/urx053/urx053.pdb",
    "test/mmcif/insertion_code/1EFW/1EFW.cif",
    "test/mmcif/insertion_code/1VVJ/1VVJ.cif",
    "test/mmcif/insertion_code/4ARC/4ARC.cif",
    "test/mmcif/nmr_structure/8if5/8if5.cif",
    "test/mmcif/other/6pom/6pom.cif",
    "test/mmcif/x-ray/3P4J/assembly-1/3p4j-assembly1.cif",
    "test/mmcif/x-ray/434D/assembly-1/434d-assembly1.cif",
    "test/mmcif/x-ray/434D/assembly-2/434d-assembly2.cif",
    "test/mmcif/x-ray/4NMG/assembly-1/4nmg-assembly1.cif",
]


def _load_suite_phase2(repo: Path) -> list[str]:
    script = repo / "test_phase2.sh"
    if not script.exists():
        return list(_SUITE_PHASE2_FALLBACK)
    text = script.read_text(encoding="utf-8", errors="replace")
    import re

    pat = re.compile(r"\btest/[A-Za-z0-9_./-]+\.(?:pdb|ent|cif)\b")
    seen: set[str] = set()
    out: list[str] = []
    for m in pat.finditer(text):
        p = m.group(0)
        if p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out or list(_SUITE_PHASE2_FALLBACK)


def _sha256_bytes(data: bytes) -> str:
    h = hashlib.sha256()
    h.update(data)
    return h.hexdigest()


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


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


@dataclass(frozen=True)
class RunSample:
    returncode: int
    wall_ms: int
    user_ms: int
    sys_ms: int
    max_rss_kb: int | None


@dataclass(frozen=True)
class EngineRun:
    name: str
    bin: str
    flags: list[str]
    samples: list[RunSample]


@dataclass(frozen=True)
class ProfileSample:
    sample: RunSample
    out_bytes: bytes
    profile: dict[str, Any]


def _summarize_ms(values: list[int]) -> dict[str, Any]:
    if not values:
        return {}
    sorted_vals = sorted(values)
    out: dict[str, Any] = {
        "n": len(values),
        "min_ms": int(sorted_vals[0]),
        "max_ms": int(sorted_vals[-1]),
        "median_ms": int(statistics.median(sorted_vals)),
        "mean_ms": int(statistics.mean(sorted_vals)),
    }
    if len(values) >= 2:
        out["stdev_ms"] = float(statistics.stdev(sorted_vals))
    return out


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _relative_to_repo_or_name(repo: Path, path: Path) -> Path:
    try:
        return path.resolve().relative_to(repo.resolve())
    except ValueError:
        return Path(path.name)


def _expected_out_path(work_dir: Path, input_arg: Path) -> Path:
    return work_dir / Path(f"{input_arg.as_posix()}.out")


def _run_one_via_wrapper(
    *,
    cmd: list[str],
    cwd: Path,
    env: dict[str, str],
    timeout_s: float | None,
    log_path: Path | None,
) -> RunSample:
    spec = {
        "cmd": cmd,
        "cwd": str(cwd),
        "env": env,
        "timeout_s": timeout_s,
        "log_path": str(log_path) if log_path is not None else None,
    }
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "_run_one"],
        input=_json_dumps(spec, indent=None),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"internal runner failed (code={proc.returncode}):\n{proc.stderr}")
    try:
        data = json.loads(proc.stdout)
    except json.JSONDecodeError as e:
        raise RuntimeError(f"internal runner produced invalid json: {e}\nstdout={proc.stdout}\nstderr={proc.stderr}") from e

    return RunSample(
        returncode=int(data.get("returncode", -1)),
        wall_ms=int(data.get("wall_ms", -1)),
        user_ms=int(data.get("user_ms", -1)),
        sys_ms=int(data.get("sys_ms", -1)),
        max_rss_kb=(int(data["max_rss_kb"]) if data.get("max_rss_kb") is not None else None),
    )


def _run_engine_case(
    *,
    engine_name: str,
    engine_bin: Path,
    input_abs: Path,
    input_arg: Path,
    flags: list[str],
    base_env: dict[str, str],
    sysroot: dict[str, str],
    warmup: int,
    runs: int,
    timeout_s: float | None,
    log_dir: Path | None,
) -> tuple[list[RunSample], bytes]:
    if warmup < 0 or runs <= 0:
        raise ValueError("warmup must be >=0 and runs must be >=1")

    repo = _repo_root()
    all_samples: list[RunSample] = []

    verify_out_bytes: bytes | None = None

    total_invocations = 1 + warmup + runs  # 1 verify run (excluded) + warmup + measured runs
    for idx in range(total_invocations):
        with tempfile.TemporaryDirectory(prefix=f"rnaview-bench-{engine_name}-") as td:
            work_dir = Path(td)
            local_input = work_dir / input_arg
            _ensure_parent(local_input)
            shutil.copy2(input_abs, local_input)

            cmd = [str(engine_bin), *flags, input_arg.as_posix()]

            env = dict(base_env)
            env["RNAVIEW"] = str(repo)
            env["PATH"] = f"{repo / 'bin'}:{env.get('PATH', '')}"
            env.update(sysroot)

            log_path: Path | None = None
            if log_dir is not None:
                log_path = log_dir / f"{engine_name}.{idx:03d}.log"
                _ensure_parent(log_path)

            sample = _run_one_via_wrapper(cmd=cmd, cwd=work_dir, env=env, timeout_s=timeout_s, log_path=log_path)

            produced_out = _expected_out_path(work_dir, input_arg)
            if not produced_out.exists():
                raise RuntimeError(f"{engine_name} produced no .out at {produced_out}")

            if sample.returncode != 0:
                raise RuntimeError(f"{engine_name} exited non-zero ({sample.returncode}); see {log_path}" if log_path else "")

            if idx == 0:
                verify_out_bytes = produced_out.read_bytes()
                continue
            if idx <= warmup:
                continue
            all_samples.append(sample)

    assert verify_out_bytes is not None
    return all_samples, verify_out_bytes


def _run_engine_profile(
    *,
    engine_name: str,
    engine_bin: Path,
    input_abs: Path,
    input_arg: Path,
    flags: list[str],
    base_env: dict[str, str],
    sysroot: dict[str, str],
    timeout_s: float | None,
    log_path: Path | None,
) -> ProfileSample:
    repo = _repo_root()
    with tempfile.TemporaryDirectory(prefix=f"rnaview-profile-{engine_name}-") as td:
        work_dir = Path(td)
        local_input = work_dir / input_arg
        _ensure_parent(local_input)
        shutil.copy2(input_abs, local_input)

        cmd = [str(engine_bin), *flags, input_arg.as_posix()]

        env = dict(base_env)
        env["RNAVIEW"] = str(repo)
        env["PATH"] = f"{repo / 'bin'}:{env.get('PATH', '')}"
        env.update(sysroot)

        profile_path = work_dir / "profile.json"
        env["RNAVIEW_PROFILE_JSON"] = str(profile_path)

        sample = _run_one_via_wrapper(cmd=cmd, cwd=work_dir, env=env, timeout_s=timeout_s, log_path=log_path)

        produced_out = _expected_out_path(work_dir, input_arg)
        if not produced_out.exists():
            raise RuntimeError(f"{engine_name} produced no .out at {produced_out}")
        if sample.returncode != 0:
            raise RuntimeError(f"{engine_name} exited non-zero ({sample.returncode}); see {log_path}" if log_path else "")

        if not profile_path.exists():
            raise RuntimeError(f"{engine_name} produced no profile json at {profile_path}")
        try:
            prof = json.loads(profile_path.read_text(encoding="utf-8", errors="replace"))
        except json.JSONDecodeError as e:  # noqa: BLE001
            raise RuntimeError(f"{engine_name} produced invalid profile json: {e}") from e

        return ProfileSample(sample=sample, out_bytes=produced_out.read_bytes(), profile=prof)


def _ms(ns: int | float | None) -> float | None:
    if ns is None:
        return None
    return float(ns) / 1_000_000.0


def _fmt_ms(ns: int | float | None) -> str:
    v = _ms(ns)
    if v is None:
        return "n/a"
    return f"{v:.3f}ms"


def _cmd_profile(args: argparse.Namespace) -> int:
    repo = _repo_root()
    input_abs = Path(args.input).resolve() if args.input else None
    if input_abs is None or not input_abs.exists():
        sys.stderr.write("missing input\n")
        return 2

    fmt = _infer_format(input_abs)
    if fmt is None:
        sys.stderr.write(f"unsupported input format: {input_abs}\n")
        return 2
    input_arg = _relative_to_repo_or_name(repo, input_abs)
    flags = _engine_flags(fmt, bool(args.ps), label=bool(args.label))

    legacy_bin = Path(args.legacy_bin).resolve() if args.legacy_bin else (repo / "bin" / "rnaview")
    rustcore_bin = Path(args.rustcore_bin).resolve() if args.rustcore_bin else (repo / "bin" / "rnaview_rustcore")
    if not legacy_bin.exists():
        sys.stderr.write(f"missing legacy binary: {legacy_bin} (build via tools/build_legacy_rnaview.sh)\n")
        return 2
    if not rustcore_bin.exists():
        sys.stderr.write(f"missing rustcore binary: {rustcore_bin} (build via tools/build_rnaview_rustcore.sh)\n")
        return 2

    timeout_s = float(args.timeout) if args.timeout is not None else None
    verify = str(args.verify)
    keep_logs = bool(args.keep_logs)
    sysroot = _sysroot_env()
    base_env = dict(os.environ)

    log_root: Path | None = None
    if keep_logs:
        log_root = Path(args.log_dir).resolve() if args.log_dir else (repo / "out" / "bench_logs")
        log_root.mkdir(parents=True, exist_ok=True)

    case_log_dir = log_root / _sha256_bytes(input_arg.as_posix().encode("utf-8"))[:12] if log_root else None
    if case_log_dir is not None:
        case_log_dir.mkdir(parents=True, exist_ok=True)

    legacy_log = (case_log_dir / "legacy.profile.log") if case_log_dir is not None else None
    rust_log = (case_log_dir / "rustcore.profile.log") if case_log_dir is not None else None

    legacy = _run_engine_profile(
        engine_name="legacy",
        engine_bin=legacy_bin,
        input_abs=input_abs,
        input_arg=input_arg,
        flags=flags,
        base_env=base_env,
        sysroot={},
        timeout_s=timeout_s,
        log_path=legacy_log,
    )
    rustcore = _run_engine_profile(
        engine_name="rustcore",
        engine_bin=rustcore_bin,
        input_abs=input_abs,
        input_arg=input_arg,
        flags=flags,
        base_env=base_env,
        sysroot=sysroot,
        timeout_s=timeout_s,
        log_path=rust_log,
    )

    verify_ok: bool | None = None
    verify_diffs: list[str] | None = None
    if verify == "none":
        verify_ok = None
    elif verify == "between":
        verify_ok = legacy.out_bytes == rustcore.out_bytes
        if not verify_ok:
            verify_diffs = _unified_diff(
                legacy.out_bytes.decode("utf-8", errors="replace"),
                rustcore.out_bytes.decode("utf-8", errors="replace"),
                fromfile="legacy.out",
                tofile="rustcore.out",
                max_lines=int(args.max_diffs),
            )
    elif verify == "golden":
        golden_path = (repo / input_arg).with_suffix(f"{input_arg.suffix}.out")
        if not golden_path.exists():
            verify_ok = None
        else:
            golden_bytes = golden_path.read_bytes()
            ok_legacy = legacy.out_bytes == golden_bytes
            ok_rustcore = rustcore.out_bytes == golden_bytes
            verify_ok = ok_legacy and ok_rustcore
            if not verify_ok:
                verify_diffs = _unified_diff(
                    golden_bytes.decode("utf-8", errors="replace"),
                    legacy.out_bytes.decode("utf-8", errors="replace"),
                    fromfile="golden.out",
                    tofile="legacy.out",
                    max_lines=int(args.max_diffs),
                )
    else:
        raise ValueError(f"unknown verify mode: {verify}")

    sys.stderr.write(f"{input_arg.as_posix()}: legacy={legacy.sample.wall_ms}ms rustcore={rustcore.sample.wall_ms}ms\n")
    if verify_ok is False and verify_diffs:
        sys.stderr.write("\n".join(verify_diffs) + "\n")

    legacy_t = legacy.profile.get("times_ns") or {}
    rust_t = rustcore.profile.get("times_ns") or {}

    stages = [
        "base_info",
        "all_pairs_total",
        "all_pairs_candidate",
        "all_pairs_check_pairs",
        "all_pairs_base_stack",
        "all_pairs_hbond_pair",
        "all_pairs_lw_pair_type",
        "best_pair_total",
        "best_pair_check_pairs",
    ]

    sys.stderr.write("stage timings (legacy vs rustcore)\n")
    for k in stages:
        l = legacy_t.get(k)
        r = rust_t.get(k)
        speed = (float(l) / float(r)) if isinstance(l, (int, float)) and isinstance(r, (int, float)) and r else None
        speed_s = f"{speed:.2f}x" if speed is not None else "n/a"
        sys.stderr.write(f"  {k}: legacy={_fmt_ms(l)} rustcore={_fmt_ms(r)} speedup={speed_s}\n")

    sub_stages = [
        "all_pairs_hbond_pair_h_catalog",
        "all_pairs_lw_get_hbond_ij",
    ]
    if any((k in legacy_t) or (k in rust_t) for k in sub_stages):
        sys.stderr.write("substage timings (legacy vs rustcore)\n")
        for k in sub_stages:
            l = legacy_t.get(k)
            r = rust_t.get(k)
            speed = (float(l) / float(r)) if isinstance(l, (int, float)) and isinstance(r, (int, float)) and r else None
            speed_s = f"{speed:.2f}x" if speed is not None else "n/a"
            sys.stderr.write(f"  {k}: legacy={_fmt_ms(l)} rustcore={_fmt_ms(r)} speedup={speed_s}\n")

    def _slowest_stage(times: dict[str, Any]) -> tuple[str, int] | None:
        best: tuple[str, int] | None = None
        for k in stages:
            v = times.get(k)
            if not isinstance(v, int):
                continue
            if best is None or v > best[1]:
                best = (k, v)
        return best

    slow_legacy = _slowest_stage(legacy_t)
    slow_rust = _slowest_stage(rust_t)
    if slow_legacy is not None:
        sys.stderr.write(f"slowest stage (legacy): {slow_legacy[0]} {_fmt_ms(slow_legacy[1])}\n")
    if slow_rust is not None:
        sys.stderr.write(f"slowest stage (rustcore): {slow_rust[0]} {_fmt_ms(slow_rust[1])}\n")

    if args.json:
        out = {
            "schema_version": 1,
            "input": str(input_abs),
            "input_arg": input_arg.as_posix(),
            "flags": flags,
            "verify": {"mode": verify, "ok": verify_ok, "diffs": verify_diffs},
            "legacy": {
                "bin": str(legacy_bin),
                "sample": legacy.sample.__dict__,
                "profile": legacy.profile,
                "out_sha256": _sha256_bytes(legacy.out_bytes),
            },
            "rustcore": {
                "bin": str(rustcore_bin),
                "sample": rustcore.sample.__dict__,
                "profile": rustcore.profile,
                "out_sha256": _sha256_bytes(rustcore.out_bytes),
            },
        }
        Path(args.json).write_text(_json_dumps(out, indent=2), encoding="utf-8")

    if verify_ok is False:
        return 1
    return 0

def _cmd_compare(args: argparse.Namespace) -> int:
    repo = _repo_root()

    if args.suite == "phase2":
        default_inputs = _load_suite_phase2(repo)
    else:
        default_inputs = []

    inputs_raw = list(args.inputs or []) or default_inputs
    inputs = _collect_inputs(inputs_raw)
    if not inputs:
        sys.stderr.write("no inputs\n")
        return 2

    legacy_bin = Path(args.legacy_bin).resolve() if args.legacy_bin else (repo / "bin" / "rnaview")
    rustcore_bin = Path(args.rustcore_bin).resolve() if args.rustcore_bin else (repo / "bin" / "rnaview_rustcore")
    if not legacy_bin.exists():
        sys.stderr.write(f"missing legacy binary: {legacy_bin} (build via tools/build_legacy_rnaview.sh)\n")
        return 2
    if not rustcore_bin.exists():
        sys.stderr.write(f"missing rustcore binary: {rustcore_bin} (build via tools/build_rnaview_rustcore.sh)\n")
        return 2

    warmup = int(args.warmup)
    runs = int(args.runs)
    timeout_s = float(args.timeout) if args.timeout is not None else None
    verify = str(args.verify)
    keep_logs = bool(args.keep_logs)

    sysroot = _sysroot_env()
    base_env = dict(os.environ)

    out: dict[str, Any] = {
        "schema_version": 1,
        "meta": {
            "repo": str(repo),
            "argv": sys.argv[1:],
            "python": sys.version.splitlines()[0],
            "platform": sys.platform,
        },
        "config": {
            "warmup": warmup,
            "runs": runs,
            "timeout_s": timeout_s,
            "verify": verify,
        },
        "engines": {
            "legacy": {"bin": str(legacy_bin)},
            "rustcore": {"bin": str(rustcore_bin)},
        },
        "cases": [],
    }

    log_root: Path | None = None
    if keep_logs:
        log_root = Path(args.log_dir).resolve() if args.log_dir else (repo / "out" / "bench_logs")
        log_root.mkdir(parents=True, exist_ok=True)

    failures = 0
    for input_abs in inputs:
        fmt = _infer_format(input_abs)
        if fmt is None:
            continue
        input_arg = _relative_to_repo_or_name(repo, input_abs)
        flags = _engine_flags(fmt, bool(args.ps), label=bool(args.label))

        case_log_dir = (log_root / _sha256_bytes(input_arg.as_posix().encode("utf-8"))[:12]) if log_root else None
        if case_log_dir is not None:
            case_log_dir.mkdir(parents=True, exist_ok=True)

        started = time.time()
        try:
            legacy_samples, legacy_out = _run_engine_case(
                engine_name="legacy",
                engine_bin=legacy_bin,
                input_abs=input_abs,
                input_arg=input_arg,
                flags=flags,
                base_env=base_env,
                sysroot={},
                warmup=warmup,
                runs=runs,
                timeout_s=timeout_s,
                log_dir=case_log_dir if keep_logs else None,
            )
            rustcore_samples, rustcore_out = _run_engine_case(
                engine_name="rustcore",
                engine_bin=rustcore_bin,
                input_abs=input_abs,
                input_arg=input_arg,
                flags=flags,
                base_env=base_env,
                sysroot=sysroot,
                warmup=warmup,
                runs=runs,
                timeout_s=timeout_s,
                log_dir=case_log_dir if keep_logs else None,
            )
        except Exception as e:  # noqa: BLE001
            failures += 1
            out["cases"].append(
                {
                    "input": str(input_abs),
                    "input_arg": input_arg.as_posix(),
                    "format": fmt,
                    "error": str(e),
                    "elapsed_ms": int((time.time() - started) * 1000),
                }
            )
            if not args.keep_going:
                break
            continue

        verify_ok: bool | None = None
        verify_diffs: list[str] | None = None
        if verify == "none":
            verify_ok = None
        elif verify == "between":
            verify_ok = legacy_out == rustcore_out
            if not verify_ok:
                verify_diffs = _unified_diff(
                    legacy_out.decode("utf-8", errors="replace"),
                    rustcore_out.decode("utf-8", errors="replace"),
                    fromfile="legacy.out",
                    tofile="rustcore.out",
                    max_lines=int(args.max_diffs),
                )
        elif verify == "golden":
            golden_path = (repo / input_arg).with_suffix(f"{input_arg.suffix}.out")
            if not golden_path.exists():
                verify_ok = None
            else:
                golden_bytes = golden_path.read_bytes()
                ok_legacy = legacy_out == golden_bytes
                ok_rustcore = rustcore_out == golden_bytes
                verify_ok = ok_legacy and ok_rustcore
                if not verify_ok:
                    verify_diffs = _unified_diff(
                        golden_bytes.decode("utf-8", errors="replace"),
                        legacy_out.decode("utf-8", errors="replace"),
                        fromfile="golden.out",
                        tofile="legacy.out",
                        max_lines=int(args.max_diffs),
                    )
        else:
            raise ValueError(f"unknown verify mode: {verify}")

        legacy_wall = [s.wall_ms for s in legacy_samples]
        rust_wall = [s.wall_ms for s in rustcore_samples]
        legacy_median = statistics.median(legacy_wall) if legacy_wall else None
        rust_median = statistics.median(rust_wall) if rust_wall else None
        speedup = (float(legacy_median) / float(rust_median)) if legacy_median and rust_median else None

        case = {
            "input": str(input_abs),
            "input_arg": input_arg.as_posix(),
            "format": fmt,
            "flags": flags,
            "verify": {"mode": verify, "ok": verify_ok, "diffs": verify_diffs},
            "elapsed_ms": int((time.time() - started) * 1000),
            "runs": {
                "legacy": {
                    "bin": str(legacy_bin),
                    "samples": [s.__dict__ for s in legacy_samples],
                    "summary": {
                        "wall_ms": _summarize_ms([s.wall_ms for s in legacy_samples]),
                        "user_ms": _summarize_ms([s.user_ms for s in legacy_samples]),
                        "sys_ms": _summarize_ms([s.sys_ms for s in legacy_samples]),
                        "max_rss_kb": _summarize_ms([s.max_rss_kb for s in legacy_samples if s.max_rss_kb is not None]),
                    },
                    "out_sha256": _sha256_bytes(legacy_out),
                    "out_bytes": len(legacy_out),
                },
                "rustcore": {
                    "bin": str(rustcore_bin),
                    "samples": [s.__dict__ for s in rustcore_samples],
                    "summary": {
                        "wall_ms": _summarize_ms([s.wall_ms for s in rustcore_samples]),
                        "user_ms": _summarize_ms([s.user_ms for s in rustcore_samples]),
                        "sys_ms": _summarize_ms([s.sys_ms for s in rustcore_samples]),
                        "max_rss_kb": _summarize_ms([s.max_rss_kb for s in rustcore_samples if s.max_rss_kb is not None]),
                    },
                    "out_sha256": _sha256_bytes(rustcore_out),
                    "out_bytes": len(rustcore_out),
                },
            },
            "comparison": {"speedup_median": speedup},
        }
        out["cases"].append(case)

        legacy_line = case["runs"]["legacy"]["summary"]["wall_ms"].get("median_ms")
        rust_line = case["runs"]["rustcore"]["summary"]["wall_ms"].get("median_ms")
        speed_line = f"{speedup:.2f}x" if speedup is not None else "n/a"
        sys.stderr.write(f"{input_arg.as_posix()}: legacy={legacy_line}ms rustcore={rust_line}ms speedup={speed_line}\n")

        if verify_ok is False:
            failures += 1
            if not args.keep_going:
                break

    out["totals"] = {"cases": len(out["cases"]), "failures": failures}

    if args.json:
        Path(args.json).write_text(_json_dumps(out, indent=2), encoding="utf-8")

    if failures:
        return 1
    return 0


def _cmd_run_one(args: argparse.Namespace) -> int:
    import resource

    spec = json.loads(sys.stdin.read() or "{}")
    cmd = list(spec.get("cmd") or [])
    cwd = spec.get("cwd")
    env = spec.get("env")
    timeout_s = spec.get("timeout_s")
    log_path = spec.get("log_path")

    if not cmd or cwd is None or env is None:
        sys.stderr.write("invalid spec\n")
        return 2

    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    log_f = None
    if log_path:
        log_f = open(log_path, "wb")  # noqa: SIM115
        stdout = log_f
        stderr = subprocess.STDOUT

    start = time.perf_counter()
    returncode = -1
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            env=dict(env),
            stdout=stdout,
            stderr=stderr,
            check=False,
            timeout=float(timeout_s) if timeout_s is not None else None,
        )
        returncode = int(proc.returncode)
    finally:
        end = time.perf_counter()
        if log_f is not None:
            log_f.close()

    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = {
        "returncode": returncode,
        "wall_ms": int((end - start) * 1000),
        "user_ms": int(ru.ru_utime * 1000),
        "sys_ms": int(ru.ru_stime * 1000),
        "max_rss_kb": int(getattr(ru, "ru_maxrss", 0)) if getattr(ru, "ru_maxrss", None) is not None else None,
    }
    sys.stdout.write(_json_dumps(result, indent=None))
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Benchmark legacy RNAVIEW vs rnaview_rustcore")
    sub = p.add_subparsers(dest="cmd", required=True)

    compare = sub.add_parser("compare", help="Compare legacy vs rustcore performance on a set of inputs")
    compare.add_argument("inputs", nargs="*", help="Input files/dirs/globs (default: --suite phase2)")
    compare.add_argument("--suite", choices=["phase2"], default="phase2", help="Default input suite if none provided")
    compare.add_argument("--legacy-bin", default=None, help="Path to legacy bin/rnaview")
    compare.add_argument("--rustcore-bin", default=None, help="Path to bin/rnaview_rustcore")
    compare.add_argument("--warmup", type=int, default=0, help="Warmup runs per case/engine (not counted)")
    compare.add_argument("--runs", type=int, default=1, help="Measured runs per case/engine")
    compare.add_argument("--timeout", type=float, default=None, help="Per-run timeout seconds")
    compare.add_argument("--ps", action="store_true", help="Enable -p (PS/XML); not recommended for benchmarking")
    compare.add_argument("--label", action="store_true", help="For mmCIF: use --label parsing (default: auth)")
    compare.add_argument(
        "--verify",
        choices=["between", "golden", "none"],
        default="between",
        help="Verify .out byte-for-byte between engines, against golden, or skip verification",
    )
    compare.add_argument("--max-diffs", type=int, default=50, help="Max unified-diff lines to keep on mismatch")
    compare.add_argument("--keep-logs", action="store_true", help="Keep stdout/stderr logs for each run")
    compare.add_argument("--log-dir", default=None, help="Directory for logs (default: repo/out/bench_logs)")
    compare.add_argument("-o", "--json", default=None, help="Write full results JSON to a file")
    compare.add_argument("--keep-going", action="store_true", help="Continue on failures/mismatches")
    compare.set_defaults(func=_cmd_compare)

    profile = sub.add_parser("profile", help="Profile stage timings for a single input")
    profile.add_argument("input", help="Input file (.pdb/.ent/.cif)")
    profile.add_argument("--legacy-bin", default=None, help="Path to legacy bin/rnaview")
    profile.add_argument("--rustcore-bin", default=None, help="Path to bin/rnaview_rustcore")
    profile.add_argument("--timeout", type=float, default=None, help="Per-run timeout seconds")
    profile.add_argument("--ps", action="store_true", help="Enable -p (PS/XML); not recommended for profiling")
    profile.add_argument("--label", action="store_true", help="For mmCIF: use --label parsing (default: auth)")
    profile.add_argument(
        "--verify",
        choices=["between", "golden", "none"],
        default="between",
        help="Verify .out byte-for-byte between engines, against golden, or skip verification",
    )
    profile.add_argument("--max-diffs", type=int, default=50, help="Max unified-diff lines to keep on mismatch")
    profile.add_argument("--keep-logs", action="store_true", help="Keep stdout/stderr logs for each run")
    profile.add_argument("--log-dir", default=None, help="Directory for logs (default: repo/out/bench_logs)")
    profile.add_argument("-o", "--json", default=None, help="Write profile JSON to a file")
    profile.set_defaults(func=_cmd_profile)

    run_one = sub.add_parser("_run_one", add_help=False)
    run_one.set_defaults(func=_cmd_run_one)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
