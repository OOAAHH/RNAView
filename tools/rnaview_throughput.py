#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _json_dumps(obj: Any, *, indent: int | None = None) -> str:
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
            for m in sorted(glob.glob(item, recursive=True)):
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


def _rel_to_repo_or_name(repo: Path, path: Path) -> Path:
    try:
        return path.resolve().relative_to(repo.resolve())
    except Exception:  # noqa: BLE001
        return Path(path.name)


def _prepare_repeated_inputs(*, repo: Path, inputs: list[Path], repeat: int, staging_dir: Path) -> list[str]:
    if repeat <= 1:
        return [str(p) for p in inputs]

    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True, exist_ok=True)

    for rep in range(repeat):
        rep_dir = staging_dir / f"rep{rep:03d}"
        for p in inputs:
            rel = _rel_to_repo_or_name(repo, p)
            dst = rep_dir / rel
            _maybe_copy(p, dst)

    return [str(staging_dir)]


def _pct(sorted_vals: list[int], pct: float) -> int | None:
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return int(sorted_vals[0])
    if pct <= 0:
        return int(sorted_vals[0])
    if pct >= 100:
        return int(sorted_vals[-1])

    k = (len(sorted_vals) - 1) * (pct / 100.0)
    f = int(k)
    c = min(f + 1, len(sorted_vals) - 1)
    if f == c:
        return int(sorted_vals[f])
    a = sorted_vals[f]
    b = sorted_vals[c]
    return int(round(a + (b - a) * (k - f)))


def _summarize_latencies(elapsed_ms: list[int]) -> dict[str, Any]:
    if not elapsed_ms:
        return {}
    vals = sorted(int(x) for x in elapsed_ms)
    out: dict[str, Any] = {
        "n": len(vals),
        "min_ms": int(vals[0]),
        "max_ms": int(vals[-1]),
        "median_ms": int(statistics.median(vals)),
        "mean_ms": float(statistics.mean(vals)),
        "p90_ms": _pct(vals, 90),
        "p95_ms": _pct(vals, 95),
        "p99_ms": _pct(vals, 99),
    }
    if len(vals) >= 2:
        out["stdev_ms"] = float(statistics.stdev(vals))
    return out


def _run_batch(
    *,
    repo: Path,
    engine: str,
    inputs: list[str],
    out_dir: Path,
    workers: int,
    extra_args: list[str],
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "batch.log"

    cmd = [
        sys.executable,
        str(repo / "tools" / "rnaview_batch.py"),
        "run",
        *inputs,
        "--out-dir",
        str(out_dir),
        "--workers",
        str(int(workers)),
        "--engine",
        str(engine),
        "--overwrite",
        "--keep-going",
    ]
    cmd.extend(extra_args)

    started = time.perf_counter()
    proc = subprocess.run(
        cmd,
        cwd=str(repo),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    ended = time.perf_counter()
    log_path.write_text(proc.stdout, encoding="utf-8", errors="replace")

    summary_path = out_dir / "summary.json"
    summary: dict[str, Any] | None = None
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8", errors="replace"))
        except Exception:  # noqa: BLE001
            summary = None

    return {
        "cmd": cmd,
        "returncode": int(proc.returncode),
        "log": str(log_path),
        "summary_path": str(summary_path),
        "summary": summary,
        "wall_ms": int((ended - started) * 1000),
    }


def _derive_metrics(summary: dict[str, Any] | None) -> dict[str, Any]:
    if not summary:
        return {}
    counts = dict(summary.get("counts") or {})
    results = list(summary.get("results") or [])

    elapsed_ms = int(summary.get("elapsed_ms") or 0)
    executed = int(counts.get("ok", 0)) + int(counts.get("failed", 0))
    elapsed_s = elapsed_ms / 1000.0 if elapsed_ms > 0 else 0.0
    throughput = (executed / elapsed_s) if elapsed_s > 0 else None

    per_job_ms: list[int] = []
    for r in results:
        if not isinstance(r, dict):
            continue
        if r.get("status") == "skipped":
            continue
        ms = r.get("elapsed_ms")
        if isinstance(ms, int):
            per_job_ms.append(ms)

    return {
        "elapsed_ms": elapsed_ms,
        "counts": counts,
        "executed_jobs": executed,
        "throughput_jobs_per_sec": throughput,
        "latency_ms": _summarize_latencies(per_job_ms),
    }


def _cmd_run(args: argparse.Namespace) -> int:
    repo = _repo_root()

    inputs: list[str]
    if args.inputs:
        inputs = args.inputs
    else:
        if args.suite == "phase2":
            inputs = _load_suite_phase2(repo)
        else:
            raise RuntimeError(f"unsupported suite: {args.suite}")

    collected = _collect_inputs(inputs)
    if not collected:
        raise RuntimeError("no inputs found")

    out_dir = Path(args.out_dir).resolve() if args.out_dir else (repo / "out" / f"throughput_{args.engine}_{int(time.time())}")
    staging_dir = Path(args.staging_dir).resolve() if args.staging_dir else (out_dir / "_inputs")

    prepared_inputs = _prepare_repeated_inputs(repo=repo, inputs=collected, repeat=int(args.repeat), staging_dir=staging_dir)

    extra: list[str] = []
    if args.legacy_bin:
        extra.extend(["--rnaview-bin", str(args.legacy_bin)])
    if args.engine == "rust" and args.rust_bin:
        extra.extend(["--rust-bin", str(args.rust_bin)])

    if args.engine == "legacy":
        if args.ps:
            extra.append("--ps")
    else:
        extra.extend(["--mmcif-parser", str(args.mmcif_parser)])
        extra.extend(["--rust-oracle", str(args.rust_oracle)])
        if args.semantics is not None:
            extra.extend(["--semantics", str(args.semantics)])
        if args.hydrogen_policy is not None:
            extra.extend(["--hydrogen-policy", str(args.hydrogen_policy)])
        if args.missing_insertion_code is not None:
            extra.extend(["--missing-insertion-code", str(args.missing_insertion_code)])
        if args.chain_id_policy is not None:
            extra.extend(["--chain-id-policy", str(args.chain_id_policy)])

    run = _run_batch(repo=repo, engine=str(args.engine), inputs=prepared_inputs, out_dir=out_dir, workers=int(args.workers), extra_args=extra)
    metrics = _derive_metrics(run.get("summary"))

    sys.stdout.write(
        _json_dumps(
            {
                "schema_version": 1,
                "kind": "throughput",
                "engine": str(args.engine),
                "workers": int(args.workers),
                "repeat": int(args.repeat),
                "input_count": len(collected),
                "prepared_inputs": prepared_inputs,
                "out_dir": str(out_dir),
                "run": run,
                "metrics": metrics,
            },
            indent=2 if args.pretty else None,
        )
    )

    # Return the underlying run code so CI can optionally treat failures as failures.
    return int(run.get("returncode", 3))


def _cmd_compare(args: argparse.Namespace) -> int:
    repo = _repo_root()

    inputs: list[str]
    if args.inputs:
        inputs = args.inputs
    else:
        if args.suite == "phase2":
            inputs = _load_suite_phase2(repo)
        else:
            raise RuntimeError(f"unsupported suite: {args.suite}")

    collected = _collect_inputs(inputs)
    if not collected:
        raise RuntimeError("no inputs found")

    out_dir = Path(args.out_dir).resolve() if args.out_dir else (repo / "out" / f"throughput_compare_{int(time.time())}")
    staging_dir = Path(args.staging_dir).resolve() if args.staging_dir else (out_dir / "_inputs")
    prepared_inputs = _prepare_repeated_inputs(repo=repo, inputs=collected, repeat=int(args.repeat), staging_dir=staging_dir)

    legacy_extra: list[str] = []
    if args.legacy_bin:
        legacy_extra.extend(["--rnaview-bin", str(args.legacy_bin)])
    if args.ps:
        legacy_extra.append("--ps")

    rust_extra: list[str] = []
    if args.rust_bin:
        rust_extra.extend(["--rust-bin", str(args.rust_bin)])
    if args.legacy_bin:
        rust_extra.extend(["--rnaview-bin", str(args.legacy_bin)])
    rust_extra.extend(["--mmcif-parser", str(args.mmcif_parser)])
    rust_extra.extend(["--rust-oracle", str(args.rust_oracle)])
    if args.semantics is not None:
        rust_extra.extend(["--semantics", str(args.semantics)])
    if args.hydrogen_policy is not None:
        rust_extra.extend(["--hydrogen-policy", str(args.hydrogen_policy)])
    if args.missing_insertion_code is not None:
        rust_extra.extend(["--missing-insertion-code", str(args.missing_insertion_code)])
    if args.chain_id_policy is not None:
        rust_extra.extend(["--chain-id-policy", str(args.chain_id_policy)])

    legacy_dir = out_dir / "legacy"
    rust_dir = out_dir / "rust"

    legacy_run = _run_batch(
        repo=repo,
        engine="legacy",
        inputs=prepared_inputs,
        out_dir=legacy_dir,
        workers=int(args.workers),
        extra_args=legacy_extra,
    )
    rust_run = _run_batch(
        repo=repo,
        engine="rust",
        inputs=prepared_inputs,
        out_dir=rust_dir,
        workers=int(args.workers),
        extra_args=rust_extra,
    )

    legacy_metrics = _derive_metrics(legacy_run.get("summary"))
    rust_metrics = _derive_metrics(rust_run.get("summary"))

    speedup: float | None = None
    l = legacy_metrics.get("throughput_jobs_per_sec")
    r = rust_metrics.get("throughput_jobs_per_sec")
    if isinstance(l, (int, float)) and isinstance(r, (int, float)) and l > 0:
        speedup = float(r) / float(l)

    sys.stdout.write(
        _json_dumps(
            {
                "schema_version": 1,
                "kind": "throughput_compare",
                "workers": int(args.workers),
                "repeat": int(args.repeat),
                "input_count": len(collected),
                "prepared_inputs": prepared_inputs,
                "out_dir": str(out_dir),
                "legacy": {"run": legacy_run, "metrics": legacy_metrics},
                "rust": {"run": rust_run, "metrics": rust_metrics},
                "speedup_rust_over_legacy": speedup,
            },
            indent=2 if args.pretty else None,
        )
    )

    # Fail only if both succeeded? Keep it simple: propagate worst returncode.
    return int(max(int(legacy_run.get("returncode", 3)), int(rust_run.get("returncode", 3))))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="High-throughput throughput benchmark (wraps tools/rnaview_batch.py)")
    sub = p.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Run one engine and report throughput/latency stats")
    run.add_argument("inputs", nargs="*", help="Input files/dirs/globs (default: --suite phase2)")
    run.add_argument("--suite", choices=["phase2"], default="phase2", help="Default input suite if no inputs provided")
    run.add_argument("--engine", choices=["legacy", "rust"], default="rust", help="Engine to benchmark")
    run.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2), help="Parallel workers")
    run.add_argument("--repeat", type=int, default=1, help="Repeat the suite by staging copies (simulates many inputs)")
    run.add_argument("--staging-dir", default=None, help="Directory for staged inputs (default: <out_dir>/_inputs)")
    run.add_argument("--out-dir", default=None, help="Output directory (default: repo/out/throughput_<engine>_<ts>)")
    run.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")

    run.add_argument("--legacy-bin", default=None, help="Override legacy bin/rnaview (used by --engine legacy; also used by --engine rust --rust-oracle legacy)")
    run.add_argument("--rust-bin", default=None, help="Override rust engine binary rnaview-hotcore (used by --engine rust)")
    run.add_argument("--ps", action="store_true", help="Legacy only: enable -p (PS/XML); not recommended for perf")

    run.add_argument("--mmcif-parser", choices=["legacy", "pdbtbx"], default="legacy", help="Rust only: mmCIF parsing mode")
    run.add_argument("--rust-oracle", choices=["legacy", "out", "compute"], default="compute", help="Rust only: oracle mode")
    run.add_argument("--semantics", choices=["legacy-v1", "science-v1"], default=None, help="Rust only: semantics preset (default: legacy-v1)")
    run.add_argument(
        "--hydrogen-policy",
        choices=["legacy-mmcif-bug", "discard-all", "keep-all"],
        default=None,
        help="Rust only: hydrogen handling policy (overrides semantics default)",
    )
    run.add_argument(
        "--missing-insertion-code",
        choices=["legacy-question-mark", "none"],
        default=None,
        help="Rust only: missing insertion-code policy (overrides semantics default)",
    )
    run.add_argument(
        "--chain-id-policy",
        choices=["legacy-1char", "unique-1char"],
        default=None,
        help="Rust only: chain-id mapping policy (overrides semantics default)",
    )
    run.set_defaults(func=_cmd_run)

    cmp = sub.add_parser("compare", help="Run legacy and rust on the same staged inputs and compare throughput")
    cmp.add_argument("inputs", nargs="*", help="Input files/dirs/globs (default: --suite phase2)")
    cmp.add_argument("--suite", choices=["phase2"], default="phase2", help="Default input suite if no inputs provided")
    cmp.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2), help="Parallel workers")
    cmp.add_argument("--repeat", type=int, default=1, help="Repeat the suite by staging copies (simulates many inputs)")
    cmp.add_argument("--staging-dir", default=None, help="Directory for staged inputs (default: <out_dir>/_inputs)")
    cmp.add_argument("--out-dir", default=None, help="Output directory base (default: repo/out/throughput_compare_<ts>)")
    cmp.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")

    cmp.add_argument("--legacy-bin", default=None, help="Override legacy bin/rnaview")
    cmp.add_argument("--rust-bin", default=None, help="Override rust engine binary (rnaview-hotcore)")
    cmp.add_argument("--ps", action="store_true", help="Legacy only: enable -p (PS/XML); not recommended for perf")

    cmp.add_argument("--mmcif-parser", choices=["legacy", "pdbtbx"], default="legacy", help="Rust: mmCIF parsing mode")
    cmp.add_argument("--rust-oracle", choices=["legacy", "out", "compute"], default="compute", help="Rust: oracle mode")
    cmp.add_argument("--semantics", choices=["legacy-v1", "science-v1"], default=None, help="Rust: semantics preset (default: legacy-v1)")
    cmp.add_argument(
        "--hydrogen-policy",
        choices=["legacy-mmcif-bug", "discard-all", "keep-all"],
        default=None,
        help="Rust: hydrogen handling policy (overrides semantics default)",
    )
    cmp.add_argument(
        "--missing-insertion-code",
        choices=["legacy-question-mark", "none"],
        default=None,
        help="Rust: missing insertion-code policy (overrides semantics default)",
    )
    cmp.add_argument(
        "--chain-id-policy",
        choices=["legacy-1char", "unique-1char"],
        default=None,
        help="Rust: chain-id mapping policy (overrides semantics default)",
    )
    cmp.set_defaults(func=_cmd_compare)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
