#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _json_dumps(obj: Any, *, indent: int | None = None) -> str:
    if indent is None:
        return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n"
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, indent=indent) + "\n"


def _relpath(repo: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(repo.resolve()).as_posix()
    except Exception:  # noqa: BLE001
        return str(path)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _load_allowlist(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {"schema_version": 1, "approved": []}
    if not path.exists():
        return {"schema_version": 1, "approved": []}
    text = path.read_text(encoding="utf-8")
    text_stripped = text.strip()
    if not text_stripped:
        return {"schema_version": 1, "approved": []}

    try:
        data = json.loads(text_stripped)
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


def _bp_key(bp: dict[str, Any]) -> tuple[int, int, str]:
    return (int(bp.get("i", -1)), int(bp.get("j", -1)), str(bp.get("kind", "")))


def _canon_bp(bp: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k in (
        "i",
        "j",
        "kind",
        "chain_i",
        "resseq_i",
        "insertion_code_i",
        "base_i",
        "base_j",
        "insertion_code_j",
        "resseq_j",
        "chain_j",
        "lw",
        "orientation",
        "syn",
        "note",
        "text",
    ):
        if k in bp:
            out[k] = bp[k]
    return out


def _multiplet_key(m: dict[str, Any]) -> tuple[tuple[int, ...], str]:
    indices = tuple(int(x) for x in (m.get("indices") or []))
    text = str(m.get("text") or "").strip()
    return (indices, text)


def _stable_id(prefix: str, input_id: str, payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    h = hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]
    return f"{input_id}|{prefix}|{h}"


def _diff_core(*, input_id: str, baseline: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    events: list[dict[str, Any]] = []

    b_stats = dict(baseline.get("stats") or {})
    c_stats = dict(candidate.get("stats") or {})
    stats_delta: dict[str, Any] = {}
    for k in ("total_pairs", "total_bases"):
        if b_stats.get(k) != c_stats.get(k):
            stats_delta[k] = {"before": b_stats.get(k), "after": c_stats.get(k)}
    b_counts = dict(b_stats.get("pair_type_counts") or {})
    c_counts = dict(c_stats.get("pair_type_counts") or {})
    if b_counts != c_counts:
        keys = sorted(set(b_counts.keys()) | set(c_counts.keys()))
        delta: dict[str, Any] = {}
        for k in keys:
            bv = b_counts.get(k)
            cv = c_counts.get(k)
            if bv != cv:
                delta[k] = {"before": bv, "after": cv}
        if delta:
            stats_delta["pair_type_counts"] = delta

    if stats_delta:
        ident = _stable_id("stats", input_id, stats_delta)
        events.append({"id": ident, "kind": "stats", "change": "changed", "delta": stats_delta})

    b_bps = { _bp_key(bp): _canon_bp(bp) for bp in (baseline.get("base_pairs") or []) }
    c_bps = { _bp_key(bp): _canon_bp(bp) for bp in (candidate.get("base_pairs") or []) }
    added_keys = sorted(set(c_bps.keys()) - set(b_bps.keys()))
    removed_keys = sorted(set(b_bps.keys()) - set(c_bps.keys()))
    common_keys = sorted(set(b_bps.keys()) & set(c_bps.keys()))

    bp_added: list[dict[str, Any]] = []
    bp_removed: list[dict[str, Any]] = []
    bp_changed: list[dict[str, Any]] = []

    for k in added_keys:
        rec = c_bps[k]
        bp_added.append(rec)
        ident = _stable_id("bp:add", input_id, {"key": k, "record": rec})
        events.append({"id": ident, "kind": "base_pair", "change": "added", "key": {"i": k[0], "j": k[1], "kind": k[2]}, "record": rec})
    for k in removed_keys:
        rec = b_bps[k]
        bp_removed.append(rec)
        ident = _stable_id("bp:rm", input_id, {"key": k, "record": rec})
        events.append({"id": ident, "kind": "base_pair", "change": "removed", "key": {"i": k[0], "j": k[1], "kind": k[2]}, "record": rec})
    for k in common_keys:
        before = b_bps[k]
        after = c_bps[k]
        if before == after:
            continue
        fields = sorted(set(before.keys()) | set(after.keys()))
        delta: dict[str, Any] = {}
        for f in fields:
            if before.get(f) != after.get(f):
                delta[f] = {"before": before.get(f), "after": after.get(f)}
        entry = {"key": {"i": k[0], "j": k[1], "kind": k[2]}, "delta": delta, "before": before, "after": after}
        bp_changed.append(entry)
        ident = _stable_id("bp:chg", input_id, {"key": k, "delta": delta})
        events.append({"id": ident, "kind": "base_pair", "change": "changed", "key": entry["key"], "delta": delta})

    b_mps = { _multiplet_key(m): dict(m) for m in (baseline.get("multiplets") or []) }
    c_mps = { _multiplet_key(m): dict(m) for m in (candidate.get("multiplets") or []) }
    mp_added_keys = sorted(set(c_mps.keys()) - set(b_mps.keys()))
    mp_removed_keys = sorted(set(b_mps.keys()) - set(c_mps.keys()))
    mp_added = [c_mps[k] for k in mp_added_keys]
    mp_removed = [b_mps[k] for k in mp_removed_keys]

    for k in mp_added_keys:
        rec = c_mps[k]
        ident = _stable_id("mp:add", input_id, {"key": k, "record": rec})
        events.append({"id": ident, "kind": "multiplet", "change": "added", "record": rec})
    for k in mp_removed_keys:
        rec = b_mps[k]
        ident = _stable_id("mp:rm", input_id, {"key": k, "record": rec})
        events.append({"id": ident, "kind": "multiplet", "change": "removed", "record": rec})

    return {
        "schema_version": 1,
        "stats_delta": stats_delta,
        "base_pairs": {
            "added": bp_added,
            "removed": bp_removed,
            "changed": bp_changed,
        },
        "multiplets": {
            "added": mp_added,
            "removed": mp_removed,
        },
        "events": events,
        "counts": {
            "stats_changed": 1 if stats_delta else 0,
            "base_pairs_added": len(bp_added),
            "base_pairs_removed": len(bp_removed),
            "base_pairs_changed": len(bp_changed),
            "multiplets_added": len(mp_added),
            "multiplets_removed": len(mp_removed),
            "events": len(events),
        },
    }


@dataclass(frozen=True)
class EngineRun:
    semantics: str
    status: str  # ok|failed
    pairs_path: Path | None = None
    out_path: Path | None = None
    log_path: Path | None = None
    elapsed_ms: int = 0
    error: str | None = None


@dataclass(frozen=True)
class CaseResult:
    input_path: Path
    input_id: str
    input_sha256: str | None
    job_id: str
    baseline: EngineRun
    candidate: EngineRun
    diff_path: Path | None
    diff_summary: dict[str, Any] | None
    unapproved_ids: list[str] | None


def _run_engine_compute(
    *,
    repo: Path,
    batch_mod: Any,
    rust_bin: Path,
    input_path: Path,
    job_dir: Path,
    semantics: str,
) -> EngineRun:
    started = time.time()

    fmt = batch_mod._infer_format(input_path)
    if fmt not in ("pdb", "cif"):
        return EngineRun(semantics=semantics, status="failed", elapsed_ms=int((time.time() - started) * 1000), error=f"unsupported input format: {input_path}")

    job_dir.mkdir(parents=True, exist_ok=True)
    pairs_path = job_dir / "pairs.json"
    out_path = job_dir / "engine.out"
    log_path = job_dir / "rust.log"

    env = dict(os.environ)
    env["RNAVIEW"] = str(repo)
    env.update(batch_mod._sysroot_env())

    cmd = [
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

    with log_path.open("wb") as log:
        proc = subprocess.run(cmd, cwd=str(job_dir), env=env, stdout=log, stderr=subprocess.STDOUT, check=False)

    if proc.returncode != 0:
        return EngineRun(
            semantics=semantics,
            status="failed",
            pairs_path=pairs_path if pairs_path.exists() else None,
            out_path=out_path if out_path.exists() else None,
            log_path=log_path,
            elapsed_ms=int((time.time() - started) * 1000),
            error=f"rust engine failed (code={proc.returncode}); see {log_path}",
        )
    if not pairs_path.exists():
        return EngineRun(
            semantics=semantics,
            status="failed",
            out_path=out_path if out_path.exists() else None,
            log_path=log_path,
            elapsed_ms=int((time.time() - started) * 1000),
            error=f"rust engine produced no pairs.json; see {log_path}",
        )
    if not out_path.exists():
        return EngineRun(
            semantics=semantics,
            status="failed",
            pairs_path=pairs_path,
            log_path=log_path,
            elapsed_ms=int((time.time() - started) * 1000),
            error=f"rust engine produced no .out; see {log_path}",
        )

    return EngineRun(
        semantics=semantics,
        status="ok",
        pairs_path=pairs_path,
        out_path=out_path,
        log_path=log_path,
        elapsed_ms=int((time.time() - started) * 1000),
    )


def _run_case(
    *,
    repo: Path,
    batch_mod: Any,
    rust_bin: Path,
    out_dir: Path,
    job_id_mode: str,
    baseline_semantics: str,
    candidate_semantics: str,
    allow_ids: set[str],
    input_path: Path,
) -> CaseResult:
    input_id = _relpath(repo, input_path)
    sha = None
    try:
        sha = _sha256(input_path)
    except Exception:  # noqa: BLE001
        sha = None

    job_id = batch_mod._job_id_for_input(input_path, job_id_mode)
    case_dir = out_dir / "cases" / job_id
    baseline_dir = case_dir / baseline_semantics
    candidate_dir = case_dir / candidate_semantics

    baseline = _run_engine_compute(
        repo=repo,
        batch_mod=batch_mod,
        rust_bin=rust_bin,
        input_path=input_path,
        job_dir=baseline_dir,
        semantics=baseline_semantics,
    )
    candidate = _run_engine_compute(
        repo=repo,
        batch_mod=batch_mod,
        rust_bin=rust_bin,
        input_path=input_path,
        job_dir=candidate_dir,
        semantics=candidate_semantics,
    )

    diff_path = case_dir / "diff.json"
    if baseline.status != "ok" or candidate.status != "ok":
        diff_obj = {
            "schema_version": 1,
            "input": input_id,
            "job_id": job_id,
            "semantics": {"baseline": baseline_semantics, "candidate": candidate_semantics},
            "status": "failed",
            "input_sha256": sha,
            "baseline": {
                "status": baseline.status,
                "pairs_json": _relpath(repo, baseline.pairs_path) if baseline.pairs_path else None,
                "out_path": _relpath(repo, baseline.out_path) if baseline.out_path else None,
                "log": _relpath(repo, baseline.log_path) if baseline.log_path else None,
                "error": baseline.error,
                "elapsed_ms": baseline.elapsed_ms,
            },
            "candidate": {
                "status": candidate.status,
                "pairs_json": _relpath(repo, candidate.pairs_path) if candidate.pairs_path else None,
                "out_path": _relpath(repo, candidate.out_path) if candidate.out_path else None,
                "log": _relpath(repo, candidate.log_path) if candidate.log_path else None,
                "error": candidate.error,
                "elapsed_ms": candidate.elapsed_ms,
            },
        }
        diff_path.write_text(_json_dumps(diff_obj, indent=2), encoding="utf-8")
        return CaseResult(
            input_path=input_path,
            input_id=input_id,
            input_sha256=sha,
            job_id=job_id,
            baseline=baseline,
            candidate=candidate,
            diff_path=diff_path,
            diff_summary=None,
            unapproved_ids=None,
        )

    baseline_pairs = json.loads(baseline.pairs_path.read_text(encoding="utf-8"))
    candidate_pairs = json.loads(candidate.pairs_path.read_text(encoding="utf-8"))
    baseline_core = baseline_pairs.get("core", {})
    candidate_core = candidate_pairs.get("core", {})
    diff = _diff_core(input_id=input_id, baseline=baseline_core, candidate=candidate_core)

    unapproved = sorted({e["id"] for e in diff.get("events", []) if isinstance(e, dict) and e.get("id") not in allow_ids})
    approved_ok = not unapproved

    diff_obj = {
        "schema_version": 1,
        "input": input_id,
        "job_id": job_id,
        "semantics": {"baseline": baseline_semantics, "candidate": candidate_semantics},
        "status": "ok" if approved_ok else "unapproved",
        "input_sha256": sha,
        "baseline": {
            "pairs_json": _relpath(repo, baseline.pairs_path) if baseline.pairs_path else None,
            "out_path": _relpath(repo, baseline.out_path) if baseline.out_path else None,
            "log": _relpath(repo, baseline.log_path) if baseline.log_path else None,
            "elapsed_ms": baseline.elapsed_ms,
            "options": baseline_pairs.get("options", {}),
        },
        "candidate": {
            "pairs_json": _relpath(repo, candidate.pairs_path) if candidate.pairs_path else None,
            "out_path": _relpath(repo, candidate.out_path) if candidate.out_path else None,
            "log": _relpath(repo, candidate.log_path) if candidate.log_path else None,
            "elapsed_ms": candidate.elapsed_ms,
            "options": candidate_pairs.get("options", {}),
        },
        "diff": diff,
        "allowlist": {
            "unapproved_ids": unapproved,
        },
    }
    diff_path.write_text(_json_dumps(diff_obj, indent=2), encoding="utf-8")

    return CaseResult(
        input_path=input_path,
        input_id=input_id,
        input_sha256=sha,
        job_id=job_id,
        baseline=baseline,
        candidate=candidate,
        diff_path=diff_path,
        diff_summary=diff.get("counts", {}),
        unapproved_ids=unapproved,
    )


def _cmd_run(args: argparse.Namespace) -> int:
    repo = _repo_root()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    batch_mod = _load_module("rnaview_batch", repo / "tools" / "rnaview_batch.py")

    inputs = batch_mod._collect_inputs(args.inputs, args.list)
    if not inputs:
        sys.stderr.write("no inputs\n")
        return 2

    baseline_semantics = str(args.baseline_semantics)
    candidate_semantics = str(args.candidate_semantics)
    if baseline_semantics not in ("legacy-v1", "science-v1") or candidate_semantics not in ("legacy-v1", "science-v1"):
        sys.stderr.write("--baseline-semantics/--candidate-semantics must be legacy-v1|science-v1\n")
        return 2

    allowlist_path: Path | None = None
    if args.allowlist is not None:
        allowlist_path = Path(args.allowlist).resolve()
    else:
        default = repo / "test" / "gate_c_allowlist.yaml"
        if default.exists():
            allowlist_path = default

    allowlist = _load_allowlist(allowlist_path)
    allow_ids = _approved_ids(allowlist)

    rust_bin = Path(args.rust_bin).resolve() if args.rust_bin else batch_mod._ensure_rust_hotcore_binary(repo)

    started = time.time()
    results: list[CaseResult] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(args.workers)) as ex:
        futs = [
            ex.submit(
                _run_case,
                repo=repo,
                batch_mod=batch_mod,
                rust_bin=rust_bin,
                out_dir=out_dir,
                job_id_mode=str(args.job_id_mode),
                baseline_semantics=baseline_semantics,
                candidate_semantics=candidate_semantics,
                allow_ids=allow_ids,
                input_path=inp,
            )
            for inp in inputs
        ]
        for fut in concurrent.futures.as_completed(futs):
            results.append(fut.result())

    def _status(r: CaseResult) -> str:
        if r.baseline.status != "ok" or r.candidate.status != "ok":
            return "failed"
        if r.unapproved_ids:
            return "unapproved"
        if r.diff_summary and int(r.diff_summary.get("events", 0)) > 0:
            return "changed"
        return "ok"

    results.sort(key=lambda r: (_status(r), r.input_id))

    counts = {"ok": 0, "changed": 0, "unapproved": 0, "failed": 0}
    for r in results:
        counts[_status(r)] += 1

    summary = {
        "schema_version": 1,
        "semantics": {"baseline": baseline_semantics, "candidate": candidate_semantics},
        "allowlist": {
            "path": _relpath(repo, allowlist_path) if allowlist_path else None,
            "schema_version": allowlist.get("schema_version", 1),
            "approved_count": len(allow_ids),
        },
        "counts": counts,
        "elapsed_ms": int((time.time() - started) * 1000),
        "results": [
            {
                "input": r.input_id,
                "input_sha256": r.input_sha256,
                "job_id": r.job_id,
                "status": _status(r),
                "case_dir": _relpath(repo, (out_dir / "cases" / r.job_id)),
                "diff_json": _relpath(repo, r.diff_path) if r.diff_path else None,
                "diff_counts": r.diff_summary,
                "unapproved_ids": r.unapproved_ids,
                "baseline": {
                    "semantics": baseline_semantics,
                    "status": r.baseline.status,
                    "pairs_json": _relpath(repo, r.baseline.pairs_path) if r.baseline.pairs_path else None,
                    "out_path": _relpath(repo, r.baseline.out_path) if r.baseline.out_path else None,
                    "log": _relpath(repo, r.baseline.log_path) if r.baseline.log_path else None,
                    "elapsed_ms": r.baseline.elapsed_ms,
                    "error": r.baseline.error,
                },
                "candidate": {
                    "semantics": candidate_semantics,
                    "status": r.candidate.status,
                    "pairs_json": _relpath(repo, r.candidate.pairs_path) if r.candidate.pairs_path else None,
                    "out_path": _relpath(repo, r.candidate.out_path) if r.candidate.out_path else None,
                    "log": _relpath(repo, r.candidate.log_path) if r.candidate.log_path else None,
                    "elapsed_ms": r.candidate.elapsed_ms,
                    "error": r.candidate.error,
                },
            }
            for r in results
        ],
    }
    (out_dir / "summary.json").write_text(_json_dumps(summary, indent=2), encoding="utf-8")

    report_lines: list[str] = []
    report_lines.append("# Gate C report\n")
    report_lines.append(f"- baseline: `{baseline_semantics}`\n")
    report_lines.append(f"- candidate: `{candidate_semantics}`\n")
    report_lines.append(f"- allowlist: `{summary['allowlist']['path']}`\n")
    report_lines.append(f"- counts: `{counts}`\n")
    report_lines.append("\n")

    def _render_case(r: CaseResult) -> str:
        st = _status(r)
        diff_counts = r.diff_summary or {}
        diff_json = _relpath(repo, r.diff_path) if r.diff_path else ""
        return f"- {st}: `{r.input_id}` job_id=`{r.job_id}` events={diff_counts.get('events', 0)} diff=`{diff_json}`\n"

    if counts["failed"]:
        report_lines.append("## Failed\n\n")
        for r in results:
            if _status(r) == "failed":
                report_lines.append(_render_case(r))
        report_lines.append("\n")
    if counts["unapproved"]:
        report_lines.append("## Unapproved\n\n")
        for r in results:
            if _status(r) == "unapproved":
                report_lines.append(_render_case(r))
        report_lines.append("\n")
    if counts["changed"]:
        report_lines.append("## Changed (approved)\n\n")
        for r in results:
            if _status(r) == "changed":
                report_lines.append(_render_case(r))
        report_lines.append("\n")

    (out_dir / "report.md").write_text("".join(report_lines), encoding="utf-8")

    sys.stderr.write(f"Gate C: ok={counts['ok']} changed={counts['changed']} unapproved={counts['unapproved']} failed={counts['failed']}\n")
    if counts["failed"] or counts["unapproved"]:
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Gate C (science mode): diff legacy-v1 vs science-v1 with allowlist")
    sub = p.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Run baseline/candidate and write diff.json + report.md + summary.json")
    run.add_argument("inputs", nargs="+", help="Input file/dir/glob (repeatable)")
    run.add_argument("--list", action="append", default=[], help="A file with one input path per line (repeatable)")
    run.add_argument("--out-dir", required=True, help="Output directory")
    run.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2), help="Parallel workers")
    run.add_argument("--job-id-mode", choices=["stem-hash", "name-hash", "stem", "name"], default="stem-hash", help="How to derive job_id from input path")
    run.add_argument("--baseline-semantics", choices=["legacy-v1", "science-v1"], default="legacy-v1", help="Baseline semantics preset")
    run.add_argument("--candidate-semantics", choices=["legacy-v1", "science-v1"], default="science-v1", help="Candidate semantics preset")
    run.add_argument("--allowlist", default=None, help="Allowlist path (default: test/gate_c_allowlist.yaml if present)")
    run.add_argument("--rust-bin", default=None, help="Path to rnaview-hotcore (default: auto-build via tools/rnaview_batch.py)")
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

