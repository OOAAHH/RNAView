#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
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


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _rel(repo: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(repo.resolve()).as_posix()
    except Exception:  # noqa: BLE001
        return str(path)


def _core_out_path(golden_root: Path, input_path: Path) -> Path:
    repo = _repo_root()
    test_root = repo / "test"
    try:
        rel = input_path.resolve().relative_to(test_root.resolve())
    except Exception as e:  # noqa: BLE001
        raise ValueError(f"input must be under test/: {input_path}") from e
    return golden_root / rel.with_suffix(rel.suffix + ".core.json")


def _canonical_core_for_regress(core: Any) -> Any:
    # Keep science-v1 golden cores focused on scientific content.
    # `out_index` is renderer/writer bookkeeping and must NOT participate in core equivalence
    # (see doc/spec.md).
    if not isinstance(core, dict):
        return core
    base_pairs = core.get("base_pairs")
    if isinstance(base_pairs, list):
        for bp in base_pairs:
            if isinstance(bp, dict):
                bp.pop("out_index", None)
    return core


@dataclass(frozen=True)
class FreezeResult:
    input: Path
    status: str  # ok|failed
    input_sha256: str | None
    core_json: Path | None
    reused_legacy_core: bool
    error: str | None = None


def _freeze_one(
    *,
    repo: Path,
    batch_mod: Any,
    out_core_mod: Any,
    rust_bin: Path,
    semantics: str,
    golden_root: Path,
    legacy_regress_index: Any | None,
    write_all: bool,
    input_path: Path,
) -> FreezeResult:
    try:
        sha = _sha256(input_path)
    except Exception:  # noqa: BLE001
        sha = None

    with tempfile.TemporaryDirectory() as td:
        tmp_out = Path(td)
        jr = batch_mod._run_one_rust(
            input_path=input_path,
            out_dir=tmp_out,
            job_id_mode="stem-hash",
            overwrite=True,
            rust_bin=rust_bin,
            mmcif_parser="legacy",
            rust_oracle="compute",
            semantics=semantics,
            hydrogen_policy=None,
            missing_insertion_code=None,
            chain_id_policy=None,
            legacy_bin=None,
            out_core_mod=out_core_mod,
            regress_index=None,
            regress_mode="core",
            max_diffs=50,
            keep_going=True,
        )
        if jr.status != "ok" or not jr.pairs_json:
            return FreezeResult(
                input=input_path,
                status="failed",
                input_sha256=sha,
                core_json=None,
                reused_legacy_core=False,
                error=jr.error or "engine failed",
            )

        pairs = json.loads(Path(jr.pairs_json).read_text(encoding="utf-8"))
        core = _canonical_core_for_regress(pairs.get("core", {}))

    expected_core_path: Path | None = None
    if legacy_regress_index is not None:
        golden = batch_mod._lookup_golden_entry(input_path, legacy_regress_index)
        if golden is not None:
            expected_core_path = golden.core_path

    reused = False
    out_path: Path
    if not write_all and expected_core_path is not None:
        golden_core = _canonical_core_for_regress(json.loads(expected_core_path.read_text(encoding="utf-8")))
        if core == golden_core:
            reused = True
            out_path = expected_core_path
        else:
            out_path = _core_out_path(golden_root, input_path)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(_json_dumps(core, indent=None), encoding="utf-8")
    else:
        out_path = _core_out_path(golden_root, input_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(_json_dumps(core, indent=None), encoding="utf-8")

    return FreezeResult(
        input=input_path,
        status="ok",
        input_sha256=sha,
        core_json=out_path,
        reused_legacy_core=reused,
    )


def _cmd_freeze(args: argparse.Namespace) -> int:
    repo = _repo_root()
    batch_mod = _load_module("rnaview_batch", repo / "tools" / "rnaview_batch.py")
    out_core_mod = _load_module("rnaview_out_core", repo / "tools" / "rnaview_out_core.py")

    inputs = batch_mod._collect_inputs(args.inputs, args.list)
    if not inputs:
        sys.stderr.write("no inputs\n")
        return 2

    semantics = str(args.semantics)
    if semantics not in ("science-v1",):
        sys.stderr.write("--semantics must be science-v1\n")
        return 2

    golden_root = Path(args.out_dir).resolve()
    golden_root.mkdir(parents=True, exist_ok=True)

    rust_bin = Path(args.rust_bin).resolve() if args.rust_bin else batch_mod._ensure_rust_hotcore_binary(repo)

    legacy_regress_index = None
    legacy_manifest_path = Path(args.legacy_manifest).resolve() if args.legacy_manifest else (repo / "test" / "golden_core" / "manifest.json")
    if legacy_manifest_path.exists():
        legacy_regress_index = batch_mod._build_regress_index(legacy_manifest_path)

    started = time.time()
    results: list[FreezeResult] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(args.workers)) as ex:
        futs = [
            ex.submit(
                _freeze_one,
                repo=repo,
                batch_mod=batch_mod,
                out_core_mod=out_core_mod,
                rust_bin=rust_bin,
                semantics=semantics,
                golden_root=golden_root,
                legacy_regress_index=legacy_regress_index,
                write_all=bool(args.write_all),
                input_path=inp,
            )
            for inp in inputs
        ]
        for fut in concurrent.futures.as_completed(futs):
            results.append(fut.result())

    results.sort(key=lambda r: (r.status, _rel(repo, r.input)))

    ok = sum(1 for r in results if r.status == "ok")
    failed = sum(1 for r in results if r.status != "ok")
    reused = sum(1 for r in results if r.reused_legacy_core)
    written = sum(1 for r in results if r.status == "ok" and not r.reused_legacy_core)

    manifest = {
        "schema_version": 1,
        "semantics": semantics,
        "notes": {
            "core_json may point to test/golden_core when identical": True,
            "write_all": bool(args.write_all),
        },
        "counts": {"ok": ok, "failed": failed, "reused_legacy_core": reused, "written": written},
        "elapsed_ms": int((time.time() - started) * 1000),
        "entries": [
            {
                "input": _rel(repo, r.input),
                "input_sha256": r.input_sha256,
                "core_json": _rel(repo, r.core_json) if r.core_json else None,
                "reused_legacy_core": bool(r.reused_legacy_core),
                "status": r.status,
                "error": r.error,
            }
            for r in results
        ],
    }

    (golden_root / "manifest.json").write_text(_json_dumps(manifest, indent=2), encoding="utf-8")

    if failed:
        for r in results:
            if r.status != "ok":
                sys.stderr.write(f"FAILED: {r.input}: {r.error}\n")
        return 1

    sys.stderr.write(f"frozen: ok={ok} written={written} reused={reused} out_dir={golden_root}\n")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Freeze science-v1 golden cores for regression")
    sub = p.add_subparsers(dest="cmd", required=True)

    freeze = sub.add_parser("freeze", help="Run science-v1 and write test/golden_science_core/manifest.json")
    freeze.add_argument("inputs", nargs="*", default=["test/pdb", "test/mmcif"], help="Input file/dir/glob (default: test/pdb test/mmcif)")
    freeze.add_argument("--list", action="append", default=[], help="A file with one input path per line (repeatable)")
    freeze.add_argument("--out-dir", default="test/golden_science_core", help="Output directory (default: test/golden_science_core)")
    freeze.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2), help="Parallel workers")
    freeze.add_argument("--semantics", default="science-v1", help="Semantics preset to freeze (default: science-v1)")
    freeze.add_argument("--rust-bin", default=None, help="Path to rnaview-hotcore (default: auto-build)")
    freeze.add_argument("--legacy-manifest", default=None, help="Legacy golden_core manifest (default: test/golden_core/manifest.json)")
    freeze.add_argument(
        "--write-all",
        action="store_true",
        help="Write all science cores under golden_science_core (default: only write diffs, reuse legacy core_json when identical)",
    )
    freeze.set_defaults(func=_cmd_freeze)

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
