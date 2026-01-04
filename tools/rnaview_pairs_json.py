#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_rnaview_out_core():
    mod_path = _repo_root() / "tools" / "rnaview_out_core.py"
    spec = importlib.util.spec_from_file_location("rnaview_out_core", mod_path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def pairs_json_from_core(
    core: dict[str, Any],
    *,
    source_path: str | None = None,
    source_format: str | None = None,
    options: dict[str, Any] | None = None,
    schema_version: int = 1,
) -> dict[str, Any]:
    out: dict[str, Any] = {"schema_version": int(schema_version), "core": core}
    if source_path is not None or source_format is not None:
        out["source"] = {}
        if source_path is not None:
            out["source"]["path"] = str(source_path)
        if source_format is not None:
            out["source"]["format"] = str(source_format)
    if options is not None:
        out["options"] = options
    return out


def pairs_json_from_out(out_path: Path) -> dict[str, Any]:
    mod = _load_rnaview_out_core()
    core = mod.extract_core(out_path)
    return pairs_json_from_core(core, source_path=str(out_path), source_format="out", options={})


def _bp_sort_key(bp: dict[str, Any]) -> tuple[Any, ...]:
    return (
        bp.get("i", -1),
        bp.get("j", -1),
        bp.get("kind", ""),
        bp.get("chain_i", ""),
        bp.get("resseq_i", -1),
        bp.get("base_i", ""),
        bp.get("base_j", ""),
        bp.get("resseq_j", -1),
        bp.get("chain_j", ""),
        bp.get("lw", ""),
        bp.get("orientation", ""),
        bp.get("syn", 0),
        bp.get("note", ""),
        bp.get("text", ""),
    )


def _format_base_pair_line(bp: dict[str, Any]) -> str:
    kind = bp.get("kind", "")
    if kind == "unknown":
        return str(bp.get("text", "")).strip()

    i = int(bp["i"])
    j = int(bp["j"])
    chain_i = str(bp.get("chain_i") or " ")
    chain_j = str(bp.get("chain_j") or " ")
    resseq_i = int(bp["resseq_i"])
    resseq_j = int(bp["resseq_j"])
    base_i = str(bp["base_i"])
    base_j = str(bp["base_j"])

    syn_count = int(bp.get("syn", 0) or 0)
    if syn_count < 0:
        syn_count = 0

    tokens: list[str] = []
    if kind == "stacked":
        tokens.extend(["syn"] * syn_count)
        tokens.append("stacked")
    else:
        lw = bp.get("lw")
        if not lw:
            raise ValueError(f"base-pair record missing lw: {bp}")
        tokens.append(str(lw))
        orientation = bp.get("orientation")
        if orientation:
            tokens.append(str(orientation))
        tokens.extend(["syn"] * syn_count)
        note = bp.get("note")
        if note:
            tokens.append(str(note))

    rest = " ".join(tokens).strip()
    return f"{i}_{j}, {chain_i}: {resseq_i} {base_i}-{base_j} {resseq_j} {chain_j}: {rest}".rstrip()


def out_core_from_pairs_json(pairs: dict[str, Any]) -> str:
    core = pairs.get("core", pairs)
    base_pairs = list(core.get("base_pairs", []))
    multiplets = list(core.get("multiplets", []))
    stats: dict[str, Any] = dict(core.get("stats", {}))

    base_pairs.sort(key=_bp_sort_key)
    multiplets.sort(key=lambda m: (m.get("indices", []), m.get("text", "")))

    lines: list[str] = []
    lines.append("BEGIN_base-pair")
    for bp in base_pairs:
        s = _format_base_pair_line(bp)
        if s:
            lines.append(s)
    lines.append("END_base-pair")
    lines.append("")
    lines.append("Summary of triplets and higher multiplets")
    lines.append("BEGIN_multiplets")
    for m in multiplets:
        indices = list(m.get("indices", []))
        text = str(m.get("text", "")).strip()
        if indices:
            idx_part = "_".join(str(int(x)) for x in indices)
            lines.append(f"{idx_part}_| {text}".rstrip())
        else:
            if text:
                lines.append(text)
    lines.append("END_multiplets")
    lines.append("")

    total_pairs = stats.get("total_pairs")
    total_bases = stats.get("total_bases")
    if total_pairs is not None and total_bases is not None:
        lines.append(f"  The total base pairs = {int(total_pairs):3d} (from {int(total_bases):4d} bases)")
        pair_type_counts = stats.get("pair_type_counts") or {}
        if pair_type_counts:
            keys = list(pair_type_counts.keys())
            lines.append("------------------------------------------------")
            lines.append(" ".join(str(k) for k in keys))
            lines.append(" ".join(str(int(pair_type_counts[k])) for k in keys))
            lines.append("------------------------------------------------")

    return "\n".join(lines) + "\n"


def _json_dumps(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n"


def _cmd_from_out(args: argparse.Namespace) -> int:
    pairs = pairs_json_from_out(Path(args.out))
    text = _json_dumps(pairs)
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


def _cmd_write_out(args: argparse.Namespace) -> int:
    pairs = json.loads(Path(args.pairs).read_text(encoding="utf-8"))
    out_text = out_core_from_pairs_json(pairs)
    if args.output:
        Path(args.output).write_text(out_text, encoding="utf-8")
    else:
        sys.stdout.write(out_text)
    return 0


def _iter_differences(a: Any, b: Any) -> Iterable[str]:
    mod = _load_rnaview_out_core()
    yield from mod._iter_differences(a, b, path="")


def _cmd_validate_golden(args: argparse.Namespace) -> int:
    repo = _repo_root()
    manifest_path = Path(args.manifest) if args.manifest else repo / "test" / "golden_core" / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    mod = _load_rnaview_out_core()

    failed: list[tuple[str, list[str]]] = []
    for entry in manifest.get("entries", []):
        out_path = repo / entry["out"]
        core_json_path = repo / entry["core_json"]
        golden_core = json.loads(core_json_path.read_text(encoding="utf-8"))

        pairs = pairs_json_from_core(
            golden_core,
            source_path=entry["out"],
            source_format="out",
            options={},
        )
        out_text = out_core_from_pairs_json(pairs)

        with tempfile.NamedTemporaryFile("w", suffix=".out", encoding="utf-8", delete=True) as tmp:
            tmp.write(out_text)
            tmp.flush()
            candidate_core = mod.extract_core(Path(tmp.name))

        if candidate_core != golden_core:
            diffs = list(_iter_differences(golden_core, candidate_core))
            failed.append((entry["out"], diffs[: args.max_diffs]))
            if not args.keep_going:
                break

    if failed:
        for path, diffs in failed:
            sys.stderr.write(f"{path}: core mismatch\n")
            for d in diffs:
                sys.stderr.write("  " + d + "\n")
        return 1

    sys.stderr.write(f"ok entries={len(manifest.get('entries', []))}\n")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="pairs.json helpers for RNAVIEW modernization")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_from_out = sub.add_parser("from-out", help="Convert RNAVIEW .out to pairs.json (core only)")
    p_from_out.add_argument("out", help="Path to RNAVIEW .out file")
    p_from_out.add_argument("-o", "--output", help="Write pairs.json to a file (default: stdout)")
    p_from_out.set_defaults(func=_cmd_from_out)

    p_write_out = sub.add_parser("write-out", help="Render a minimal RNAVIEW .out(core) from pairs.json")
    p_write_out.add_argument("pairs", help="Path to pairs.json")
    p_write_out.add_argument("-o", "--output", help="Write .out(core) to a file (default: stdout)")
    p_write_out.set_defaults(func=_cmd_write_out)

    p_validate = sub.add_parser(
        "validate-golden",
        help="Validate pairs.json -> .out(core) writer against frozen golden cores",
    )
    p_validate.add_argument(
        "--manifest",
        default=None,
        help="Path to golden_core/manifest.json (default: test/golden_core/manifest.json)",
    )
    p_validate.add_argument("--max-diffs", type=int, default=20, help="Max diff lines per file")
    p_validate.add_argument("--keep-going", action="store_true", help="Report all mismatches")
    p_validate.set_defaults(func=_cmd_validate_golden)

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
