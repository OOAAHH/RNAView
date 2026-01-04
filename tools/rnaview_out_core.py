#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


_RE_BEGIN_BASE_PAIR = re.compile(r"^\s*BEGIN_base-pair\s*$")
_RE_END_BASE_PAIR = re.compile(r"^\s*END_base-pair\s*$")
_RE_BEGIN_MULTIPLETS = re.compile(r"^\s*BEGIN_multiplets\s*$")
_RE_END_MULTIPLETS = re.compile(r"^\s*END_multiplets\s*$")
_RE_TOTAL = re.compile(
    r"^\s*The total base pairs\s*=\s*(\d+)\s*\(from\s*(\d+)\s*bases\)\s*$",
    re.IGNORECASE,
)
_RE_SEPARATOR = re.compile(r"^\s*-{5,}\s*$")

_RE_PAIR_LINE = re.compile(
    r"^\s*(\d+)_(\d+),\s*(.):\s*([0-9]+)\s+(\S)-(\S)\s+([0-9]+)\s+(.):\s*(.*?)\s*$"
)


def _norm_ws(text: str) -> str:
    return " ".join(text.strip().split())


def _norm_orientation(token: str | None) -> str | None:
    if token is None:
        return None
    t = token.strip().lower()
    if t == "cis":
        return "cis"
    if t.startswith("tran"):
        return "tran"
    return t


@dataclass(frozen=True)
class BasePairRecord:
    i: int
    j: int
    chain_i: str
    resseq_i: int
    base_i: str
    base_j: str
    resseq_j: int
    chain_j: str
    kind: str  # "pair" | "stacked" | "unknown"
    lw: str | None = None
    orientation: str | None = None
    syn: int = 0
    note: str | None = None
    text: str | None = None  # only for kind=="unknown"

    def to_json(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "i": self.i,
            "j": self.j,
            "chain_i": self.chain_i,
            "resseq_i": self.resseq_i,
            "base_i": self.base_i,
            "base_j": self.base_j,
            "resseq_j": self.resseq_j,
            "chain_j": self.chain_j,
            "kind": self.kind,
        }
        if self.kind == "unknown":
            out["text"] = self.text or ""
            return out
        if self.lw is not None:
            out["lw"] = self.lw
        if self.orientation is not None:
            out["orientation"] = self.orientation
        if self.syn:
            out["syn"] = self.syn
        if self.note is not None:
            out["note"] = self.note
        return out


@dataclass(frozen=True)
class MultipletRecord:
    indices: tuple[int, ...]
    text: str

    def to_json(self) -> dict[str, Any]:
        return {"indices": list(self.indices), "text": self.text}


def _parse_pair_rest(rest: str) -> tuple[str, str | None, int, str | None, str]:
    text = _norm_ws(rest)
    if not text:
        return "", None, 0, None, ""

    tokens = text.split()

    if tokens and tokens[-1].lower() == "stacked":
        syn_count = sum(1 for t in tokens if t.lower() == "syn")
        return "stacked", None, syn_count, None, "stacked"

    lw = tokens[0] if tokens else ""
    orientation = None
    syn_count = 0
    note_tokens: list[str] = []

    for t in tokens[1:]:
        tl = t.lower()
        if tl in ("cis", "tran", "trans"):
            orientation = _norm_orientation(tl)
            continue
        if tl == "syn":
            syn_count += 1
            continue
        note_tokens.append(t)

    note = " ".join(note_tokens).strip() or None
    return "pair", orientation, syn_count, note, lw


def _parse_base_pair_line(line: str) -> BasePairRecord | None:
    m = _RE_PAIR_LINE.match(line)
    if not m:
        return None

    i = int(m.group(1))
    j = int(m.group(2))
    chain_i = m.group(3)
    resseq_i = int(m.group(4))
    base_i = m.group(5)
    base_j = m.group(6)
    resseq_j = int(m.group(7))
    chain_j = m.group(8)
    rest = m.group(9)

    kind, orientation, syn_count, note, lw = _parse_pair_rest(rest)
    if kind == "stacked":
        return BasePairRecord(
            i=i,
            j=j,
            chain_i=chain_i,
            resseq_i=resseq_i,
            base_i=base_i,
            base_j=base_j,
            resseq_j=resseq_j,
            chain_j=chain_j,
            kind="stacked",
            syn=syn_count,
        )

    return BasePairRecord(
        i=i,
        j=j,
        chain_i=chain_i,
        resseq_i=resseq_i,
        base_i=base_i,
        base_j=base_j,
        resseq_j=resseq_j,
        chain_j=chain_j,
        kind="pair",
        lw=lw or None,
        orientation=orientation,
        syn=syn_count,
        note=note,
    )


def _extract_block(lines: list[str], begin_re: re.Pattern[str], end_re: re.Pattern[str]) -> list[str]:
    in_block = False
    out: list[str] = []
    for line in lines:
        if not in_block:
            if begin_re.match(line):
                in_block = True
            continue
        if end_re.match(line):
            break
        out.append(line.rstrip("\n"))
    return out


def _extract_stats(lines: list[str]) -> dict[str, Any]:
    stats: dict[str, Any] = {}

    total_idx = None
    total_pairs = None
    total_bases = None
    for idx, line in enumerate(lines):
        m = _RE_TOTAL.match(line.rstrip("\n"))
        if m:
            total_pairs = int(m.group(1))
            total_bases = int(m.group(2))
            total_idx = idx
            break

    if total_pairs is None or total_bases is None or total_idx is None:
        return stats

    stats["total_pairs"] = total_pairs
    stats["total_bases"] = total_bases

    pair_type_counts: dict[str, int] = {}
    pending_header: list[str] | None = None

    for line in lines[total_idx + 1 :]:
        line = line.rstrip("\n")
        if _RE_SEPARATOR.match(line):
            continue
        tokens = line.split()
        if not tokens:
            continue
        if any("--" in t for t in tokens):
            pending_header = tokens
            continue
        if pending_header is None:
            continue
        if all(re.fullmatch(r"-?\d+", t) for t in tokens):
            if len(tokens) == len(pending_header):
                for key, val in zip(pending_header, tokens):
                    pair_type_counts[key] = int(val)
            pending_header = None

    if pair_type_counts:
        stats["pair_type_counts"] = dict(sorted(pair_type_counts.items()))
    return stats


def extract_core(out_path: Path) -> dict[str, Any]:
    lines = out_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)

    base_pair_lines = _extract_block(lines, _RE_BEGIN_BASE_PAIR, _RE_END_BASE_PAIR)
    base_pairs: list[BasePairRecord] = []
    for raw in base_pair_lines:
        if not raw.strip():
            continue
        parsed = _parse_base_pair_line(raw)
        if parsed is None:
            base_pairs.append(
                BasePairRecord(
                    i=-1,
                    j=-1,
                    chain_i="",
                    resseq_i=-1,
                    base_i="",
                    base_j="",
                    resseq_j=-1,
                    chain_j="",
                    kind="unknown",
                    text=_norm_ws(raw),
                )
            )
        else:
            base_pairs.append(parsed)

    multiplet_lines = _extract_block(lines, _RE_BEGIN_MULTIPLETS, _RE_END_MULTIPLETS)
    multiplets: list[MultipletRecord] = []
    for raw in multiplet_lines:
        s = raw.strip()
        if not s:
            continue
        if "|" in s:
            left, right = s.split("|", 1)
            left = left.rstrip("_").strip()
            idxs: list[int] = []
            for part in left.split("_"):
                part = part.strip()
                if not part:
                    continue
                if part.isdigit():
                    idxs.append(int(part))
            multiplets.append(MultipletRecord(indices=tuple(idxs), text=_norm_ws(right)))
        else:
            multiplets.append(MultipletRecord(indices=tuple(), text=_norm_ws(s)))

    stats = _extract_stats(lines)

    def bp_sort_key(bp: BasePairRecord) -> tuple[Any, ...]:
        return (
            bp.i,
            bp.j,
            bp.kind,
            bp.chain_i,
            bp.resseq_i,
            bp.base_i,
            bp.base_j,
            bp.resseq_j,
            bp.chain_j,
            bp.lw or "",
            bp.orientation or "",
            bp.syn,
            bp.note or "",
            bp.text or "",
        )

    core = {
        "base_pairs": [bp.to_json() for bp in sorted(base_pairs, key=bp_sort_key)],
        "multiplets": [m.to_json() for m in sorted(multiplets, key=lambda x: (x.indices, x.text))],
        "stats": stats,
    }
    return core


def _has_base_pair_block(out_path: Path) -> bool:
    try:
        for line in out_path.open("r", encoding="utf-8", errors="replace"):
            if _RE_BEGIN_BASE_PAIR.match(line):
                return True
    except OSError:
        return False
    return False


def _iter_differences(a: Any, b: Any, path: str = "") -> Iterable[str]:
    if a == b:
        return
    if isinstance(a, dict) and isinstance(b, dict):
        a_keys = set(a.keys())
        b_keys = set(b.keys())
        for k in sorted(a_keys - b_keys):
            yield f"{path}/{k}: only in left"
        for k in sorted(b_keys - a_keys):
            yield f"{path}/{k}: only in right"
        for k in sorted(a_keys & b_keys):
            yield from _iter_differences(a[k], b[k], f"{path}/{k}")
        return
    if isinstance(a, list) and isinstance(b, list):
        yield f"{path}: list length {len(a)} != {len(b)}"
        for i, (ai, bi) in enumerate(zip(a, b)):
            if ai != bi:
                yield from _iter_differences(ai, bi, f"{path}[{i}]")
                break
        return
    yield f"{path}: {a!r} != {b!r}"


def _cmd_extract(args: argparse.Namespace) -> int:
    core = extract_core(Path(args.out))
    sys.stdout.write(json.dumps(core, sort_keys=True, ensure_ascii=False, separators=(",", ":")))
    sys.stdout.write("\n")
    return 0


def _cmd_compare(args: argparse.Namespace) -> int:
    left = extract_core(Path(args.left))
    right = extract_core(Path(args.right))
    if left == right:
        return 0

    diffs = list(_iter_differences(left, right, path=""))
    for line in diffs[: args.max_diffs]:
        sys.stderr.write(line + "\n")
    if len(diffs) > args.max_diffs:
        sys.stderr.write(f"... and {len(diffs) - args.max_diffs} more\n")
    return 1


def _cmd_freeze(args: argparse.Namespace) -> int:
    repo_root = Path(__file__).resolve().parents[1]

    root = Path(args.root)
    if not root.is_absolute():
        root = repo_root / root
    root = root.resolve()

    out_dir = Path(args.out_dir) if args.out_dir else root / "golden_core"
    if not out_dir.is_absolute():
        out_dir = repo_root / out_dir
    out_dir = out_dir.resolve()

    exclude_suffixes = tuple(args.exclude_suffix)
    allow_unknown = bool(args.allow_unknown)
    allow_missing_stats = bool(args.allow_missing_stats)
    keep_going = bool(args.keep_going)
    dry_run = bool(args.dry_run)

    candidates: list[Path] = []
    for out_file in sorted(root.rglob("*.out")):
        if any(out_file.name.endswith(sfx) for sfx in exclude_suffixes):
            continue
        if not _has_base_pair_block(out_file):
            continue
        candidates.append(out_file)

    errors: list[str] = []
    manifest_entries: list[dict[str, Any]] = []

    if not dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    def _rel(p: Path) -> str:
        try:
            return p.relative_to(repo_root).as_posix()
        except ValueError:
            return str(p)

    for out_file in candidates:
        core = extract_core(out_file)

        unknown_count = sum(1 for bp in core["base_pairs"] if bp.get("kind") == "unknown")
        if unknown_count and not allow_unknown:
            errors.append(f"{out_file}: unknown base-pair lines ({unknown_count})")
            if not keep_going:
                break
            continue

        stats = core.get("stats", {})
        required_keys = ("total_pairs", "total_bases", "pair_type_counts")
        missing = [k for k in required_keys if k not in stats]
        if missing and not allow_missing_stats:
            errors.append(f"{out_file}: missing stats fields {missing}")
            if not keep_going:
                break
            continue

        rel_out = out_file.relative_to(root)
        rel_core = rel_out.with_suffix(".core.json")
        core_path = out_dir / rel_core

        if args.verbose:
            sys.stderr.write(f"{out_file} -> {core_path}\n")

        if not dry_run:
            core_path.parent.mkdir(parents=True, exist_ok=True)
            core_path.write_text(
                json.dumps(core, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n",
                encoding="utf-8",
            )

        manifest_entries.append(
            {
                "out": _rel(out_file),
                "core_json": _rel(core_path),
                "counts": {
                    "base_pairs": len(core["base_pairs"]),
                    "multiplets": len(core["multiplets"]),
                    "unknown_base_pairs": unknown_count,
                },
            }
        )

    manifest_entries.sort(key=lambda x: (x["out"], x["core_json"]))

    if not dry_run:
        manifest_path = out_dir / "manifest.json"
        manifest_path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "root": _rel(root),
                    "exclude_suffix": list(exclude_suffixes),
                    "entries": manifest_entries,
                },
                sort_keys=True,
                ensure_ascii=False,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    if errors:
        for line in errors[:50]:
            sys.stderr.write(line + "\n")
        if len(errors) > 50:
            sys.stderr.write(f"... and {len(errors) - 50} more\n")
        return 1

    sys.stderr.write(f"frozen={len(manifest_entries)} out_dir={out_dir}\n")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extract and compare RNAVIEW .out core sections "
            "(base pairs, multiplets, and summary statistics)."
        )
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_extract = sub.add_parser("extract", help="Print canonical JSON for core sections")
    p_extract.add_argument("out", help="Path to RNAVIEW .out file")
    p_extract.set_defaults(func=_cmd_extract)

    p_compare = sub.add_parser("compare", help="Compare two .out files semantically")
    p_compare.add_argument("left", help="Golden/expected .out")
    p_compare.add_argument("right", help="Candidate .out")
    p_compare.add_argument("--max-diffs", type=int, default=50, help="Max diff lines to print")
    p_compare.set_defaults(func=_cmd_compare)

    p_freeze = sub.add_parser(
        "freeze",
        help="Extract core JSON for a directory of .out files (used to build golden baselines)",
    )
    p_freeze.add_argument("root", help="Directory to scan (e.g. test)")
    p_freeze.add_argument(
        "--out-dir",
        default=None,
        help="Where to write *.core.json and manifest.json (default: <root>/golden_core)",
    )
    p_freeze.add_argument(
        "--exclude-suffix",
        action="append",
        default=["_sort.out", "_patt.out", ".out.out"],
        help="Skip .out files whose name ends with this suffix (repeatable)",
    )
    p_freeze.add_argument(
        "--allow-unknown",
        action="store_true",
        help="Allow non-matching base-pair lines (kind=unknown) in extracted core",
    )
    p_freeze.add_argument(
        "--allow-missing-stats",
        action="store_true",
        help="Allow missing total_pairs/total_bases/pair_type_counts in extracted core",
    )
    p_freeze.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue even if some files fail validation (still exits non-zero)",
    )
    p_freeze.add_argument("--dry-run", action="store_true", help="List files but do not write output")
    p_freeze.add_argument("--verbose", action="store_true", help="Print per-file mapping to stderr")
    p_freeze.set_defaults(func=_cmd_freeze)

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
