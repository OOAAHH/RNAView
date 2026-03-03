#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


def _json_dumps(obj: Any, *, indent: int | None = 2) -> str:
    if indent is None:
        return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n"
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, indent=indent) + "\n"


def _iter_rows(path: Path) -> Iterable[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield {k: (v or "") for k, v in row.items()}


def _is_true(v: str) -> bool:
    return v.strip() == "1"


def _as_int(v: str) -> int:
    try:
        return int(v.strip())
    except Exception:  # noqa: BLE001
        return 0


def _cmd_stats(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest).expanduser().resolve()
    if not manifest.exists():
        sys.stderr.write(f"missing manifest: {manifest}\n")
        return 2

    n_total = 0
    n_downloaded = 0
    counts = Counter()
    dna_chain_count = Counter()
    rna_chain_count = Counter()
    protein_chain_count = Counter()

    for row in _iter_rows(manifest):
        n_total += 1
        if row.get("status") == "downloaded":
            n_downloaded += 1

        for k in (
            "has_rna",
            "has_dna",
            "has_protein",
            "has_other_polymer",
            "is_pure_rna",
            "is_rna_complex",
            "is_prot_rna",
            "is_dna_rna",
            "is_prot_dna_rna",
        ):
            if _is_true(row.get(k, "")):
                counts[k] += 1

        dna_chain_count[_as_int(row.get("dna_chain_count", ""))] += 1
        rna_chain_count[_as_int(row.get("rna_chain_count", ""))] += 1
        protein_chain_count[_as_int(row.get("protein_chain_count", ""))] += 1

    out = {
        "schema_version": 1,
        "manifest": str(manifest),
        "rows": {"total": n_total, "downloaded": n_downloaded},
        "counts": dict(counts),
        "distributions": {
            "dna_chain_count": {str(k): v for k, v in sorted(dna_chain_count.items())},
            "rna_chain_count": {str(k): v for k, v in sorted(rna_chain_count.items())},
            "protein_chain_count": {str(k): v for k, v in sorted(protein_chain_count.items())},
        },
    }
    sys.stdout.write(_json_dumps(out, indent=2))
    return 0


@dataclass(frozen=True)
class SelectedRow:
    pdb_id: str
    mmcif_path: Path
    size_bytes: int | None


def _matches_kind(row: dict[str, str], kind: str) -> bool:
    if kind == "dna-rna":
        return _is_true(row.get("is_dna_rna", ""))
    if kind == "prot-dna-rna":
        return _is_true(row.get("is_prot_dna_rna", ""))
    if kind == "has-dna":
        return _is_true(row.get("has_dna", ""))
    raise ValueError(f"unknown kind: {kind}")


def _cmd_select(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest).expanduser().resolve()
    if not manifest.exists():
        sys.stderr.write(f"missing manifest: {manifest}\n")
        return 2

    kind = str(args.kind).strip().lower()
    want_status = str(args.status).strip()

    rows: list[SelectedRow] = []
    for row in _iter_rows(manifest):
        if want_status and row.get("status") != want_status:
            continue
        if not _matches_kind(row, kind):
            continue
        mmcif_path = (row.get("mmcif_path") or "").strip()
        if not mmcif_path:
            continue
        p = Path(mmcif_path)
        if args.existing_only and not p.exists():
            continue
        size = None
        if p.exists():
            try:
                size = p.stat().st_size
            except OSError:
                size = None
        rows.append(SelectedRow(pdb_id=(row.get("pdb_id") or "").strip(), mmcif_path=p, size_bytes=size))

    sort_key = str(args.sort).strip().lower()
    if sort_key == "none":
        pass
    elif sort_key == "size":
        rows.sort(key=lambda r: (r.size_bytes is None, r.size_bytes or 0, r.pdb_id, str(r.mmcif_path)))
    else:
        sys.stderr.write("--sort must be none|size\n")
        return 2

    if args.limit is not None:
        rows = rows[: int(args.limit)]

    fmt = str(args.format).strip().lower()
    out_path = Path(args.output).expanduser().resolve() if args.output else None
    out_f = out_path.open("w", encoding="utf-8") if out_path else sys.stdout
    try:
        if fmt == "paths":
            for r in rows:
                out_f.write(str(r.mmcif_path) + "\n")
            return 0
        if fmt == "tsv":
            out_f.write("pdb_id\tmmcif_path\tsize_bytes\n")
            for r in rows:
                out_f.write(f"{r.pdb_id}\t{r.mmcif_path}\t{'' if r.size_bytes is None else r.size_bytes}\n")
            return 0
        sys.stderr.write("--format must be paths|tsv\n")
        return 2
    finally:
        if out_path:
            out_f.close()


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Utilities for the mmCIF manifest TSV used by RNAView modernization")
    sub = p.add_subparsers(dest="cmd", required=True)

    ps = sub.add_parser("stats", help="Print aggregate stats for a manifest TSV (JSON)")
    ps.add_argument("--manifest", required=True, help="Path to pdb_rna_until_*.tsv (tab-separated)")
    ps.set_defaults(func=_cmd_stats)

    psel = sub.add_parser("select", help="Select rows from a manifest TSV and output a list of mmCIF paths")
    psel.add_argument("--manifest", required=True, help="Path to pdb_rna_until_*.tsv (tab-separated)")
    psel.add_argument(
        "--kind",
        required=True,
        choices=("dna-rna", "prot-dna-rna", "has-dna"),
        help="Which subset to select (based on boolean columns in the TSV)",
    )
    psel.add_argument(
        "--status",
        default="downloaded",
        help="Filter by status column (default: downloaded); pass empty string to disable",
    )
    psel.add_argument("--existing-only", action="store_true", help="Only keep rows whose mmcif_path exists on disk")
    psel.add_argument("--sort", default="size", choices=("none", "size"), help="Sort output rows (default: size)")
    psel.add_argument("--limit", type=int, default=None, help="Keep only the first N rows after sorting")
    psel.add_argument("--format", default="paths", choices=("paths", "tsv"), help="Output format (default: paths)")
    psel.add_argument("-o", "--output", default=None, help="Write output to a file (default: stdout)")
    psel.set_defaults(func=_cmd_select)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())

