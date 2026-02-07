from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from ._hotcore import run_hotcore


class UnsupportedInputError(ValueError):
    pass


def _infer_format(path: Path) -> str:
    ext = path.suffix.lower()
    if ext == ".cif":
        return "cif"
    if ext in (".pdb", ".ent"):
        return "pdb"
    raise UnsupportedInputError(f"unsupported input: {path}")


@dataclass(frozen=True)
class AnalyzeOutputs:
    pairs_json: Path
    out: Path
    xml: Path | None = None
    ps: Path | None = None
    wrl: Path | None = None


def from_structure(
    input_path: str | os.PathLike[str],
    *,
    out_pairs_json: str | os.PathLike[str],
    out_out: str | os.PathLike[str],
    semantics: str = "legacy-v1",
) -> None:
    input_path = Path(input_path).resolve()
    out_pairs_json = Path(out_pairs_json).resolve()
    out_out = Path(out_out).resolve()
    out_pairs_json.parent.mkdir(parents=True, exist_ok=True)
    out_out.parent.mkdir(parents=True, exist_ok=True)

    fmt = _infer_format(input_path)
    run_hotcore(
        [
            "from-structure",
            str(input_path),
            "--format",
            fmt,
            "--oracle",
            "compute",
            "--mmcif-parser",
            "legacy",
            "--semantics",
            semantics,
            "--hydrogen-policy",
            "legacy-mmcif-bug",
            "--missing-insertion-code",
            "legacy-question-mark",
            "--chain-id-policy",
            "legacy-1char",
            "-o",
            str(out_pairs_json),
            "--emit-out",
            str(out_out),
        ]
    )


def render_2d(
    pairs_json: str | os.PathLike[str],
    *,
    source_path: str | os.PathLike[str],
    out_xml: str | os.PathLike[str] | None = None,
    out_ps: str | os.PathLike[str] | None = None,
) -> None:
    pairs_json = Path(pairs_json).resolve()
    source_path = Path(source_path).resolve()
    out_xml_p = Path(out_xml).resolve() if out_xml is not None else None
    out_ps_p = Path(out_ps).resolve() if out_ps is not None else None
    if out_xml_p is None and out_ps_p is None:
        return
    if out_xml_p is not None:
        out_xml_p.parent.mkdir(parents=True, exist_ok=True)
    if out_ps_p is not None:
        out_ps_p.parent.mkdir(parents=True, exist_ok=True)

    args = ["render-2d", str(pairs_json), "--source", str(source_path)]
    if out_xml_p is not None:
        args += ["--out-xml", str(out_xml_p)]
    if out_ps_p is not None:
        args += ["--out-ps", str(out_ps_p)]
    run_hotcore(args)


def render_wrl(
    pairs_json: str | os.PathLike[str],
    *,
    source_path: str | os.PathLike[str],
    out_wrl: str | os.PathLike[str],
) -> None:
    pairs_json = Path(pairs_json).resolve()
    source_path = Path(source_path).resolve()
    out_wrl_p = Path(out_wrl).resolve()
    out_wrl_p.parent.mkdir(parents=True, exist_ok=True)
    run_hotcore(["render-wrl", str(pairs_json), "--source", str(source_path), "-o", str(out_wrl_p)])


def analyze(
    input_path: str | os.PathLike[str],
    *,
    out_dir: str | os.PathLike[str],
    formats: Iterable[str] = ("pairs", "out", "xml", "ps", "wrl"),
    semantics: str = "legacy-v1",
) -> AnalyzeOutputs:
    input_path = Path(input_path).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    formats_set = {str(f).strip().lower() for f in formats if str(f).strip()}

    pairs_json = out_dir / "pairs.json"
    out_out = out_dir / "engine.out"
    from_structure(input_path, out_pairs_json=pairs_json, out_out=out_out, semantics=semantics)

    out_xml = out_dir / f"{input_path.name}.xml" if "xml" in formats_set else None
    out_ps = out_dir / f"{input_path.name}.ps" if "ps" in formats_set else None
    if out_xml is not None or out_ps is not None:
        render_2d(pairs_json, source_path=input_path, out_xml=out_xml, out_ps=out_ps)

    out_wrl = out_dir / f"{input_path.name}.wrl" if "wrl" in formats_set else None
    if out_wrl is not None:
        render_wrl(pairs_json, source_path=input_path, out_wrl=out_wrl)

    return AnalyzeOutputs(
        pairs_json=pairs_json,
        out=out_out,
        xml=out_xml,
        ps=out_ps,
        wrl=out_wrl,
    )

