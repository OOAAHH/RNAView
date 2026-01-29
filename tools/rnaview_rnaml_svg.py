#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path


RENDERER_VERSION = "rnaml-svg-v2"


def _hsb_to_rgb(h: float, s: float, v: float) -> tuple[float, float, float]:
    h = h % 1.0
    if s <= 0.0:
        return (v, v, v)
    i = int(h * 6.0)
    f = (h * 6.0) - i
    p = v * (1.0 - s)
    q = v * (1.0 - s * f)
    t = v * (1.0 - s * (1.0 - f))
    i = i % 6
    if i == 0:
        return (v, t, p)
    if i == 1:
        return (q, v, p)
    if i == 2:
        return (p, v, t)
    if i == 3:
        return (p, q, v)
    if i == 4:
        return (t, p, v)
    return (v, p, q)


def _rgb_hex(rgb: tuple[float, float, float]) -> str:
    r, g, b = rgb
    ri = max(0, min(255, int(round(r * 255.0))))
    gi = max(0, min(255, int(round(g * 255.0))))
    bi = max(0, min(255, int(round(b * 255.0))))
    return f"#{ri:02x}{gi:02x}{bi:02x}"


_BASE_COLORS = {
    # Copied from BASEPARS/ps_image.par (sethsbcolor).
    "A": _rgb_hex(_hsb_to_rgb(0.0000, 1.00, 1.00)),
    "C": _rgb_hex(_hsb_to_rgb(0.1667, 1.00, 1.00)),
    "G": _rgb_hex(_hsb_to_rgb(0.3333, 1.00, 1.00)),
    "U": _rgb_hex(_hsb_to_rgb(0.5000, 1.00, 1.00)),
    "T": _rgb_hex(_hsb_to_rgb(0.6667, 1.00, 1.00)),
    "I": _rgb_hex(_hsb_to_rgb(0.3333, 1.00, 0.57)),
    "X": _rgb_hex(_hsb_to_rgb(0.0000, 0.00, 0.00)),
}


def _fmt3(x: float) -> str:
    if abs(x) < 0.0005:
        x = 0.0
    return f"{x:.3f}"


@dataclass(frozen=True)
class BaseKey:
    mol: str
    pos: int


@dataclass(frozen=True)
class Base2D:
    key: BaseKey
    base: str
    x: float
    y: float


@dataclass(frozen=True)
class BasePair2D:
    a: BaseKey
    b: BaseKey
    edge_i: str
    edge_j: str
    orient: str


def _find_first(root: ET.Element, path: str) -> ET.Element | None:
    try:
        return root.find(path)
    except Exception:  # noqa: BLE001
        return None


def _findall(root: ET.Element, path: str) -> list[ET.Element]:
    try:
        return list(root.findall(path))
    except Exception:  # noqa: BLE001
        return []


def _mol_sort_key(mol_id: str) -> tuple[int, int | str]:
    try:
        return (0, int(mol_id))
    except ValueError:
        return (1, mol_id)


def _base_key_sort_key(key: BaseKey) -> tuple[int, int | str, int]:
    a, b = _mol_sort_key(key.mol)
    return (a, b, key.pos)


def _parse_base_key(base_id: ET.Element | None, *, default_mol: str | None) -> BaseKey | None:
    if base_id is None:
        return None
    mol = default_mol
    mol_node = _find_first(base_id, "./molecule-id")
    if mol_node is not None:
        ref = mol_node.get("ref")
        if ref:
            mol = ref
    if not mol:
        return None
    pos_text = _find_first(base_id, "./position")
    if pos_text is None:
        return None
    try:
        pos = int((pos_text.text or "").strip())
    except ValueError:
        return None
    return BaseKey(mol=mol, pos=pos)


def _parse_rnaml(xml_text: str) -> tuple[dict[str, list[Base2D]], list[BasePair2D]]:
    root = ET.fromstring(xml_text)

    # RNAML positions are 1..N per molecule (chain). Never merge across molecules.
    base_types: dict[BaseKey, str] = {}
    coords: dict[BaseKey, tuple[float, float]] = {}

    mol_elems = list(_findall(root, "./molecule"))
    mol_ids = [m.get("id") or "?" for m in mol_elems]

    for mol, mol_id in sorted(zip(mol_elems, mol_ids), key=lambda t: _mol_sort_key(t[1])):
        for base in _findall(mol, ".//structure//model//base"):
            pos_text = _find_first(base, "./position")
            typ_text = _find_first(base, "./base-type")
            if pos_text is None or typ_text is None:
                continue
            try:
                pos = int((pos_text.text or "").strip())
            except ValueError:
                continue
            bt = (typ_text.text or "").strip()
            if not bt:
                continue
            base_types[BaseKey(mol=mol_id, pos=pos)] = bt

        for node in _findall(mol, ".//secondary-structure-display//ss-base-coord"):
            pos_node = _find_first(node, "./base-id/position")
            coord_node = _find_first(node, "./coordinates")
            if pos_node is None or coord_node is None:
                continue
            try:
                pos = int((pos_node.text or "").strip())
            except ValueError:
                continue
            coord_text = (coord_node.text or "").strip()
            parts = coord_text.split()
            if len(parts) < 2:
                continue
            try:
                x = float(parts[0])
                y = float(parts[1])
            except ValueError:
                continue
            coords[BaseKey(mol=mol_id, pos=pos)] = (x, y)

    bases_by_mol: dict[str, list[Base2D]] = {}
    for key, (x, y) in sorted(coords.items(), key=lambda kv: _base_key_sort_key(kv[0])):
        bt = base_types.get(key, "?")
        bases_by_mol.setdefault(key.mol, []).append(Base2D(key=key, base=bt, x=x, y=y))

    bps: list[BasePair2D] = []

    # 1) Intra-molecule base pairs (typically lack <molecule-id> in <base-id>).
    for mol, mol_id in sorted(zip(mol_elems, mol_ids), key=lambda t: _mol_sort_key(t[1])):
        for bp in _findall(mol, ".//str-annotation//base-pair"):
            a_key = _parse_base_key(_find_first(bp, "./base-id-5p/base-id"), default_mol=mol_id)
            b_key = _parse_base_key(_find_first(bp, "./base-id-3p/base-id"), default_mol=mol_id)
            if a_key is None or b_key is None:
                continue
            edge_i = (bp.findtext("./edge-5p") or "").strip()
            edge_j = (bp.findtext("./edge-3p") or "").strip()
            orient = (bp.findtext("./bond-orientation") or "").strip()
            bps.append(BasePair2D(a=a_key, b=b_key, edge_i=edge_i, edge_j=edge_j, orient=orient))

    # 2) Interactions section (usually includes <molecule-id ref="...">).
    for bp in _findall(root, ".//interactions//str-annotation//base-pair"):
        a_key = _parse_base_key(_find_first(bp, "./base-id-5p/base-id"), default_mol=None)
        b_key = _parse_base_key(_find_first(bp, "./base-id-3p/base-id"), default_mol=None)
        if a_key is None or b_key is None:
            continue
        edge_i = (bp.findtext("./edge-5p") or "").strip()
        edge_j = (bp.findtext("./edge-3p") or "").strip()
        orient = (bp.findtext("./bond-orientation") or "").strip()
        bps.append(BasePair2D(a=a_key, b=b_key, edge_i=edge_i, edge_j=edge_j, orient=orient))

    # Deterministic ordering: endpoints first (molecule+position), then type fields.
    def bp_sort_key(r: BasePair2D) -> tuple[object, ...]:
        a = _base_key_sort_key(r.a)
        b = _base_key_sort_key(r.b)
        lo, hi = (a, b) if a <= b else (b, a)
        return (lo, hi, r.edge_i, r.edge_j, r.orient)

    bps.sort(key=bp_sort_key)
    return bases_by_mol, bps


def svg_from_rnaml(xml_text: str) -> str:
    bases_by_mol, bps = _parse_rnaml(xml_text)
    all_bases = [b for mol, bases in sorted(bases_by_mol.items(), key=lambda t: _mol_sort_key(t[0])) for b in bases]
    if not all_bases:
        raise ValueError("RNAML has no <secondary-structure-display>/<ss-base-coord> coordinates")

    by_key = {b.key: b for b in all_bases}
    xs = [b.x for b in all_bases]
    ys = [b.y for b in all_bases]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    margin = 10.0
    width = (max_x - min_x) + 2 * margin
    height = (max_y - min_y) + 2 * margin

    def map_xy(x: float, y: float) -> tuple[float, float]:
        # PS/RNAML uses y-up; SVG uses y-down. Flip to match legacy view.
        xm = (x - min_x) + margin
        ym = (max_y - y) + margin
        return xm, ym

    lines: list[str] = []
    lines.append('<?xml version="1.0" encoding="UTF-8"?>')
    lines.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 {_fmt3(width)} {_fmt3(height)}">'
    )
    lines.append("  <defs>")
    lines.append("    <style><![CDATA[")
    lines.append("      text.base { font-family: Times; font-weight: bold; font-size: 13px; text-anchor: middle; dominant-baseline: middle; }")
    lines.append("      line.backbone { stroke: #444444; stroke-width: 1; }")
    lines.append("      line.bp { stroke: #000000; stroke-width: 1; }")
    lines.append("      line.bp.tran { stroke-dasharray: 2 4; }")
    lines.append("      line.bp.unknown { stroke: #888888; stroke-dasharray: 1 3; }")
    lines.append("    ]]></style>")
    lines.append("  </defs>")

    # Backbone: connect consecutive bases if both have coordinates.
    lines.append('  <g id="backbone">')
    for mol_id, bases in sorted(bases_by_mol.items(), key=lambda t: _mol_sort_key(t[0])):
        lines.append(f'    <g id="backbone-mol-{mol_id}">')
        for a, b in zip(bases, bases[1:]):
            ax, ay = map_xy(a.x, a.y)
            bx, by = map_xy(b.x, b.y)
            lines.append(
                f'      <line class="backbone" x1="{_fmt3(ax)}" y1="{_fmt3(ay)}" x2="{_fmt3(bx)}" y2="{_fmt3(by)}" />'
            )
        lines.append("    </g>")
    lines.append("  </g>")

    # Base pairs.
    lines.append('  <g id="base-pairs">')
    for bp in bps:
        bi = by_key.get(bp.a)
        bj = by_key.get(bp.b)
        if bi is None or bj is None:
            continue
        x1, y1 = map_xy(bi.x, bi.y)
        x2, y2 = map_xy(bj.x, bj.y)

        cls = "bp"
        if bp.orient.lower().startswith("t"):
            cls = "bp tran"
        elif bp.orient.lower().startswith("c"):
            cls = "bp cis"
        else:
            cls = "bp unknown"

        title = f"{bp.a.mol}:{bp.a.pos}-{bp.b.mol}:{bp.b.pos} {bp.edge_i}{bp.edge_j}{bp.orient}".strip()
        lines.append(
            f'    <line class="{cls}" x1="{_fmt3(x1)}" y1="{_fmt3(y1)}" x2="{_fmt3(x2)}" y2="{_fmt3(y2)}"><title>{title}</title></line>'
        )
    lines.append("  </g>")

    # Bases (letters).
    lines.append('  <g id="bases">')
    for b in all_bases:
        x, y = map_xy(b.x, b.y)
        bt = b.base or "?"
        color = _BASE_COLORS.get(bt.upper(), _BASE_COLORS["X"])
        lines.append(
            f'    <text class="base" x="{_fmt3(x)}" y="{_fmt3(y)}" fill="{color}">{bt}</text>'
        )
    lines.append("  </g>")

    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Convert RNAVIEW RNAML(XML) to deterministic SVG")
    parser.add_argument("xml", help="Path to RNAML (.xml)")
    parser.add_argument("-o", "--output", required=True, help="Write SVG to this path")
    args = parser.parse_args(argv)

    xml_path = Path(args.xml)
    out_path = Path(args.output)
    text = xml_path.read_text(encoding="utf-8", errors="replace")
    svg = svg_from_rnaml(text)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(svg, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
