#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, replace
from pathlib import Path


RENDERER_VERSION = "ps-svg-v1"


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


def _fmt2(x: float) -> str:
    if abs(x) < 0.0005:
        x = 0.0
    return f"{x:.2f}"


def _fmt3(x: float) -> str:
    if abs(x) < 0.0005:
        x = 0.0
    return f"{x:.3f}"


@dataclass(frozen=True)
class _ColorMacro:
    color_hex: str
    does_fill: bool


@dataclass
class _PsToSvgState:
    llx: float
    lly: float
    urx: float
    ury: float

    tx: float = 0.0
    ty: float = 0.0

    color_hex: str = "#000000"
    line_width: float = 1.0
    dasharray: tuple[float, float] | None = None

    font_family: str = "Times-Bold"
    font_size: float = 13.0

    current_point: tuple[float, float] | None = None
    pending_moveto: tuple[float, float] | None = None

    path_d: str | None = None

    center_text: bool = False


def _strip_ps_comments(line: str) -> str:
    if "%" not in line:
        return line
    out: list[str] = []
    in_str = False
    esc = False
    for ch in line:
        if in_str:
            out.append(ch)
            if esc:
                esc = False
                continue
            if ch == "\\":
                esc = True
            elif ch == ")":
                in_str = False
            continue
        if ch == "%":
            break
        out.append(ch)
        if ch == "(":
            in_str = True
    return "".join(out)


def _tokenize_line(line: str) -> list[str]:
    s = _strip_ps_comments(line).strip()
    if not s:
        return []
    if s.startswith("%%"):
        return []
    tokens: list[str] = []
    i = 0
    n = len(s)
    while i < n:
        while i < n and s[i].isspace():
            i += 1
        if i >= n:
            break
        if s[i] == "(":
            i += 1
            buf: list[str] = []
            esc = False
            while i < n:
                ch = s[i]
                i += 1
                if esc:
                    buf.append(ch)
                    esc = False
                    continue
                if ch == "\\":
                    esc = True
                    continue
                if ch == ")":
                    break
                buf.append(ch)
            tokens.append(f"({''.join(buf)})")
            continue
        j = i
        while j < n and not s[j].isspace():
            j += 1
        tokens.append(s[i:j])
        i = j
    return tokens


_RE_BBOX = re.compile(r"^%%BoundingBox:\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*$")
_RE_SAT = re.compile(r"^/(MINOR_SAT|MAJOR_SAT)\s+([0-9]*\.?[0-9]+)\s+def\s*$")
_RE_WIDTH = re.compile(r"^/(W\d+|Dw)\s*{\s*([0-9]*\.?[0-9]+)\s+setlinewidth\s*}\s+bind\s+def\s*$")
_RE_COLOR_DEF = re.compile(r"^/([A-Za-z][A-Za-z0-9]*)\s*{\s*(.*?)\s*}\s+bind\s+def\s*$")


def _parse_ps_preamble(ps_text: str) -> tuple[tuple[float, float, float, float], dict[str, float], dict[str, float], dict[str, _ColorMacro]]:
    llx = lly = urx = ury = None  # type: ignore[assignment]
    vars_: dict[str, float] = {}
    widths: dict[str, float] = {}
    colors: dict[str, _ColorMacro] = {}

    for raw in ps_text.splitlines():
        line = raw.rstrip("\n")
        m = _RE_BBOX.match(line)
        if m:
            llx, lly, urx, ury = (float(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4)))
            continue

        m = _RE_SAT.match(line.strip())
        if m:
            vars_[m.group(1)] = float(m.group(2))
            continue

        m = _RE_WIDTH.match(line.strip())
        if m:
            widths[m.group(1)] = float(m.group(2))
            continue

        m = _RE_COLOR_DEF.match(line.strip())
        if not m:
            continue
        name = m.group(1)
        body = m.group(2).strip()
        body_tokens = body.split()
        if not body_tokens:
            continue
        does_fill = body_tokens[-1] == "fill"
        if does_fill:
            body_tokens = body_tokens[:-1]
        if len(body_tokens) >= 4 and body_tokens[-1] == "sethsbcolor":
            try:
                h = float(body_tokens[-4])
                s_token = body_tokens[-3]
                v_token = body_tokens[-2]
                s = float(vars_.get(s_token, s_token))
                v = float(vars_.get(v_token, v_token))
            except ValueError:
                continue
            colors[name] = _ColorMacro(color_hex=_rgb_hex(_hsb_to_rgb(h, s, v)), does_fill=does_fill)
            continue
        if len(body_tokens) >= 4 and body_tokens[-1] == "setrgbcolor":
            try:
                r = float(body_tokens[-4])
                g = float(body_tokens[-3])
                b = float(body_tokens[-2])
            except ValueError:
                continue
            colors[name] = _ColorMacro(color_hex=_rgb_hex((r, g, b)), does_fill=does_fill)
            continue
        if len(body_tokens) >= 2 and body_tokens[-1] == "setgray":
            try:
                gray = float(body_tokens[-2])
            except ValueError:
                continue
            colors[name] = _ColorMacro(color_hex=_rgb_hex((gray, gray, gray)), does_fill=does_fill)
            continue

    if llx is None or lly is None or urx is None or ury is None:
        raise ValueError("PS missing %%BoundingBox header")
    return (llx, lly, urx, ury), vars_, widths, colors


def _ps_to_svg_xy(state: _PsToSvgState, x: float, y: float) -> tuple[float, float]:
    xp = x + state.tx
    yp = y + state.ty
    # SVG uses y-down; map PS (y-up) into viewBox with origin at (llx, lly).
    return (xp - state.llx, state.ury - yp)


def _path_from_polygon(state: _PsToSvgState, points: list[tuple[float, float]]) -> str:
    if not points:
        return ""
    x0, y0 = _ps_to_svg_xy(state, points[0][0], points[0][1])
    parts = [f"M{_fmt2(x0)},{_fmt2(y0)}"]
    for x, y in points[1:]:
        xs, ys = _ps_to_svg_xy(state, x, y)
        parts.append(f"L{_fmt2(xs)},{_fmt2(ys)}")
    parts.append("Z")
    return " ".join(parts)


def _path_from_circle(state: _PsToSvgState, cx: float, cy: float, r: float) -> str:
    # Two arcs make a full circle.
    cxs, cys = _ps_to_svg_xy(state, cx, cy)
    rs = r
    x0 = cxs + rs
    y0 = cys
    x1 = cxs - rs
    y1 = cys
    return (
        f"M{_fmt2(x0)},{_fmt2(y0)} "
        f"A{_fmt2(rs)},{_fmt2(rs)} 0 1 0 {_fmt2(x1)},{_fmt2(y1)} "
        f"A{_fmt2(rs)},{_fmt2(rs)} 0 1 0 {_fmt2(x0)},{_fmt2(y0)} Z"
    )


def _svg_path_element(
    *,
    d: str,
    fill: str | None,
    stroke: str | None,
    stroke_width: float | None,
    dasharray: tuple[float, float] | None,
) -> str:
    attrs: list[str] = [f'd="{d}"']
    attrs.append(f'fill="{fill}"' if fill else 'fill="none"')
    if stroke:
        attrs.append(f'stroke="{stroke}"')
        if stroke_width is not None:
            attrs.append(f'stroke-width="{_fmt3(stroke_width)}"')
        attrs.append('stroke-linecap="round"')
        attrs.append('stroke-linejoin="round"')
        if dasharray:
            a, b = dasharray
            attrs.append(f'stroke-dasharray="{_fmt2(a)} {_fmt2(b)}"')
    else:
        attrs.append('stroke="none"')
    return f"<path {' '.join(attrs)} />"


def _svg_line_element(
    *,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    stroke: str,
    stroke_width: float,
    dasharray: tuple[float, float] | None,
) -> str:
    attrs: list[str] = [
        f'x1="{_fmt2(x1)}"',
        f'y1="{_fmt2(y1)}"',
        f'x2="{_fmt2(x2)}"',
        f'y2="{_fmt2(y2)}"',
        f'stroke="{stroke}"',
        f'stroke-width="{_fmt3(stroke_width)}"',
        'fill="none"',
        'stroke-linecap="round"',
        'stroke-linejoin="round"',
    ]
    if dasharray:
        a, b = dasharray
        attrs.append(f'stroke-dasharray="{_fmt2(a)} {_fmt2(b)}"')
    return f"<line {' '.join(attrs)} />"


def _svg_text_element(*, x: float, y: float, text: str, fill: str, font_family: str, font_size: float) -> str:
    safe = (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )
    return (
        f'<text x="{_fmt2(x)}" y="{_fmt2(y)}" fill="{fill}" '
        f'font-family="{font_family}" font-size="{_fmt2(font_size)}" '
        'text-anchor="middle" dominant-baseline="middle">'
        f"{safe}</text>"
    )


def svg_from_ps(ps_text: str) -> str:
    (llx, lly, urx, ury), _vars, widths, colors = _parse_ps_preamble(ps_text)
    width = urx - llx
    height = ury - lly

    state = _PsToSvgState(llx=llx, lly=lly, urx=urx, ury=ury)
    gsave_stack: list[_PsToSvgState] = []
    operand: list[object] = []
    out_elems: list[str] = []
    skipping_def = False

    def pop_num() -> float:
        if not operand:
            raise ValueError("PS stack underflow (number)")
        v = operand.pop()
        if isinstance(v, (int, float)):
            return float(v)
        raise ValueError(f"PS expected number, got: {v!r}")

    def pop_str() -> str:
        if not operand:
            raise ValueError("PS stack underflow (string)")
        v = operand.pop()
        if isinstance(v, str):
            return v
        raise ValueError(f"PS expected string, got: {v!r}")

    def do_newpath() -> None:
        state.path_d = None
        state.pending_moveto = None

    def do_fill() -> None:
        if not state.path_d:
            return
        out_elems.append(
            _svg_path_element(d=state.path_d, fill=state.color_hex, stroke=None, stroke_width=None, dasharray=None)
        )
        state.path_d = None
        state.pending_moveto = None

    def do_stroke() -> None:
        if not state.path_d:
            return
        out_elems.append(
            _svg_path_element(
                d=state.path_d,
                fill=None,
                stroke=state.color_hex,
                stroke_width=state.line_width,
                dasharray=state.dasharray,
            )
        )
        state.path_d = None
        state.pending_moveto = None

    for raw_line in ps_text.splitlines():
        line_stripped = raw_line.strip()
        if skipping_def:
            if "def" in line_stripped:
                skipping_def = False
            continue
        if line_stripped.startswith("/") and "{" in line_stripped and "findfont" not in line_stripped:
            if "def" not in line_stripped:
                skipping_def = True
            continue
        if line_stripped.startswith("/") and "def" in line_stripped and "findfont" not in line_stripped:
            continue

        tokens = _tokenize_line(raw_line)
        if not tokens:
            continue

        for tok in tokens:
            if tok.startswith("/") and len(tok) > 1:
                operand.append(tok[1:])
                continue

            if tok.startswith("(") and tok.endswith(")"):
                operand.append(tok[1:-1])
                continue

            try:
                operand.append(float(tok))
                continue
            except ValueError:
                pass

            if tok in widths:
                state.line_width = widths[tok]
                continue

            cm = colors.get(tok)
            if cm is not None:
                state.color_hex = cm.color_hex
                if cm.does_fill:
                    do_fill()
                continue

            if tok in ("NP", "newpath"):
                do_newpath()
                continue

            if tok == "gsave":
                gsave_stack.append(replace(state))
                continue

            if tok == "grestore":
                if not gsave_stack:
                    continue
                state = gsave_stack.pop()
                continue

            if tok == "translate":
                ty = pop_num()
                tx = pop_num()
                state.tx += tx
                state.ty += ty
                continue

            if tok == "setlinewidth":
                state.line_width = pop_num()
                continue

            if tok == "sethsbcolor":
                v = pop_num()
                s = pop_num()
                h = pop_num()
                state.color_hex = _rgb_hex(_hsb_to_rgb(h, s, v))
                continue

            if tok == "setrgbcolor":
                b = pop_num()
                g = pop_num()
                r = pop_num()
                state.color_hex = _rgb_hex((r, g, b))
                continue

            if tok == "setgray":
                g = pop_num()
                state.color_hex = _rgb_hex((g, g, g))
                continue

            if tok == "moveto":
                y = pop_num()
                x = pop_num()
                state.current_point = (x, y)
                state.pending_moveto = (x, y)
                continue

            if tok == "lineto":
                y = pop_num()
                x = pop_num()
                if state.pending_moveto is None:
                    state.pending_moveto = state.current_point
                if state.pending_moveto is not None:
                    x0, y0 = state.pending_moveto
                    x0s, y0s = _ps_to_svg_xy(state, x0, y0)
                    x1s, y1s = _ps_to_svg_xy(state, x, y)
                    if state.path_d is None:
                        state.path_d = f"M{_fmt2(x0s)},{_fmt2(y0s)} L{_fmt2(x1s)},{_fmt2(y1s)}"
                    else:
                        state.path_d = f"{state.path_d} L{_fmt2(x1s)},{_fmt2(y1s)}"
                    state.pending_moveto = (x0, y0)
                state.current_point = (x, y)
                continue

            if tok == "closepath":
                if state.path_d is not None:
                    state.path_d = f"{state.path_d} Z"
                continue

            if tok == "stroke":
                do_stroke()
                continue

            if tok == "fill":
                do_fill()
                continue

            if tok == "center":
                state.center_text = True
                continue

            if tok == "show":
                text = pop_str()
                if state.current_point is None:
                    continue
                x, y = state.current_point
                xs, ys = _ps_to_svg_xy(state, x, y)
                out_elems.append(
                    _svg_text_element(
                        x=xs, y=ys, text=text, fill=state.color_hex, font_family=state.font_family, font_size=state.font_size
                    )
                )
                state.center_text = False
                state.pending_moveto = None
                continue

            if tok == "findfont":
                # No-op placeholder (the literal font name is on the stack).
                continue

            if tok == "scalefont":
                size = pop_num()
                font = pop_str()
                operand.append((font, size))
                continue

            if tok == "setfont":
                if not operand:
                    continue
                v = operand.pop()
                if isinstance(v, tuple) and len(v) == 2:
                    font, size = v
                    if isinstance(font, str) and isinstance(size, (int, float)):
                        state.font_family = font
                        state.font_size = float(size)
                continue

            if tok == "LINE":
                y2 = pop_num()
                x2 = pop_num()
                y1 = pop_num()
                x1 = pop_num()
                x1s, y1s = _ps_to_svg_xy(state, x1, y1)
                x2s, y2s = _ps_to_svg_xy(state, x2, y2)
                out_elems.append(
                    _svg_line_element(
                        x1=x1s,
                        y1=y1s,
                        x2=x2s,
                        y2=y2s,
                        stroke=state.color_hex,
                        stroke_width=state.line_width,
                        dasharray=None,
                    )
                )
                continue

            if tok == "DASHLINE":
                y2 = pop_num()
                x2 = pop_num()
                y1 = pop_num()
                x1 = pop_num()
                x1s, y1s = _ps_to_svg_xy(state, x1, y1)
                x2s, y2s = _ps_to_svg_xy(state, x2, y2)
                out_elems.append(
                    _svg_line_element(
                        x1=x1s,
                        y1=y1s,
                        x2=x2s,
                        y2=y2s,
                        stroke=state.color_hex,
                        stroke_width=state.line_width,
                        dasharray=(2.0, 4.0),
                    )
                )
                continue

            if tok == "CIRCLE":
                r = pop_num()
                cy = pop_num()
                cx = pop_num()
                state.path_d = _path_from_circle(state, cx, cy, r)
                state.pending_moveto = None
                continue

            if tok == "TRIANGLE":
                y3 = pop_num()
                x3 = pop_num()
                y2 = pop_num()
                x2 = pop_num()
                y1 = pop_num()
                x1 = pop_num()
                state.path_d = _path_from_polygon(state, [(x1, y1), (x2, y2), (x3, y3)])
                state.pending_moveto = None
                continue

            if tok == "SQUARE":
                y4 = pop_num()
                x4 = pop_num()
                y3 = pop_num()
                x3 = pop_num()
                y2 = pop_num()
                x2 = pop_num()
                y1 = pop_num()
                x1 = pop_num()
                state.path_d = _path_from_polygon(state, [(x1, y1), (x2, y2), (x3, y3), (x4, y4)])
                state.pending_moveto = None
                continue

            if tok == "R6":
                pts: list[tuple[float, float]] = []
                for _ in range(6):
                    y = pop_num()
                    x = pop_num()
                    pts.append((x, y))
                pts.reverse()
                state.path_d = _path_from_polygon(state, pts)
                state.pending_moveto = None
                continue

            if tok == "R9":
                pts = []
                for _ in range(9):
                    y = pop_num()
                    x = pop_num()
                    pts.append((x, y))
                pts.reverse()
                state.path_d = _path_from_polygon(state, pts)
                state.pending_moveto = None
                continue

            if tok == "showpage":
                continue

            # Unknown token: ignore (legacy PS has lots of macro defs we don't need to evaluate).

    svg_lines: list[str] = []
    svg_lines.append('<?xml version="1.0" encoding="UTF-8"?>')
    svg_lines.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 {_fmt2(width)} {_fmt2(height)}">'
    )
    svg_lines.extend(f"  {e}" for e in out_elems)
    svg_lines.append("</svg>")
    return "\n".join(svg_lines) + "\n"


def _cmd_convert(args: argparse.Namespace) -> int:
    ps_path = Path(args.input).resolve()
    svg_path = Path(args.output).resolve()
    svg_text = svg_from_ps(ps_path.read_text(encoding="utf-8", errors="replace"))
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path.write_text(svg_text, encoding="utf-8")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Deterministic legacy RNAVIEW PostScript -> SVG converter")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("convert", help="Convert a legacy RNAVIEW .ps to .svg")
    c.add_argument("--input", required=True, help="Input .ps file")
    c.add_argument("--output", required=True, help="Output .svg file")
    c.set_defaults(func=_cmd_convert)

    ns = p.parse_args(argv)
    return int(ns.func(ns))


if __name__ == "__main__":
    raise SystemExit(main())
