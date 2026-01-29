#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import json
import math
import struct
from dataclasses import dataclass
from pathlib import Path


RENDERER_VERSION = "vrml-gltf-v1"


def _fmt_float(x: float, *, ndigits: int = 6) -> float:
    return float(round(x, ndigits))


def _axis_angle_to_quat(ax: float, ay: float, az: float, angle: float) -> tuple[float, float, float, float]:
    n = math.sqrt(ax * ax + ay * ay + az * az)
    if n <= 0.0:
        return (0.0, 0.0, 0.0, 1.0)
    ax /= n
    ay /= n
    az /= n
    s = math.sin(angle * 0.5)
    c = math.cos(angle * 0.5)
    return (_fmt_float(ax * s), _fmt_float(ay * s), _fmt_float(az * s), _fmt_float(c))


def _tokenize_vrml_body(text: str) -> list[str]:
    out: list[str] = []
    i = 0
    while i < len(text):
        while i < len(text) and text[i].isspace():
            i += 1
        if i >= len(text):
            break
        if text[i] == '"':
            i += 1
            start = i
            while i < len(text) and text[i] != '"':
                i += 1
            out.append(text[start:i])
            if i < len(text) and text[i] == '"':
                i += 1
            continue
        start = i
        while i < len(text) and not text[i].isspace():
            i += 1
        out.append(text[start:i])
    return out


def _to_float(tok: str) -> float:
    return float(tok.strip().rstrip(","))


def _to_int(tok: str) -> int:
    return int(tok.strip().rstrip(","))


@dataclass(frozen=True)
class VrmlInstance:
    kind: str
    fields: dict[str, object]


@dataclass(frozen=True)
class VrmlLineSet:
    points: list[tuple[float, float, float]]
    segments: list[list[int]]
    color: tuple[float, float, float]
    transparency: float


def _parse_instance_line(line: str) -> VrmlInstance | None:
    s = line.strip()
    for kind in ("Label", "Sphere", "Cyl", "Box", "Cone"):
        if not s.startswith(kind):
            continue
        # Most legacy instances are single-line: Kind { ... }.
        lb = s.find("{")
        rb = s.rfind("}")
        if lb < 0 or rb < 0 or rb <= lb:
            return None
        body = s[lb + 1 : rb].strip()
        tokens = _tokenize_vrml_body(body)
        fields: dict[str, object] = {}
        i = 0
        while i < len(tokens):
            key = tokens[i]
            i += 1
            if key == "p":
                fields["p"] = (_to_float(tokens[i]), _to_float(tokens[i + 1]), _to_float(tokens[i + 2]))
                i += 3
            elif key == "r":
                fields["r"] = (
                    _to_float(tokens[i]),
                    _to_float(tokens[i + 1]),
                    _to_float(tokens[i + 2]),
                    _to_float(tokens[i + 3]),
                )
                i += 4
            elif key == "s":
                fields["s"] = (_to_float(tokens[i]), _to_float(tokens[i + 1]), _to_float(tokens[i + 2]))
                i += 3
            elif key in ("dc", "ec", "sc"):
                fields[key] = (_to_float(tokens[i]), _to_float(tokens[i + 1]), _to_float(tokens[i + 2]))
                i += 3
            elif key in ("rad", "sz", "sh", "tr"):
                fields[key] = _to_float(tokens[i])
                i += 1
            elif key == "o":
                fields["o"] = (_to_float(tokens[i]), _to_float(tokens[i + 1]), _to_float(tokens[i + 2]))
                i += 3
            elif key == "c":
                fields["c"] = tokens[i]
                i += 1
            else:
                # Unknown key: skip one token as a best-effort fallback.
                if i < len(tokens):
                    i += 1
        return VrmlInstance(kind=kind, fields=fields)
    return None


def _parse_indexed_line_sets(lines: list[str], *, start_at: int) -> list[VrmlLineSet]:
    out: list[VrmlLineSet] = []

    i = start_at
    in_shape = False
    shape_depth = 0
    in_indexed = False
    indexed_depth = 0

    diffuse: tuple[float, float, float] | None = None
    transparency: float | None = None
    points: list[tuple[float, float, float]] = []
    coord_index: list[int] = []

    def flush() -> None:
        nonlocal diffuse, transparency, points, coord_index
        if not points or not coord_index:
            diffuse = None
            transparency = None
            points = []
            coord_index = []
            return
        segs: list[list[int]] = []
        cur: list[int] = []
        for idx in coord_index:
            if idx == -1:
                if len(cur) >= 2:
                    segs.append(cur)
                cur = []
                continue
            cur.append(idx)
        if len(cur) >= 2:
            segs.append(cur)

        out.append(
            VrmlLineSet(
                points=points,
                segments=segs,
                color=diffuse or (0.7, 0.7, 0.7),
                transparency=transparency if transparency is not None else 0.0,
            )
        )
        diffuse = None
        transparency = None
        points = []
        coord_index = []

    def add_points_from_text(text: str) -> None:
        nonlocal points
        toks = text.replace(",", " ").split()
        vals: list[float] = []
        for t in toks:
            try:
                vals.append(float(t))
            except ValueError:
                continue
        for j in range(0, len(vals) - 2, 3):
            points.append((_fmt_float(vals[j]), _fmt_float(vals[j + 1]), _fmt_float(vals[j + 2])))

    def add_indices_from_text(text: str) -> None:
        nonlocal coord_index
        toks = text.replace(",", " ").split()
        for t in toks:
            try:
                coord_index.append(int(t))
            except ValueError:
                continue

    while i < len(lines):
        line = lines[i]
        s = line.strip()

        if not in_shape and s.startswith("Shape {"):
            in_shape = True
            shape_depth = line.count("{") - line.count("}")
            diffuse = None
            transparency = None
            points = []
            coord_index = []
            in_indexed = False
            indexed_depth = 0
            i += 1
            continue

        if in_shape:
            shape_depth += line.count("{") - line.count("}")

            if s.startswith("diffuseColor "):
                parts = s.split()
                if len(parts) >= 4:
                    try:
                        diffuse = (_fmt_float(float(parts[1])), _fmt_float(float(parts[2])), _fmt_float(float(parts[3])))
                    except ValueError:
                        pass
            elif s.startswith("transparency "):
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        transparency = _fmt_float(float(parts[1]))
                    except ValueError:
                        pass

            if not in_indexed and "geometry IndexedLineSet" in s:
                in_indexed = True
                indexed_depth = line.count("{") - line.count("}")
                i += 1
                continue

            if in_indexed:
                indexed_depth += line.count("{") - line.count("}")

                if "point [" in s:
                    buf = s.split("point [", 1)[1]
                    while True:
                        if "]" in buf:
                            before = buf.split("]", 1)[0]
                            add_points_from_text(before)
                            break
                        add_points_from_text(buf)
                        i += 1
                        if i >= len(lines):
                            break
                        buf = lines[i].strip()

                if "coordIndex [" in s:
                    buf = s.split("coordIndex [", 1)[1]
                    while True:
                        if "]" in buf:
                            before = buf.split("]", 1)[0]
                            add_indices_from_text(before)
                            break
                        add_indices_from_text(buf)
                        i += 1
                        if i >= len(lines):
                            break
                        buf = lines[i].strip()

                if indexed_depth <= 0:
                    in_indexed = False

            if shape_depth <= 0:
                in_shape = False
                flush()

        i += 1

    return out


def _pad4(data: bytes) -> bytes:
    pad = (-len(data)) % 4
    if not pad:
        return data
    return data + (b"\x00" * pad)


def _pack_f32_vec3(values: list[tuple[float, float, float]]) -> bytes:
    out = bytearray()
    for x, y, z in values:
        out.extend(struct.pack("<fff", float(x), float(y), float(z)))
    return bytes(out)


def _pack_u16(values: list[int]) -> bytes:
    out = bytearray()
    for v in values:
        out.extend(struct.pack("<H", int(v)))
    return bytes(out)


def _min_max_vec3(values: list[tuple[float, float, float]]) -> tuple[list[float], list[float]]:
    xs = [v[0] for v in values]
    ys = [v[1] for v in values]
    zs = [v[2] for v in values]
    return [min(xs), min(ys), min(zs)], [max(xs), max(ys), max(zs)]


def _unit_box() -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float]], list[int]]:
    # A VRML Box defaults to size 2x2x2, centered at origin.
    # Use 24 vertices (4 per face) so normals are per-face.
    v: list[tuple[float, float, float]] = []
    n: list[tuple[float, float, float]] = []
    idx: list[int] = []

    faces = [
        ((0, 0, 1), [(-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]),
        ((0, 0, -1), [(1, -1, -1), (-1, -1, -1), (-1, 1, -1), (1, 1, -1)]),
        ((1, 0, 0), [(1, -1, 1), (1, -1, -1), (1, 1, -1), (1, 1, 1)]),
        ((-1, 0, 0), [(-1, -1, -1), (-1, -1, 1), (-1, 1, 1), (-1, 1, -1)]),
        ((0, 1, 0), [(-1, 1, 1), (1, 1, 1), (1, 1, -1), (-1, 1, -1)]),
        ((0, -1, 0), [(-1, -1, -1), (1, -1, -1), (1, -1, 1), (-1, -1, 1)]),
    ]

    for normal, quad in faces:
        base = len(v)
        for vx, vy, vz in quad:
            v.append((float(vx), float(vy), float(vz)))
            n.append((float(normal[0]), float(normal[1]), float(normal[2])))
        idx.extend([base + 0, base + 1, base + 2, base + 2, base + 3, base + 0])

    return v, n, idx


def _unit_cylinder(*, segments: int = 16) -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float]], list[int]]:
    # VRML Cylinder defaults: radius=1, height=2, centered at origin (y in [-1, 1]).
    v: list[tuple[float, float, float]] = []
    n: list[tuple[float, float, float]] = []
    idx: list[int] = []

    # Side vertices.
    for k in range(segments):
        a = (2.0 * math.pi * k) / segments
        x = math.cos(a)
        z = math.sin(a)
        v.append((x, -1.0, z))
        v.append((x, 1.0, z))
        n.append((x, 0.0, z))
        n.append((x, 0.0, z))

    for k in range(segments):
        i0 = 2 * k
        i1 = 2 * k + 1
        i2 = 2 * ((k + 1) % segments)
        i3 = 2 * ((k + 1) % segments) + 1
        idx.extend([i0, i2, i1, i1, i2, i3])

    # Caps.
    top_center = len(v)
    v.append((0.0, 1.0, 0.0))
    n.append((0.0, 1.0, 0.0))
    bot_center = len(v)
    v.append((0.0, -1.0, 0.0))
    n.append((0.0, -1.0, 0.0))

    for k in range(segments):
        k_next = (k + 1) % segments
        top_k = 2 * k + 1
        top_next = 2 * k_next + 1
        idx.extend([top_center, top_k, top_next])

        bot_k = 2 * k
        bot_next = 2 * k_next
        idx.extend([bot_center, bot_next, bot_k])

    return v, n, idx


def _unit_cone(*, segments: int = 16) -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float]], list[int]]:
    # VRML Cone defaults: bottomRadius=1, height=2, apex at y=+1, base at y=-1.
    v: list[tuple[float, float, float]] = []
    n: list[tuple[float, float, float]] = []
    idx: list[int] = []

    apex = (0.0, 1.0, 0.0)
    base_y = -1.0
    # Side faces.
    apex_index = 0
    v.append(apex)
    n.append((0.0, 0.0, 0.0))  # placeholder

    rim_start = len(v)
    for k in range(segments):
        a = (2.0 * math.pi * k) / segments
        x = math.cos(a)
        z = math.sin(a)
        v.append((x, base_y, z))
        # Approx side normal.
        # The cone's side normal has a y component; approximate for stable rendering.
        ny = 0.5
        nn = math.sqrt(x * x + ny * ny + z * z)
        n.append((x / nn, ny / nn, z / nn))

    # Fix apex normal as average (0,1,0) for stability.
    n[apex_index] = (0.0, 1.0, 0.0)

    for k in range(segments):
        a = rim_start + k
        b = rim_start + ((k + 1) % segments)
        idx.extend([apex_index, a, b])

    # Base cap.
    base_center = len(v)
    v.append((0.0, base_y, 0.0))
    n.append((0.0, -1.0, 0.0))
    for k in range(segments):
        a = rim_start + k
        b = rim_start + ((k + 1) % segments)
        idx.extend([base_center, b, a])

    return v, n, idx


def _unit_sphere(*, stacks: int = 8, slices: int = 16) -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float]], list[int]]:
    v: list[tuple[float, float, float]] = []
    n: list[tuple[float, float, float]] = []
    idx: list[int] = []

    for i in range(stacks + 1):
        phi = math.pi * i / stacks
        y = math.cos(phi)
        r = math.sin(phi)
        for j in range(slices + 1):
            theta = 2.0 * math.pi * j / slices
            x = r * math.cos(theta)
            z = r * math.sin(theta)
            v.append((_fmt_float(x), _fmt_float(y), _fmt_float(z)))
            n.append((_fmt_float(x), _fmt_float(y), _fmt_float(z)))

    def vid(i: int, j: int) -> int:
        return i * (slices + 1) + j

    for i in range(stacks):
        for j in range(slices):
            a = vid(i, j)
            b = vid(i + 1, j)
            c = vid(i + 1, j + 1)
            d = vid(i, j + 1)
            idx.extend([a, b, d, b, c, d])

    return v, n, idx


@dataclass
class _Geom:
    buffer: bytes
    accessors: dict[str, int]


class _GltfBuilder:
    def __init__(self) -> None:
        self.asset: dict[str, object] = {"version": "2.0", "generator": f"rnaview_vrml_gltf {RENDERER_VERSION}"}
        self.buffers: list[dict[str, object]] = []
        self.buffer_views: list[dict[str, object]] = []
        self.accessors: list[dict[str, object]] = []
        self.materials: list[dict[str, object]] = []
        self.meshes: list[dict[str, object]] = []
        self.nodes: list[dict[str, object]] = []

        self._material_map: dict[tuple[float, float, float, float], int] = {}
        self._geom_by_kind: dict[str, _Geom] = {}
        self._mesh_by_kind_material: dict[tuple[str, int], int] = {}

    def material(self, *, color: tuple[float, float, float], alpha: float) -> int:
        r, g, b = (_fmt_float(color[0], ndigits=3), _fmt_float(color[1], ndigits=3), _fmt_float(color[2], ndigits=3))
        a = _fmt_float(alpha, ndigits=3)
        key = (r, g, b, a)
        if key in self._material_map:
            return self._material_map[key]
        mat: dict[str, object] = {
            "pbrMetallicRoughness": {
                "baseColorFactor": [r, g, b, a],
                "metallicFactor": 0.0,
                "roughnessFactor": 1.0,
            }
        }
        if a < 1.0:
            mat["alphaMode"] = "BLEND"
        idx = len(self.materials)
        self.materials.append(mat)
        self._material_map[key] = idx
        return idx

    def _add_buffer(self, data: bytes) -> int:
        uri = "data:application/octet-stream;base64," + base64.b64encode(data).decode("ascii")
        idx = len(self.buffers)
        self.buffers.append({"byteLength": len(data), "uri": uri})
        return idx

    def _add_view(self, *, buffer: int, offset: int, length: int, target: int | None) -> int:
        view: dict[str, object] = {"buffer": buffer, "byteOffset": offset, "byteLength": length}
        if target is not None:
            view["target"] = target
        idx = len(self.buffer_views)
        self.buffer_views.append(view)
        return idx

    def _add_accessor(
        self,
        *,
        view: int,
        component_type: int,
        count: int,
        type_name: str,
        minv: list[float] | None = None,
        maxv: list[float] | None = None,
    ) -> int:
        acc: dict[str, object] = {
            "bufferView": view,
            "componentType": component_type,
            "count": count,
            "type": type_name,
        }
        if minv is not None:
            acc["min"] = minv
        if maxv is not None:
            acc["max"] = maxv
        idx = len(self.accessors)
        self.accessors.append(acc)
        return idx

    def _ensure_unit_geom(self, kind: str) -> _Geom:
        if kind in self._geom_by_kind:
            return self._geom_by_kind[kind]

        if kind == "Box":
            pos, nor, ind = _unit_box()
        elif kind == "Cyl":
            pos, nor, ind = _unit_cylinder()
        elif kind == "Cone":
            pos, nor, ind = _unit_cone()
        elif kind == "Sphere":
            pos, nor, ind = _unit_sphere()
        else:
            raise ValueError(f"unknown unit geom kind: {kind}")

        pos_bytes = _pack_f32_vec3(pos)
        nor_bytes = _pack_f32_vec3(nor)
        ind_bytes = _pack_u16(ind)

        blob = bytearray()
        off_pos = 0
        blob.extend(_pad4(pos_bytes))
        off_nor = len(blob)
        blob.extend(_pad4(nor_bytes))
        off_ind = len(blob)
        blob.extend(_pad4(ind_bytes))

        buf_idx = self._add_buffer(bytes(blob))
        pos_view = self._add_view(buffer=buf_idx, offset=off_pos, length=len(pos_bytes), target=34962)
        nor_view = self._add_view(buffer=buf_idx, offset=off_nor, length=len(nor_bytes), target=34962)
        ind_view = self._add_view(buffer=buf_idx, offset=off_ind, length=len(ind_bytes), target=34963)

        minv, maxv = _min_max_vec3(pos)
        pos_acc = self._add_accessor(view=pos_view, component_type=5126, count=len(pos), type_name="VEC3", minv=minv, maxv=maxv)
        nor_acc = self._add_accessor(view=nor_view, component_type=5126, count=len(nor), type_name="VEC3")
        ind_acc = self._add_accessor(view=ind_view, component_type=5123, count=len(ind), type_name="SCALAR")

        geom = _Geom(buffer=bytes(blob), accessors={"POSITION": pos_acc, "NORMAL": nor_acc, "INDICES": ind_acc})
        self._geom_by_kind[kind] = geom
        return geom

    def mesh_for_unit(self, *, kind: str, material: int) -> int:
        key = (kind, material)
        if key in self._mesh_by_kind_material:
            return self._mesh_by_kind_material[key]
        geom = self._ensure_unit_geom(kind)
        prim = {
            "attributes": {"POSITION": geom.accessors["POSITION"], "NORMAL": geom.accessors["NORMAL"]},
            "indices": geom.accessors["INDICES"],
            "material": material,
        }
        mesh_idx = len(self.meshes)
        self.meshes.append({"primitives": [prim]})
        self._mesh_by_kind_material[key] = mesh_idx
        return mesh_idx

    def mesh_for_lines(self, *, positions: list[tuple[float, float, float]], material: int) -> int:
        pos_bytes = _pack_f32_vec3(positions)
        blob = _pad4(pos_bytes)
        buf_idx = self._add_buffer(blob)
        pos_view = self._add_view(buffer=buf_idx, offset=0, length=len(pos_bytes), target=34962)
        minv, maxv = _min_max_vec3(positions)
        pos_acc = self._add_accessor(view=pos_view, component_type=5126, count=len(positions), type_name="VEC3", minv=minv, maxv=maxv)

        prim = {"attributes": {"POSITION": pos_acc}, "mode": 3, "material": material}
        mesh_idx = len(self.meshes)
        self.meshes.append({"primitives": [prim]})
        return mesh_idx

    def node(
        self,
        *,
        mesh: int | None = None,
        translation: tuple[float, float, float] | None = None,
        rotation: tuple[float, float, float, float] | None = None,
        scale: tuple[float, float, float] | None = None,
        extras: dict[str, object] | None = None,
    ) -> int:
        n: dict[str, object] = {}
        if mesh is not None:
            n["mesh"] = mesh
        if translation is not None:
            n["translation"] = [_fmt_float(translation[0]), _fmt_float(translation[1]), _fmt_float(translation[2])]
        if rotation is not None:
            n["rotation"] = [_fmt_float(rotation[0]), _fmt_float(rotation[1]), _fmt_float(rotation[2]), _fmt_float(rotation[3])]
        if scale is not None:
            n["scale"] = [_fmt_float(scale[0]), _fmt_float(scale[1]), _fmt_float(scale[2])]
        if extras:
            n["extras"] = extras
        idx = len(self.nodes)
        self.nodes.append(n)
        return idx


def gltf_from_vrml(vrml_text: str) -> str:
    lines = vrml_text.splitlines()
    start_at = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("Collision {"):
            start_at = i
            break

    line_sets = _parse_indexed_line_sets(lines, start_at=start_at)

    instances: list[VrmlInstance] = []
    for line in lines[start_at:]:
        inst = _parse_instance_line(line)
        if inst is not None:
            instances.append(inst)

    b = _GltfBuilder()

    # Root node (scene 0).
    root = b.node(extras={"source": "vrml", "renderer_version": RENDERER_VERSION})
    children: list[int] = []

    def add_child(n: int) -> None:
        children.append(n)

    # Lines first (stable and cheap).
    for ls in line_sets:
        alpha = 1.0 - float(ls.transparency)
        mat = b.material(color=ls.color, alpha=alpha)
        for seg in ls.segments:
            pts: list[tuple[float, float, float]] = []
            for idx in seg:
                if 0 <= idx < len(ls.points):
                    pts.append(ls.points[idx])
            if len(pts) < 2:
                continue
            mesh = b.mesh_for_lines(positions=pts, material=mat)
            n = b.node(mesh=mesh, extras={"kind": "IndexedLineSet"})
            add_child(n)

    # Instances.
    for inst in instances:
        f = inst.fields
        if inst.kind == "Label":
            p = f.get("p", (0.0, 0.0, 0.0))
            dc = f.get("dc", (0.0, 0.0, 0.0))
            tr = float(f.get("tr", 0.0))
            extras: dict[str, object] = {"kind": "Label"}
            if "c" in f:
                extras["text"] = str(f["c"])
            if "sz" in f:
                extras["size"] = float(f["sz"])
            if "o" in f:
                extras["offset"] = list(f["o"])  # type: ignore[arg-type]
            extras["color"] = [float(dc[0]), float(dc[1]), float(dc[2]), float(1.0 - tr)]  # type: ignore[index]
            n = b.node(translation=p, extras=extras)  # type: ignore[arg-type]
            add_child(n)
            continue

        dc = f.get("dc", (1.0, 0.0, 0.0))
        tr = float(f.get("tr", 0.0))
        mat = b.material(color=dc, alpha=1.0 - tr)  # type: ignore[arg-type]

        p = f.get("p", (0.0, 0.0, 0.0))
        if inst.kind == "Sphere":
            rad = float(f.get("rad", 1.0))
            mesh = b.mesh_for_unit(kind="Sphere", material=mat)
            n = b.node(mesh=mesh, translation=p, scale=(rad, rad, rad), extras={"kind": "Sphere"})  # type: ignore[arg-type]
            add_child(n)
            continue

        r = f.get("r", (0.0, 0.0, 1.0, 0.0))
        quat = _axis_angle_to_quat(float(r[0]), float(r[1]), float(r[2]), float(r[3]))  # type: ignore[index]
        s = f.get("s", (1.0, 1.0, 1.0))

        mesh = b.mesh_for_unit(kind=inst.kind, material=mat)
        n = b.node(mesh=mesh, translation=p, rotation=quat, scale=s, extras={"kind": inst.kind})  # type: ignore[arg-type]
        add_child(n)

    b.nodes[root]["children"] = children

    gltf: dict[str, object] = {
        "asset": b.asset,
        "scene": 0,
        "scenes": [{"nodes": [root]}],
        "nodes": b.nodes,
    }
    if b.meshes:
        gltf["meshes"] = b.meshes
    if b.materials:
        gltf["materials"] = b.materials
    if b.accessors:
        gltf["accessors"] = b.accessors
    if b.buffer_views:
        gltf["bufferViews"] = b.buffer_views
    if b.buffers:
        gltf["buffers"] = b.buffers

    return json.dumps(gltf, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Convert RNAVIEW legacy VRML (.wrl) to deterministic glTF (.gltf)")
    parser.add_argument("wrl", help="Path to legacy VRML (.wrl)")
    parser.add_argument("-o", "--output", required=True, help="Write glTF to this path")
    args = parser.parse_args(argv)

    wrl_path = Path(args.wrl)
    out_path = Path(args.output)
    text = wrl_path.read_text(encoding="utf-8", errors="replace")
    gltf = gltf_from_vrml(text)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(gltf, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

