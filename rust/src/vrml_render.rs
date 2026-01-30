use crate::semantics::StructurePolicies;
use crate::structure::parse_structure_bases_with_atoms_with_policies;
use crate::{BasePair, PairsJson};
use std::fs;
use std::path::Path;

const XEPS: f64 = 1.0e-7;

#[derive(Debug, Clone)]
struct Residue3d {
    chain: char,
    resname_field: [char; 3], // legacy 3-char residue name field (right-aligned)
    c4: [f64; 3],
}

fn residue_name_field_legacy(resname: &str) -> [char; 3] {
    let mut s = resname.trim().to_ascii_uppercase();
    if s.len() == 2 && s.starts_with('D') {
        let mut chars = s.chars();
        let _ = chars.next();
        if let Some(second) = chars.next() {
            if matches!(second, 'A' | 'T' | 'G' | 'C') {
                s = second.to_string();
            }
        }
    }
    if s.len() > 3 {
        s = s.chars().take(3).collect();
    }
    let field = format!("{s:>3}");
    let mut it = field.chars();
    [
        it.next().unwrap_or(' '),
        it.next().unwrap_or(' '),
        it.next().unwrap_or(' '),
    ]
}

fn veclen(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn vec_norm(mut v: [f64; 3]) -> [f64; 3] {
    let l = veclen(v);
    if l > 0.0 {
        v[0] /= l;
        v[1] /= l;
        v[2] /= l;
    }
    v
}

fn identity_matrix() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

fn deg2rad(ang_deg: f64) -> f64 {
    ang_deg * std::f64::consts::PI / 180.0
}

fn rad2deg(ang_rad: f64) -> f64 {
    ang_rad * 180.0 / std::f64::consts::PI
}

fn dot2ang(dotval: f64) -> f64 {
    if dotval >= 1.0 {
        0.0
    } else if dotval <= -1.0 {
        180.0
    } else {
        rad2deg(dotval.acos())
    }
}

fn magang(va: [f64; 3], vb: [f64; 3]) -> f64 {
    if veclen(va) < XEPS || veclen(vb) < XEPS {
        return 0.0;
    }
    let va = vec_norm(va);
    let vb = vec_norm(vb);
    dot2ang(dot(va, vb))
}

fn arb_rotation(mut axis: [f64; 3], ang_deg: f64) -> [[f64; 3]; 3] {
    let vlen = veclen(axis);
    if vlen < XEPS {
        return identity_matrix();
    }
    axis[0] /= vlen;
    axis[1] /= vlen;
    axis[2] /= vlen;

    let ang = deg2rad(ang_deg);
    let c = ang.cos();
    let s = ang.sin();
    let dc = 1.0 - c;

    [
        [
            c + dc * axis[0] * axis[0],
            axis[0] * axis[1] * dc - axis[2] * s,
            axis[0] * axis[2] * dc + axis[1] * s,
        ],
        [
            axis[0] * axis[1] * dc + axis[2] * s,
            c + dc * axis[1] * axis[1],
            axis[1] * axis[2] * dc - axis[0] * s,
        ],
        [
            axis[0] * axis[2] * dc - axis[1] * s,
            axis[1] * axis[2] * dc + axis[0] * s,
            c + dc * axis[2] * axis[2],
        ],
    ]
}

fn rotate_point(rot: [[f64; 3]; 3], p: [f64; 3]) -> [f64; 3] {
    [dot(p, rot[0]), dot(p, rot[1]), dot(p, rot[2])]
}

fn covariance_matrix(points: &[[f64; 3]]) -> [[f64; 3]; 3] {
    let n = points.len();
    if n == 0 {
        return [[0.0; 3]; 3];
    }
    let mut mean = [0.0; 3];
    for p in points {
        mean[0] += p[0];
        mean[1] += p[1];
        mean[2] += p[2];
    }
    let nf = n as f64;
    mean[0] /= nf;
    mean[1] /= nf;
    mean[2] /= nf;

    let mut cov = [[0.0; 3]; 3];
    for p in points {
        let dx = [p[0] - mean[0], p[1] - mean[1], p[2] - mean[2]];
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += dx[i] * dx[j];
            }
        }
    }

    let denom = (n as f64) - 1.0;
    if denom.abs() < XEPS {
        return cov;
    }
    for i in 0..3 {
        for j in 0..3 {
            cov[i][j] /= denom;
        }
    }
    cov
}

fn rotate_mtx_entry(mtx: &mut [[f64; 3]; 3], i: usize, j: usize, k: usize, l: usize, s: f64, tau: f64) {
    let g = mtx[i][j];
    let h = mtx[k][l];
    mtx[i][j] = g - s * (h + g * tau);
    mtx[k][l] = h + s * (g - h * tau);
}

fn eigsrt(d: &mut [f64; 3], v: &mut [[f64; 3]; 3]) {
    for i in 0..2 {
        let mut k = i;
        let mut p = d[i];
        for j in (i + 1)..3 {
            if d[j] < p {
                p = d[j];
                k = j;
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for row in 0..3 {
                v[row].swap(i, k);
            }
        }
    }
}

fn jacobi(mut a: [[f64; 3]; 3]) -> Result<([f64; 3], [[f64; 3]; 3]), String> {
    let mut v = identity_matrix();
    let mut d = [a[0][0], a[1][1], a[2][2]];
    let mut b = d;
    let mut z = [0.0f64; 3];

    for iter in 0..100 {
        let mut sm = 0.0;
        for ip in 0..2 {
            for iq in (ip + 1)..3 {
                sm += a[ip][iq].abs();
            }
        }
        if sm < XEPS {
            eigsrt(&mut d, &mut v);
            return Ok((d, v));
        }

        let tresh = if iter < 3 { 0.2 * sm / 9.0 } else { 0.0 };
        for ip in 0..2 {
            for iq in (ip + 1)..3 {
                let g = 100.0 * a[ip][iq].abs();
                if iter > 3 && (d[ip].abs() + g) == d[ip].abs() && (d[iq].abs() + g) == d[iq].abs() {
                    a[ip][iq] = 0.0;
                    continue;
                }
                if a[ip][iq].abs() <= tresh {
                    continue;
                }

                let h = d[iq] - d[ip];
                let mut t = if (h.abs() + g) == h.abs() {
                    a[ip][iq] / h
                } else {
                    let theta = 0.5 * h / a[ip][iq];
                    let mut t = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
                    if theta < 0.0 {
                        t = -t;
                    }
                    t
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                let tau = s / (1.0 + c);
                let h = t * a[ip][iq];
                z[ip] -= h;
                z[iq] += h;
                d[ip] -= h;
                d[iq] += h;
                a[ip][iq] = 0.0;

                for j in 0..ip {
                    rotate_mtx_entry(&mut a, j, ip, j, iq, s, tau);
                }
                for j in (ip + 1)..iq {
                    rotate_mtx_entry(&mut a, ip, j, j, iq, s, tau);
                }
                for j in (iq + 1)..3 {
                    rotate_mtx_entry(&mut a, ip, j, iq, j, s, tau);
                }
                for j in 0..3 {
                    rotate_mtx_entry(&mut v, j, ip, j, iq, s, tau);
                }
            }
        }
        for ip in 0..3 {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }

    Err("too many iterations in jacobi".to_string())
}

fn ls_plane_normal(points: &[[f64; 3]]) -> [f64; 3] {
    if points.len() < 3 {
        return [0.0, 0.0, 1.0];
    }

    let cov = covariance_matrix(points);
    let Ok((_d, v)) = jacobi(cov) else {
        return [0.0, 0.0, 1.0];
    };

    let mut nml = false;
    for i in 0..3 {
        for j in 0..3 {
            let ident = if i == j { 1.0 } else { 0.0 };
            if (v[i][j] - ident).abs() > XEPS {
                nml = true;
                break;
            }
        }
        if nml {
            break;
        }
    }

    let mut normal = if nml {
        [v[0][0], v[1][0], v[2][0]]
    } else {
        [0.0, 0.0, 1.0]
    };

    if normal[2] < 0.0 {
        normal[0] = -normal[0];
        normal[1] = -normal[1];
        normal[2] = -normal[2];
    }
    normal
}

fn rot_2_lsplane_matrix(points: &[[f64; 3]]) -> [[f64; 3]; 3] {
    let zphy = [0.0, 0.0, 1.0];
    let normal = ls_plane_normal(points);
    let hinge = cross(normal, zphy);
    let ang = magang(zphy, normal);
    arb_rotation(hinge, ang)
}

fn pair_type3_from_pairs(bp: &BasePair) -> Option<[char; 3]> {
    if bp.kind == "stacked" {
        return None;
    }
    let lw = bp.lw.as_deref()?.trim();
    let mut lw_chars = lw.chars();
    let a = lw_chars.next()?;
    let _slash = lw_chars.next()?;
    let b = lw_chars.next()?;
    let orient = bp.orientation.as_deref().unwrap_or("cis").trim().to_ascii_lowercase();
    let c = if orient.starts_with("tran") { 't' } else { 'c' };
    Some([a, b, c])
}

#[derive(Debug, Clone)]
struct VrmlWriter {
    out: String,
    indent_flag: bool,
    indent: usize,
    max_buflen: usize,
    buflen: usize,
    needs_delimiter: bool,
}

impl VrmlWriter {
    fn new() -> Self {
        Self {
            out: String::new(),
            indent_flag: true,
            indent: 0,
            max_buflen: 70,
            buflen: 0,
            needs_delimiter: false,
        }
    }

    fn finish(self) -> String {
        self.out
    }

    fn check(&mut self) {
        if self.indent_flag && self.buflen >= self.max_buflen {
            self.newline();
        }
    }

    fn put_delimiter(&mut self) {
        if self.needs_delimiter {
            self.out.push(' ');
            self.buflen += 1;
            self.needs_delimiter = false;
        }
    }

    fn pad_blanks(&mut self, n: usize) {
        if self.indent_flag && n > 0 {
            for _ in 0..n {
                self.out.push(' ');
            }
            self.buflen += n;
        }
    }

    fn newline(&mut self) {
        if self.buflen != 0 {
            self.out.push('\n');
        }
        self.buflen = 0;
        self.needs_delimiter = false;
        self.pad_blanks(self.indent);
    }

    fn add_indent(&mut self, incr: i32) {
        if incr == 0 {
            self.indent = 0;
            return;
        }
        self.indent = (self.indent as i32 + incr).max(0) as usize;
    }

    fn s_newline(&mut self, s: &str) {
        self.check();
        self.put_delimiter();
        self.out.push_str(s);
        self.out.push('\n');
        self.buflen = 0;
        self.needs_delimiter = false;
        self.pad_blanks(self.indent);
    }

    fn s(&mut self, s: &str) {
        self.check();
        self.put_delimiter();
        self.out.push_str(s);
        self.buflen += s.len();
        let last = s.as_bytes().last().copied().unwrap_or(b' ');
        self.needs_delimiter = (last as char).is_ascii_alphanumeric();
    }

    fn node(&mut self, node: &str) {
        self.check();
        self.put_delimiter();
        self.out.push_str(node);
        self.out.push_str(" {");
        self.buflen += node.len() + 2;
        self.add_indent(2);
        self.newline();
    }

    fn finish_node(&mut self) {
        self.add_indent(-2);
        self.newline();
        self.out.push('}');
        self.buflen += 1;
        self.needs_delimiter = false;
    }

    fn list(&mut self, list: &str) {
        self.check();
        self.put_delimiter();
        self.out.push_str(list);
        self.out.push_str(" [");
        self.buflen += list.len() + 2;
        self.add_indent(2);
        self.newline();
    }

    fn finish_list(&mut self) {
        self.add_indent(-2);
        self.newline();
        self.out.push(']');
        self.buflen += 1;
        self.needs_delimiter = false;
    }

    fn i(&mut self, i: i32) {
        self.check();
        self.put_delimiter();
        let s = i.to_string();
        self.out.push_str(&s);
        self.buflen += s.len();
        self.needs_delimiter = true;
    }

    fn qs(&mut self, s: &str) {
        self.check();
        self.put_delimiter();
        self.out.push('"');
        self.out.push_str(s);
        self.out.push('"');
        self.buflen += s.len() + 2;
        self.needs_delimiter = true;
    }

    fn v3(&mut self, p: [f64; 3]) {
        self.check();
        if self.indent_flag {
            if self.needs_delimiter {
                self.out.push_str(", ");
                self.buflen += 2;
            }
        } else if self.needs_delimiter {
            self.out.push(' ');
            self.buflen += 1;
        }
        let s = format!("{:.2} {:.2} {:.2}", p[0], p[1], p[2]);
        self.out.push_str(&s);
        self.buflen += s.len();
        self.needs_delimiter = true;
    }

    fn rgb_color(&mut self, rgb: [f64; 3]) {
        self.check();
        if self.indent_flag {
            if self.needs_delimiter {
                self.out.push_str(", ");
                self.buflen += 2;
            }
        } else if self.needs_delimiter {
            self.out.push(' ');
            self.buflen += 1;
        }
        let s = format!("{} {} {}", fmt_g2(rgb[0]), fmt_g2(rgb[1]), fmt_g2(rgb[2]));
        self.out.push_str(&s);
        self.buflen += s.len();
        self.needs_delimiter = true;
    }
}

fn fmt_g2(x: f64) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    if x == 1.0 {
        return "1".to_string();
    }
    let abs = x.abs();
    if abs == 0.0 {
        return "0".to_string();
    }
    let exp = abs.log10().floor() as i32;
    let sig_digits: i32 = 2;
    let decimals = (sig_digits - 1 - exp).max(0) as usize;
    let scale = 10_f64.powi(decimals as i32);
    let rounded = (x * scale).round() / scale;
    let mut s = format!("{:.*}", decimals, rounded);
    if let Some(dot) = s.find('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.len() == dot + 1 {
            s.pop();
        }
    }
    if s.is_empty() {
        "0".to_string()
    } else {
        s
    }
}

fn vrml_material_color(chain_idx: usize) -> [f64; 3] {
    const COLORS: [[f64; 3]; 20] = [
        [0.263, 0.882, 0.341],
        [0.713, 0.337, 0.875],
        [1.000, 1.000, 0.500],
        [0.537, 0.000, 0.537],
        [0.500, 0.000, 0.000],
        [0.800, 0.500, 1.000],
        [0.086, 0.337, 0.282],
        [1.000, 0.500, 0.500],
        [0.000, 1.000, 1.000],
        [0.208, 0.486, 0.784],
        [0.000, 0.800, 0.000],
        [0.800, 0.800, 0.000],
        [1.000, 1.000, 1.000],
        [0.500, 1.000, 1.000],
        [0.767, 0.000, 0.767],
        [0.400, 0.400, 0.400],
        [1.000, 1.000, 0.700],
        [0.970, 0.970, 0.970],
        [0.000, 0.000, 1.000],
        [0.000, 0.000, 0.000],
    ];
    if chain_idx < 20 {
        COLORS[chain_idx]
    } else {
        [0.0, 0.0, 0.0]
    }
}

fn output_appearance(w: &mut VrmlWriter, chain_idx: usize) {
    // vrml_cond_newline
    if w.indent_flag && w.buflen != w.indent {
        w.out.push('\n');
        w.buflen = 0;
        w.needs_delimiter = false;
        w.pad_blanks(w.indent);
    }

    w.node("appearance Appearance");
    w.s("material");
    if chain_idx != 0 {
        vrml_material(w, chain_idx, 1.0);
    } else {
        vrml_material(w, chain_idx, 0.5);
    }
    w.finish_node();
}

fn vrml_material(w: &mut VrmlWriter, chain_idx: usize, tr: f64) {
    w.node("Material");
    w.check();
    if chain_idx < 20 {
        let rgb = vrml_material_color(chain_idx);
        let s = format!(
            "diffuseColor {} {} {}",
            fmt_g2(rgb[0]),
            fmt_g2(rgb[1]),
            fmt_g2(rgb[2])
        );
        w.out.push_str(&s);
        w.buflen += s.len();
    } else {
        w.out.push_str("diffuseColor 0 0 0 ");
        w.buflen += "diffuseColor 0 0 0 ".len();
    }

    w.check();
    w.newline();
    let s = format!("transparency {}", fmt_g2(tr));
    w.out.push_str(&s);
    w.buflen += s.len();
    w.finish_node();
}

fn vrml_draw_chain(w: &mut VrmlWriter, chain_idx: usize, chain_ranges: &[(usize, usize)], c4xyz: &[[f64; 3]]) {
    w.newline();
    w.node("Shape");
    output_appearance(w, chain_idx);
    w.newline();
    w.node("geometry IndexedLineSet");
    w.s("coord");
    w.node("Coordinate");
    w.list("point");

    let (start, end) = chain_ranges[chain_idx - 1];
    for i in start..=end {
        w.v3(c4xyz[i]);
    }

    w.finish_list();
    w.finish_node();
    w.newline();
    w.list("coordIndex");
    let resnum = (end as i32) - (start as i32);
    for i in 0..=resnum {
        w.i(i);
    }
    w.i(-1);
    w.finish_list();
    w.newline();
    w.s_newline("colorPerVertex FALSE");
    w.node("color Color");
    w.list("color");
    let rgb = vrml_material_color(chain_idx);
    w.rgb_color(rgb);
    w.finish_list();
    w.finish_node();
    w.newline();
    w.list("colorIndex");
    w.i(0);
    w.finish_list();
    w.finish_node();
    w.finish_node();
}

fn label_residue(out: &mut String, residues: &[Residue3d], start: usize, end: usize) {
    out.push('\n');
    for idx in start..=end {
        let r = &residues[idx];
        let c = r.resname_field[2];
        let res_color = match c {
            'A' => "dc  1.000, 0.000, 0.000",
            'U' => "dc  0.000, 1.000, 1.000",
            'G' => "dc  0.000, 1.000, 0.000",
            'C' => "dc  0.000, 0.000, 1.000",
            'I' => "dc  0.086, 0.337, 0.282",
            'T' => "dc  0.537, 0.000, 0.537",
            _ => "dc  0.970, 0.970, 0.970",
        };
        let label: String = r.resname_field.iter().collect();
        out.push_str(&format!(
            "Label {{ p  {:.2} {:.2} {:.2} c \"{}\" {}}}\n",
            r.c4[0], r.c4[1], r.c4[2], label, res_color
        ));
    }
}

fn draw_interact(out: &mut String, k1: usize, k2: usize, pair_type: [char; 3], c4xyz: &[[f64; 3]]) {
    let mut xyz1 = c4xyz[k1];
    let mut xyz2 = c4xyz[k2];
    let xyz0 = [(xyz1[0] + xyz2[0]) * 0.5, (xyz1[1] + xyz2[1]) * 0.5, (xyz1[2] + xyz2[2]) * 0.5];
    let xyz01 = [(xyz1[0] + xyz0[0]) * 0.5, (xyz1[1] + xyz0[1]) * 0.5, (xyz1[2] + xyz0[2]) * 0.5];
    let xyz02 = [(xyz2[0] + xyz0[0]) * 0.5, (xyz2[1] + xyz0[1]) * 0.5, (xyz2[2] + xyz0[2]) * 0.5];
    let xyz = [
        [(xyz1[0] + xyz01[0]) * 0.5, (xyz1[1] + xyz01[1]) * 0.5, (xyz1[2] + xyz01[2]) * 0.5],
        [(xyz01[0] + xyz0[0]) * 0.5, (xyz01[1] + xyz0[1]) * 0.5, (xyz01[2] + xyz0[2]) * 0.5],
        [(xyz0[0] + xyz02[0]) * 0.5, (xyz0[1] + xyz02[1]) * 0.5, (xyz0[2] + xyz02[2]) * 0.5],
        [(xyz02[0] + xyz2[0]) * 0.5, (xyz02[1] + xyz2[1]) * 0.5, (xyz02[2] + xyz2[2]) * 0.5],
    ];

    let (dc, tr) = if pair_type[2] == 't' {
        (" dc 1 1 0", 0.6)
    } else {
        (" dc 1 0 0", 0.0)
    };

    let vec = [xyz2[0] - xyz1[0], xyz2[1] - xyz1[1], xyz2[2] - xyz1[2]];
    let d = veclen(vec);
    let vector = if d > 0.0 {
        [vec[0] / d, vec[1] / d, vec[2] / d]
    } else {
        [0.0, 1.0, 0.0]
    };

    let yphy = [0.0, 1.0, 0.0];
    let angle = dot(vector, yphy).clamp(-1.0, 1.0).acos();
    let rot = if angle > 0.0 { vec_norm(cross(yphy, vector)) } else { yphy };

    let box_scale = [d / 18.0, d / 18.0, d / 18.0];
    let cone_scale = [d / 15.0, d / 15.0, d / 15.0];
    let mut cyl_scale = [d / 70.0, d / 2.2, d / 70.0];
    let radius = d / 15.0;

    fn output_shape(out: &mut String, kind: &str, tran: [f64; 3]) {
        out.push_str(&format!("{} {{p {:.2} {:.2} {:.2} ", kind, tran[0], tran[1], tran[2]));
    }
    fn output_rot(out: &mut String, rot: [f64; 3], angle: f64) {
        out.push_str(&format!(" r  {:.2} {:.2} {:.2} {:.2} ", rot[0], rot[1], rot[2], angle));
    }
    fn output_scale(out: &mut String, scale: [f64; 3]) {
        out.push_str(&format!(" s  {:.2} {:.2} {:.2} ", scale[0], scale[1], scale[2]));
    }
    fn output_radius(out: &mut String, rad: f64) {
        out.push_str(&format!(" rad  {:.2} ", rad));
    }
    fn output_color(out: &mut String, dc: &str) {
        out.push_str(&format!(" {} ", dc));
    }
    fn output_tr(out: &mut String, tr: f64) {
        out.push_str(&format!(" tr  {:.2}}}\n ", tr));
    }

    if pair_type[0] == '.' && pair_type[1] == '.' {
        cyl_scale[1] *= 0.2;
        for i in 0..4 {
            output_shape(out, "Cyl", xyz[i]);
            output_rot(out, rot, angle);
            output_scale(out, cyl_scale);
            output_tr(out, 0.0);
        }
        return;
    }

    if (pair_type[0] == '-' || pair_type[0] == '+') && (pair_type[1] == '-' || pair_type[1] == '+') {
        output_shape(out, "Cyl", xyz0);
        output_rot(out, rot, angle);
        output_scale(out, cyl_scale);
        output_tr(out, 0.0);
    } else if pair_type[0] != '.' && pair_type[1] != '.' {
        cyl_scale[0] = d / 80.0;
        cyl_scale[1] = d / 2.2;
        cyl_scale[2] = d / 80.0;
        output_shape(out, "Cyl", xyz0);
        output_rot(out, rot, angle);
        output_scale(out, cyl_scale);
        output_tr(out, 0.6);
    } else {
        cyl_scale[1] *= 0.2;
        for i in 0..4 {
            output_shape(out, "Cyl", xyz[i]);
            output_rot(out, rot, angle);
            output_scale(out, cyl_scale);
            output_tr(out, 0.0);
        }
    }

    // Edge decorations (match legacy draw_interact switch).
    match (pair_type[0], pair_type[1]) {
        ('W', 'W') => {
            output_shape(out, "Sphere", xyz0);
            output_radius(out, radius);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('H', 'H') => {
            output_shape(out, "Box", xyz0);
            output_rot(out, rot, angle);
            output_scale(out, box_scale);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('S', 'S') => {
            output_shape(out, "Cone", xyz0);
            output_rot(out, rot, angle);
            output_scale(out, cone_scale);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('W', 'H') => {
            output_shape(out, "Sphere", xyz01);
            output_radius(out, radius);
            output_color(out, dc);
            output_tr(out, tr);
            output_shape(out, "Box", xyz02);
            output_rot(out, rot, angle);
            output_scale(out, box_scale);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('H', 'W') => {
            output_shape(out, "Box", xyz01);
            output_rot(out, rot, angle);
            output_scale(out, box_scale);
            output_color(out, dc);
            output_tr(out, tr);
            output_shape(out, "Sphere", xyz02);
            output_radius(out, radius);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('W', 'S') => {
            output_shape(out, "Sphere", xyz01);
            output_radius(out, radius);
            output_color(out, dc);
            output_tr(out, tr);
            output_shape(out, "Cone", xyz02);
            output_rot(out, rot, angle);
            output_scale(out, cone_scale);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('S', 'W') => {
            output_shape(out, "Cone", xyz01);
            output_rot(out, rot, angle);
            output_scale(out, cone_scale);
            output_color(out, dc);
            output_tr(out, tr);
            output_shape(out, "Sphere", xyz02);
            output_radius(out, radius);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('H', 'S') => {
            output_shape(out, "Box", xyz01);
            output_rot(out, rot, angle);
            output_scale(out, box_scale);
            output_color(out, dc);
            output_tr(out, tr);
            output_shape(out, "Cone", xyz02);
            output_rot(out, rot, angle);
            output_scale(out, cone_scale);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('S', 'H') => {
            output_shape(out, "Cone", xyz01);
            output_rot(out, rot, angle);
            output_scale(out, cone_scale);
            output_color(out, dc);
            output_tr(out, tr);
            output_shape(out, "Box", xyz02);
            output_rot(out, rot, angle);
            output_scale(out, box_scale);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('W', '.') => {
            output_shape(out, "Sphere", xyz01);
            output_radius(out, radius);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('.', 'W') => {
            output_shape(out, "Sphere", xyz02);
            output_radius(out, radius);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('.', 'S') => {
            output_shape(out, "Cone", xyz02);
            output_rot(out, rot, angle);
            output_scale(out, cone_scale);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('S', '.') => {
            output_shape(out, "Cone", xyz01);
            output_rot(out, rot, angle);
            output_scale(out, cone_scale);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('H', '.') => {
            output_shape(out, "Box", xyz01);
            output_rot(out, rot, angle);
            output_scale(out, box_scale);
            output_color(out, dc);
            output_tr(out, tr);
        }
        ('.', 'H') => {
            output_shape(out, "Box", xyz02);
            output_rot(out, rot, angle);
            output_scale(out, box_scale);
            output_color(out, dc);
            output_tr(out, tr);
        }
        _ => {}
    }
}

fn render_vrml_body(residues: &[Residue3d], base_pairs: &[BasePair]) -> String {
    let num_residue = residues.len() - 1; // 1-based

    let c4xyz: Vec<[f64; 3]> = residues.iter().map(|r| r.c4).collect();

    let mut chain_ranges: Vec<(usize, usize)> = Vec::new();
    if num_residue > 0 {
        let mut start = 1usize;
        for idx in 2..=num_residue {
            if residues[idx - 1].chain != residues[idx].chain {
                chain_ranges.push((start, idx - 1));
                start = idx;
            }
        }
        chain_ranges.push((start, num_residue));
    }

    let mut w = VrmlWriter::new();
    w.newline();
    w.node("Collision");
    w.s_newline("collide FALSE");
    w.s("children [");
    w.add_indent(2);
    w.newline();

    for (i, _) in chain_ranges.iter().enumerate() {
        let chain_idx = i + 1;
        vrml_draw_chain(&mut w, chain_idx, &chain_ranges, &c4xyz);
        // label_residue writes directly to string with legacy formatting.
        label_residue(&mut w.out, residues, chain_ranges[i].0, chain_ranges[i].1);
    }

    let mut base_pairs_in_out_order: Vec<&BasePair> = base_pairs.iter().collect();
    base_pairs_in_out_order.sort_by_key(|bp| bp.out_index.unwrap_or(i32::MAX));
    for bp in base_pairs_in_out_order {
        if bp.kind != "pair" {
            continue;
        }
        if bp.note.as_deref().is_some_and(|s| s.contains('!')) {
            continue;
        }
        let Some(pt) = pair_type3_from_pairs(bp) else {
            continue;
        };
        let k1 = bp.i as usize;
        let k2 = bp.j as usize;
        if k1 == 0 || k2 == 0 || k1 > num_residue || k2 > num_residue {
            continue;
        }
        draw_interact(&mut w.out, k1, k2, pt, &c4xyz);
    }

    w.finish_list();
    w.finish_node();
    w.newline();
    w.node("NavigationInfo");
    w.s("speed");
    w.s("4");
    w.newline();
    w.list("type");
    w.qs("EXAMINE");
    w.qs("FLY");
    w.finish_list();
    w.finish_node();
    w.newline();
    w.finish()
}

pub fn render_vrml_from_pairs_json(
    pairs: &PairsJson,
    structure_path: &Path,
    structure_policies: &StructurePolicies,
    basepars_vrml_path: &Path,
) -> Result<String, String> {
    let residues_atoms = parse_structure_bases_with_atoms_with_policies(structure_path, structure_policies)
        .map_err(|e| e.to_string())?;
    if residues_atoms.is_empty() {
        return Err("structure has no nucleic residues".to_string());
    }

    let mut plane_points: Vec<[f64; 3]> = Vec::new();
    for r in &residues_atoms {
        for a in &r.atoms {
            if matches!(a.name.as_str(), "P" | "O5'" | "C5'" | "C4'" | "C3'" | "O3'") {
                plane_points.push([a.x, a.y, a.z]);
            }
        }
    }
    let rotmat = rot_2_lsplane_matrix(&plane_points);

    let mut residues: Vec<Residue3d> = Vec::with_capacity(residues_atoms.len() + 1);
    residues.push(Residue3d {
        chain: ' ',
        resname_field: [' '; 3],
        c4: [0.0; 3],
    });

    let mut avg = [0.0f64; 3];
    for r in residues_atoms {
        let mut c4 = None;
        for a in &r.atoms {
            if a.name == "C4'" {
                c4 = Some([a.x, a.y, a.z]);
                break;
            }
        }
        let c4 = rotate_point(rotmat, c4.unwrap_or([0.0, 0.0, 0.0]));
        avg[0] += c4[0];
        avg[1] += c4[1];
        avg[2] += c4[2];
        residues.push(Residue3d {
            chain: r.chain,
            resname_field: residue_name_field_legacy(&r.resname),
            c4,
        });
    }

    let n = (residues.len() - 1) as f64;
    if n > 0.0 {
        avg[0] /= n;
        avg[1] /= n;
        avg[2] /= n;
    }
    for r in residues.iter_mut().skip(1) {
        r.c4[0] -= avg[0];
        r.c4[1] -= avg[1];
        r.c4[2] -= avg[2];
        r.c4[2] -= 50.0;
    }

    let mut out = String::new();
    out.push_str("#VRML V2.0 utf8\n\n");
    let par = fs::read_to_string(basepars_vrml_path).map_err(|e| e.to_string())?;
    out.push_str(&par);

    // Legacy does not update wrapping state while copying the .par file; start fresh for the dynamic body.
    let body = render_vrml_body(&residues, &pairs.core.base_pairs);
    out.push_str(&body);
    Ok(out)
}
