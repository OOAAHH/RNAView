use std::path::{Path, PathBuf};

const XEPS: f64 = 1.0e-7;
const PI: f64 = std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct AtomName4(pub(crate) [u8; 4]);

impl AtomName4 {
    pub(crate) const fn from_bytes(b: [u8; 4]) -> Self {
        Self(b)
    }

    pub(crate) fn from_str4(s: &str) -> Result<Self, String> {
        let b = s.as_bytes();
        if b.len() != 4 {
            return Err(format!("AtomName4 must be 4 bytes, got {s:?}"));
        }
        Ok(Self([b[0], b[1], b[2], b[3]]))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ResName3(pub(crate) [u8; 3]);

impl ResName3 {
    pub(crate) const fn from_bytes(b: [u8; 3]) -> Self {
        Self(b)
    }

    pub(crate) fn from_str3(s: &str) -> Result<Self, String> {
        let b = s.as_bytes();
        if b.len() != 3 {
            return Err(format!("ResName3 must be 3 bytes, got {s:?}"));
        }
        Ok(Self([b[0], b[1], b[2]]))
    }
}

pub(crate) fn dot(va: &[f64; 4], vb: &[f64; 4]) -> f64 {
    va[1] * vb[1] + va[2] * vb[2] + va[3] * vb[3]
}

pub(crate) fn cross(va: &[f64; 4], vb: &[f64; 4]) -> [f64; 4] {
    [
        0.0,
        va[2] * vb[3] - va[3] * vb[2],
        va[3] * vb[1] - va[1] * vb[3],
        va[1] * vb[2] - va[2] * vb[1],
    ]
}

pub(crate) fn veclen(va: &[f64; 4]) -> f64 {
    dot(va, va).sqrt()
}

pub(crate) fn vec_norm(va: &mut [f64; 4]) {
    let vlen = veclen(va);
    if vlen > XEPS {
        va[1] /= vlen;
        va[2] /= vlen;
        va[3] /= vlen;
    }
}

pub(crate) fn rad2deg(ang: f64) -> f64 {
    ang * 180.0 / PI
}

pub(crate) fn dot2ang(dotval: f64) -> f64 {
    if dotval >= 1.0 {
        0.0
    } else if dotval <= -1.0 {
        180.0
    } else {
        rad2deg(dotval.acos())
    }
}

pub(crate) fn dswap(a: &mut f64, b: &mut f64) {
    let t = *a;
    *a = *b;
    *b = t;
}

fn identity_matrix(n: usize) -> Vec<Vec<f64>> {
    let mut out = vec![vec![0.0; n + 1]; n + 1];
    for i in 1..=n {
        out[i][i] = 1.0;
    }
    out
}

fn rotate(a: &mut [Vec<f64>], i: usize, j: usize, k: usize, l: usize, s: f64, tau: f64) {
    let g = a[i][j];
    let h = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);
}

fn eigsrt(d: &mut [f64], v: &mut [Vec<f64>], n: usize) {
    for i in 1..n {
        let mut k = i;
        let mut p = d[i];
        for j in (i + 1)..=n {
            if d[j] < p {
                p = d[j];
                k = j;
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for j in 1..=n {
                let tmp = v[j][i];
                v[j][i] = v[j][k];
                v[j][k] = tmp;
            }
        }
    }
}

fn jacobi(a: &mut [Vec<f64>], n: usize, d: &mut [f64], v: &mut [Vec<f64>]) -> Result<(), String> {
    if n == 0 {
        return Err("jacobi: n=0".to_string());
    }
    if d.len() < n + 1 {
        return Err("jacobi: d too small".to_string());
    }
    if v.len() < n + 1 || v.iter().any(|row| row.len() < n + 1) {
        return Err("jacobi: v too small".to_string());
    }

    let mut b = vec![0.0; n + 1];
    let mut z = vec![0.0; n + 1];
    let ident = identity_matrix(n);
    for i in 1..=n {
        for j in 1..=n {
            v[i][j] = ident[i][j];
        }
    }
    for ip in 1..=n {
        b[ip] = a[ip][ip];
        d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }

    for iter in 1..=100 {
        let mut sm = 0.0;
        for ip in 1..=(n - 1) {
            for iq in (ip + 1)..=n {
                sm += a[ip][iq].abs();
            }
        }
        if sm < XEPS {
            eigsrt(d, v, n);
            return Ok(());
        }
        let tresh = if iter < 4 { 0.2 * sm / ((n * n) as f64) } else { 0.0 };

        for ip in 1..=(n - 1) {
            for iq in (ip + 1)..=n {
                let g = 100.0 * a[ip][iq].abs();
                if iter > 4 && (d[ip].abs() + g) == d[ip].abs() && (d[iq].abs() + g) == d[iq].abs() {
                    a[ip][iq] = 0.0;
                    continue;
                }
                if a[ip][iq].abs() <= tresh {
                    continue;
                }
                let mut h = d[iq] - d[ip];
                let t = if (h.abs() + g) == h.abs() {
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
                h = t * a[ip][iq];
                z[ip] -= h;
                z[iq] += h;
                d[ip] -= h;
                d[iq] += h;
                a[ip][iq] = 0.0;
                for j in 1..=(ip - 1) {
                    rotate(a, j, ip, j, iq, s, tau);
                }
                for j in (ip + 1)..=(iq - 1) {
                    rotate(a, ip, j, j, iq, s, tau);
                }
                for j in (iq + 1)..=n {
                    rotate(a, ip, j, iq, j, s, tau);
                }
                for j in 1..=n {
                    rotate(v, j, ip, j, iq, s, tau);
                }
            }
        }

        for ip in 1..=n {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }

    Err("jacobi: too many iterations".to_string())
}

fn ave_dmatrix(a: &[Vec<f64>], nr: usize, nc: usize) -> Result<Vec<f64>, String> {
    if nr == 0 {
        return Err("ave_dmatrix: nr=0".to_string());
    }
    let mut out = vec![0.0; nc + 1];
    for i in 1..=nc {
        let mut sum = 0.0;
        for j in 1..=nr {
            sum += a[j][i];
        }
        out[i] = sum / (nr as f64);
    }
    Ok(out)
}

fn cov_matrix(a: &[Vec<f64>], b: &[Vec<f64>], nr: usize, nc: usize) -> Result<Vec<Vec<f64>>, String> {
    if nr < 2 {
        return Err("cov_matrix: nr<2".to_string());
    }
    let ave_a = ave_dmatrix(a, nr, nc)?;
    let ave_b = ave_dmatrix(b, nr, nc)?;

    let mut ta_x_b = vec![vec![0.0; nc + 1]; nc + 1];
    for i in 1..=nc {
        for j in 1..=nc {
            let mut sum = 0.0;
            for k in 1..=nr {
                sum += a[k][i] * b[k][j];
            }
            ta_x_b[i][j] = sum;
        }
    }

    let mut cmtx = vec![vec![0.0; nc + 1]; nc + 1];
    for i in 1..=nc {
        for j in 1..=nc {
            cmtx[i][j] = (ta_x_b[i][j] - ave_a[i] * ave_b[j] * (nr as f64)) / ((nr - 1) as f64);
        }
    }
    Ok(cmtx)
}

pub(crate) fn ls_fitting(
    sxyz: &[Vec<f64>],
    exyz: &[Vec<f64>],
    n: usize,
) -> Result<([f64; 4], [[f64; 4]; 4]), String> {
    if n < 3 {
        return Err("too few atoms for least-squares fitting".to_string());
    }
    let u = cov_matrix(sxyz, exyz, n, 3)?;

    let mut nmat = vec![vec![0.0; 5]; 5];
    nmat[1][1] = u[1][1] + u[2][2] + u[3][3];
    nmat[2][2] = u[1][1] - u[2][2] - u[3][3];
    nmat[3][3] = -u[1][1] + u[2][2] - u[3][3];
    nmat[4][4] = -u[1][1] - u[2][2] + u[3][3];
    nmat[1][2] = u[2][3] - u[3][2];
    nmat[2][1] = nmat[1][2];
    nmat[1][3] = u[3][1] - u[1][3];
    nmat[3][1] = nmat[1][3];
    nmat[1][4] = u[1][2] - u[2][1];
    nmat[4][1] = nmat[1][4];
    nmat[2][3] = u[1][2] + u[2][1];
    nmat[3][2] = nmat[2][3];
    nmat[2][4] = u[3][1] + u[1][3];
    nmat[4][2] = nmat[2][4];
    nmat[3][4] = u[2][3] + u[3][2];
    nmat[4][3] = nmat[3][4];

    let mut v = vec![vec![0.0; 5]; 5];
    let mut d = vec![0.0; 5];
    jacobi(&mut nmat, 4, &mut d, &mut v)?;

    for i in 1..=4 {
        for j in 1..=4 {
            nmat[i][j] = v[i][4] * v[j][4];
        }
    }

    let mut r = [[0.0; 4]; 4];
    r[1][1] = nmat[1][1] + nmat[2][2] - nmat[3][3] - nmat[4][4];
    r[1][2] = 2.0 * (nmat[2][3] - nmat[1][4]);
    r[1][3] = 2.0 * (nmat[2][4] + nmat[1][3]);
    r[2][1] = 2.0 * (nmat[3][2] + nmat[1][4]);
    r[2][2] = nmat[1][1] - nmat[2][2] + nmat[3][3] - nmat[4][4];
    r[2][3] = 2.0 * (nmat[3][4] - nmat[1][2]);
    r[3][1] = 2.0 * (nmat[4][2] - nmat[1][3]);
    r[3][2] = 2.0 * (nmat[4][3] + nmat[1][2]);
    r[3][3] = nmat[1][1] - nmat[2][2] - nmat[3][3] + nmat[4][4];

    let ave_sxyz = ave_dmatrix(sxyz, n, 3)?;
    let ave_exyz = ave_dmatrix(exyz, n, 3)?;
    let ave_s = [0.0, ave_sxyz[1], ave_sxyz[2], ave_sxyz[3]];

    let mut orgi = [0.0; 4];
    for i in 1..=3 {
        let row = [0.0, r[i][1], r[i][2], r[i][3]];
        orgi[i] = ave_exyz[i] - dot(&ave_s, &row);
    }
    Ok((orgi, r))
}

pub(crate) fn get_bdir(filename: &str) -> Result<PathBuf, String> {
    if Path::new(filename).exists() {
        return Ok(PathBuf::from("./"));
    }

    if let Some(root) = std::env::var_os("RNAVIEW") {
        let mut p = PathBuf::from(root);
        p.push("BASEPARS");
        p.push("");
        return Ok(p);
    }

    if let Some(home) = std::env::var_os("HOME").or_else(|| std::env::var_os("HOMEDRIVE")) {
        let mut p = PathBuf::from(home);
        p.push("RNAVIEW");
        p.push("BASEPARS");
        p.push("");
        return Ok(p);
    }

    Err("cannot locate base geometry and parameter files (RNAVIEW env not set)".to_string())
}

#[derive(Debug, Clone)]
pub(crate) struct StdBase {
    pub(crate) atom_names: Vec<AtomName4>, // 1-based
    pub(crate) xyz: Vec<[f64; 4]>,         // 1-based
}

fn read_pdb_ref(path: &Path) -> Result<StdBase, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("read_pdb_ref failed for {path:?}: {e}"))?;

    let mut atom_names: Vec<AtomName4> = vec![AtomName4([0; 4])];
    let mut xyz: Vec<[f64; 4]> = vec![[0.0; 4]];

    for (line_no, line) in text.lines().enumerate() {
        if !line.starts_with("ATOM") {
            continue;
        }
        let b = line.as_bytes();
        let name = b
            .get(12..16)
            .ok_or_else(|| format!("read_pdb_ref: short ATOM line {}: {line:?}", line_no + 1))?;
        let name4 = AtomName4::from_bytes([name[0], name[1], name[2], name[3]]);
        let (x, y, z) = if let Some(coord) = line.get(30..) {
            let mut it = coord.split_whitespace();
            let x = it
                .next()
                .ok_or_else(|| format!("read_pdb_ref: missing x {}: {line:?}", line_no + 1))?
                .parse::<f64>()
                .map_err(|_| format!("read_pdb_ref: invalid x {}: {line:?}", line_no + 1))?;
            let y = it
                .next()
                .ok_or_else(|| format!("read_pdb_ref: missing y {}: {line:?}", line_no + 1))?
                .parse::<f64>()
                .map_err(|_| format!("read_pdb_ref: invalid y {}: {line:?}", line_no + 1))?;
            let z = it
                .next()
                .ok_or_else(|| format!("read_pdb_ref: missing z {}: {line:?}", line_no + 1))?
                .parse::<f64>()
                .map_err(|_| format!("read_pdb_ref: invalid z {}: {line:?}", line_no + 1))?;
            (x, y, z)
        } else {
            let toks: Vec<&str> = line.split_whitespace().collect();
            if toks.len() < 3 {
                return Err(format!(
                    "read_pdb_ref: missing xyz fields {}: {line:?}",
                    line_no + 1
                ));
            }
            let x = toks[toks.len() - 3]
                .parse::<f64>()
                .map_err(|_| format!("read_pdb_ref: invalid x {}: {line:?}", line_no + 1))?;
            let y = toks[toks.len() - 2]
                .parse::<f64>()
                .map_err(|_| format!("read_pdb_ref: invalid y {}: {line:?}", line_no + 1))?;
            let z = toks[toks.len() - 1]
                .parse::<f64>()
                .map_err(|_| format!("read_pdb_ref: invalid z {}: {line:?}", line_no + 1))?;
            (x, y, z)
        };
        atom_names.push(name4);
        xyz.push([0.0, x, y, z]);
    }

    Ok(StdBase { atom_names, xyz })
}

pub(crate) fn load_reference_bases(bdir: &Path) -> Result<[StdBase; 7], String> {
    let refs = ['A', 'G', 'C', 'U', 'T', 'I', 'P'];
    let mut out: Vec<StdBase> = Vec::with_capacity(refs.len());
    for base in refs {
        let path = bdir.join(format!("Atomic_{base}.pdb"));
        out.push(read_pdb_ref(&path)?);
    }
    out.try_into()
        .map_err(|_| "unexpected reference base count".to_string())
}

pub(crate) fn ref_idx(resnam: u8) -> Option<usize> {
    match resnam {
        b'A' | b'a' => Some(0),
        b'G' | b'g' => Some(1),
        b'C' | b'c' => Some(2),
        b'U' | b'u' => Some(3),
        b'T' | b't' => Some(4),
        b'I' => Some(5),
        b'P' => Some(6),
        _ => None,
    }
}

pub(crate) fn find_1st_atom(
    needle: AtomName4,
    haystack: &[AtomName4],
    nb: usize,
    ne: usize,
) -> Option<usize> {
    if nb == 0 || ne == 0 || nb > ne || ne >= haystack.len() {
        return None;
    }
    for i in nb..=ne {
        if haystack[i] == needle {
            return Some(i);
        }
    }
    None
}

pub(crate) fn read_bprs_from_misc_rna_par() -> Result<[f64; 7], String> {
    let bdir = get_bdir("Atomic_A.pdb")?;
    let path = bdir.join("misc_rna.par");
    let text = std::fs::read_to_string(&path)
        .map_err(|e| format!("read misc_rna.par failed for {path:?}: {e}"))?;
    let mut bprs = [0.0f64; 7];
    let mut i = 1usize;
    for line in text.lines() {
        if i > 6 {
            break;
        }
        let tok = line.split_whitespace().next().unwrap_or("");
        if tok.is_empty() {
            continue;
        }
        bprs[i] = tok
            .parse::<f64>()
            .map_err(|_| format!("invalid misc_rna.par value at line {i}: {line:?}"))?;
        i += 1;
    }
    if i <= 6 {
        return Err(format!("misc_rna.par too short: expected 6 values, got {}", i - 1));
    }
    bprs[4] = 90.0 - (bprs[4].abs() - 90.0).abs();
    if bprs[6] > 12.0 {
        bprs[6] = 12.0;
    }
    Ok(bprs)
}

#[derive(Debug, Clone)]
pub(crate) struct BaseInfo {
    pub(crate) bprs: [f64; 7],
    pub(crate) orien: Vec<[f64; 10]>,
    pub(crate) org: Vec<[f64; 4]>,
    pub(crate) nxyz: Vec<[f64; 4]>,
}

pub(crate) fn base_frame(
    num_residue: usize,
    bseq: &[u8],
    seidx: &[[usize; 3]],
    ry: &[i32],
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
    std_bases: &[StdBase; 7],
) -> Result<(Vec<[f64; 10]>, Vec<[f64; 4]>), String> {
    const RING_ATOMS: [AtomName4; 9] = [
        AtomName4::from_bytes(*b" C4 "),
        AtomName4::from_bytes(*b" N3 "),
        AtomName4::from_bytes(*b" C2 "),
        AtomName4::from_bytes(*b" N1 "),
        AtomName4::from_bytes(*b" C6 "),
        AtomName4::from_bytes(*b" C5 "),
        AtomName4::from_bytes(*b" N7 "),
        AtomName4::from_bytes(*b" C8 "),
        AtomName4::from_bytes(*b" N9 "),
    ];

    let mut orien: Vec<[f64; 10]> = vec![[0.0; 10]; num_residue + 1];
    let mut org: Vec<[f64; 4]> = vec![[0.0; 4]; num_residue + 1];

    for i in 1..=num_residue {
        if ry.get(i).copied().unwrap_or(-1) < 0 {
            continue;
        }
        let resnam = *bseq.get(i).unwrap_or(&b'?');
        let Some(ii) = ref_idx(resnam) else {
            continue;
        };
        let std = &std_bases[ii];
        let std_natom = std.atom_names.len().saturating_sub(1);

        let ib = seidx
            .get(i)
            .and_then(|v| v.get(1).copied())
            .unwrap_or(0);
        let ie = seidx
            .get(i)
            .and_then(|v| v.get(2).copied())
            .unwrap_or(0);
        if ib == 0 || ie == 0 || ib > ie {
            continue;
        }

        let ring_atom_num = if ry[i] == 1 { 9usize } else { 6usize };
        let mut nmatch = 0usize;
        let mut e_ring: Vec<Vec<f64>> = vec![vec![0.0; 4]; ring_atom_num + 1];
        let mut s_ring: Vec<Vec<f64>> = vec![vec![0.0; 4]; ring_atom_num + 1];

        for atom in RING_ATOMS.iter().take(ring_atom_num) {
            let exp_katom = find_1st_atom(*atom, atom_name, ib, ie);
            let std_katom = find_1st_atom(*atom, &std.atom_names, 1, std_natom);
            if let (Some(exp), Some(st)) = (exp_katom, std_katom) {
                nmatch += 1;
                for k in 1..=3 {
                    e_ring[nmatch][k] = xyz[exp][k];
                    s_ring[nmatch][k] = std.xyz[st][k];
                }
            }
        }

        if nmatch < 3 {
            return Err(format!("base_frame: too few matched ring atoms for residue {i}"));
        }

        let (orgi, r) = ls_fitting(&s_ring, &e_ring, nmatch)?;
        for j in 1..=3 {
            org[i][j] = orgi[j];
            orien[i][j] = r[j][1];
            orien[i][j + 3] = r[j][2];
            orien[i][j + 6] = r[j][3];
        }
    }

    Ok((orien, org))
}

pub(crate) fn compute_base_info(
    num_residue: usize,
    bseq: &[u8],
    seidx: &[[usize; 3]],
    ry: &[i32],
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
) -> Result<BaseInfo, String> {
    let bprs = read_bprs_from_misc_rna_par()?;

    let bdir = get_bdir("Atomic_A.pdb")?;
    let std_bases = load_reference_bases(&bdir)?;

    let (orien, org) = base_frame(num_residue, bseq, seidx, ry, atom_name, xyz, &std_bases)?;

    let mut nxyz: Vec<[f64; 4]> = vec![[0.0; 4]; num_residue + 1];
    let n9 = AtomName4::from_bytes(*b" N9 ");
    let n1 = AtomName4::from_bytes(*b" N1 ");
    let c5 = AtomName4::from_bytes(*b" C5 ");

    for i in 1..=num_residue {
        let ib = seidx
            .get(i)
            .and_then(|v| v.get(1).copied())
            .unwrap_or(0);
        let ie = seidx
            .get(i)
            .and_then(|v| v.get(2).copied())
            .unwrap_or(0);
        if ib == 0 || ie == 0 || ib > ie {
            continue;
        }
        if ry.get(i).copied().unwrap_or(-1) < 0 {
            continue;
        }

        let target = if ry[i] == 1 {
            n9
        } else if bseq[i] == b'P' || bseq[i] == b'p' {
            c5
        } else {
            n1
        };
        let Some(idx) = find_1st_atom(target, atom_name, ib, ie) else {
            return Err(format!("compute_base_info: missing N9/N1 anchor atom for residue {i}"));
        };
        nxyz[i][1] = xyz[idx][1];
        nxyz[i][2] = xyz[idx][2];
        nxyz[i][3] = xyz[idx][3];
    }

    Ok(BaseInfo { bprs, orien, org, nxyz })
}
