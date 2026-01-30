use crate::legacy_alg::{dot2ang, find_1st_atom, AtomName4, ResName3};
use crate::out_full::{OutBasePairLine, OutPairKind};

const BUF512: usize = 512;
const NP: usize = 101;
const MFACTOR: f64 = 10000.0;
const BOND_UPPER_LIMIT: f64 = 2.5;
const EMPTY_NUMBER: f64 = -9999.99;
const WC_DORG: f64 = 2.5;
const XEPS: f64 = 1.0e-7;

fn is_ascii_digit(b: u8) -> bool {
    (b'0'..=b'9').contains(&b)
}

fn to_ascii_upper(b: u8) -> u8 {
    if (b'a'..=b'z').contains(&b) {
        b - 32
    } else {
        b
    }
}

fn is_in(set: &[u8], b: u8) -> bool {
    set.contains(&b)
}

fn dot3(a: &[f64; 4], b: &[f64; 4]) -> f64 {
    a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
}

fn cross3(a: &[f64; 4], b: &[f64; 4]) -> [f64; 4] {
    [
        0.0,
        a[2] * b[3] - a[3] * b[2],
        a[3] * b[1] - a[1] * b[3],
        a[1] * b[2] - a[2] * b[1],
    ]
}

fn veclen3(v: &[f64; 4]) -> f64 {
    dot3(v, v).sqrt()
}

fn vec_norm3(v: &mut [f64; 4]) {
    let len = veclen3(v);
    if len > XEPS {
        v[1] /= len;
        v[2] /= len;
        v[3] /= len;
    }
}

fn magang(mut a: [f64; 4], mut b: [f64; 4]) -> f64 {
    if veclen3(&a) < XEPS || veclen3(&b) < XEPS {
        return 0.0;
    }
    vec_norm3(&mut a);
    vec_norm3(&mut b);
    dot2ang(dot3(&a, &b))
}

fn vec_orth(v: &mut [f64; 4], mut vref: [f64; 4]) {
    vec_norm3(&mut vref);
    let d = dot3(v, &vref);
    for i in 1..=3 {
        v[i] -= d * vref[i];
    }
    vec_norm3(v);
}

fn vec_ang(va: [f64; 4], vb: [f64; 4], vref: [f64; 4]) -> f64 {
    let mut a = va;
    let mut b = vb;
    vec_orth(&mut a, vref);
    vec_orth(&mut b, vref);
    let mut ang = magang(a, b);
    let vc = cross3(&a, &b);
    if dot3(&vc, &vref) < 0.0 {
        ang = -ang;
    }
    ang
}

fn torsion(d: &[[f64; 4]; 5]) -> f64 {
    let mut vec3 = [[0.0f64; 4]; 4];
    for i in 1..=3 {
        for j in 1..=3 {
            vec3[i][j] = if i == 1 {
                d[i][j] - d[i + 1][j]
            } else {
                d[i + 1][j] - d[i][j]
            };
        }
        let dij = veclen3(&vec3[i]);
        if dij > BOND_UPPER_LIMIT {
            return EMPTY_NUMBER;
        }
    }
    vec_ang(vec3[1], vec3[3], vec3[2])
}

fn as_chars3(bytes: [u8; 3]) -> String {
    String::from_utf8_lossy(&bytes).to_string()
}

fn as_chars4(bytes: [u8; 4]) -> String {
    String::from_utf8_lossy(&bytes).to_string()
}

fn resname_eq(a: ResName3, b: [u8; 3]) -> bool {
    a.0 == b
}

fn atom_eq(a: AtomName4, b: [u8; 4]) -> bool {
    a.0 == b
}

fn pair_code_from_type_field(type_field: &str) -> [u8; 3] {
    let b = type_field.as_bytes();
    if b.len() < 5 {
        return [b'?', b'?', b'?'];
    }
    [b[0], b[2], b[4]]
}

fn ring_center(
    i: usize,
    seidx: &[[usize; 3]],
    bseq: &[u8],
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
) -> [f64; 4] {
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

    let mut xyz_c = [0.0f64; 4];
    let base = *bseq.get(i).unwrap_or(&b'?') as char;
    let natm = if matches!(base, 'A' | 'G' | 'a' | 'g' | 'I' | 'i') {
        9usize
    } else {
        6usize
    };

    let ib = seidx.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0);
    let ie = seidx.get(i).and_then(|v| v.get(2).copied()).unwrap_or(0);
    if ib == 0 || ie == 0 || ib > ie {
        return xyz_c;
    }

    let mut n = 0i64;
    for atom in RING_ATOMS.iter().take(natm) {
        if let Some(j) = find_1st_atom(*atom, atom_name, ib, ie) {
            n += 1;
            for k in 1..=3 {
                xyz_c[k] += xyz[j][k];
            }
        }
    }
    if n > 0 {
        for k in 1..=3 {
            xyz_c[k] /= n as f64;
        }
    }
    xyz_c
}

fn base_stack(
    i: usize,
    j: usize,
    bseq: &[u8],
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
    rtn_val: &[f64; 21],
) -> bool {
    let xyz_c1 = ring_center(i, seidx, bseq, atom_name, xyz);
    let xyz_c2 = ring_center(j, seidx, bseq, atom_name, xyz);
    let mut cc_vec = [0.0f64; 4];
    for k in 1..=3 {
        cc_vec[k] = xyz_c2[k] - xyz_c1[k];
    }
    let cc_dist = veclen3(&cc_vec);

    let bi = to_ascii_upper(*bseq.get(i).unwrap_or(&b'?'));
    let bj = to_ascii_upper(*bseq.get(j).unwrap_or(&b'?'));
    let is_pyr_i = matches!(bi, b'C' | b'U' | b'T');
    let is_pyr_j = matches!(bj, b'C' | b'U' | b'T');

    let dv = rtn_val[2];
    if is_pyr_i && is_pyr_j {
        if dv > 2.3 {
            cc_dist < 5.6
        } else if dv <= 2.3 && dv >= 1.9 {
            cc_dist < 5.5
        } else {
            cc_dist < 5.0
        }
    } else {
        if dv > 2.3 {
            cc_dist < 5.8
        } else if dv <= 2.3 && dv >= 1.9 {
            cc_dist < 5.7
        } else {
            cc_dist < 5.4
        }
    }
}

fn H_catalog(i: usize, m: usize, bseq: &[u8], atom_name: &[AtomName4]) -> (i64, i64) {
    let mut without_h = 0i64;
    let mut with_h = 0i64;

    let base = to_ascii_upper(*bseq.get(i).unwrap_or(&b'?'));
    let a = atom_name.get(m).copied().unwrap_or(AtomName4([0; 4]));

    let eq = |s: [u8; 4]| atom_eq(a, s);

    match base {
        b'A' => {
            if eq(*b" O3'") || eq(*b" O4'") || eq(*b" O5'") || eq(*b" O1P") || eq(*b" O2P") || eq(*b" N9 ")
                || eq(*b" N7 ") || eq(*b" N3 ")
            {
                without_h = 1;
            }
            if eq(*b" N1 ") || eq(*b" N6 ") || eq(*b" C8 ") || eq(*b" C2 ") || eq(*b" O2'") {
                with_h = 1;
            }
        }
        b'G' => {
            if eq(*b" O3'") || eq(*b" O4'") || eq(*b" O5'") || eq(*b" O1P") || eq(*b" O2P") || eq(*b" N9 ")
                || eq(*b" N7 ") || eq(*b" O6 ") || eq(*b" N3 ")
            {
                without_h = 1;
            }
            if eq(*b" N1 ") || eq(*b" N2 ") || eq(*b" C8 ") || eq(*b" O2'") {
                with_h = 1;
            }
        }
        b'I' => {
            if eq(*b" O3'") || eq(*b" O4'") || eq(*b" O5'") || eq(*b" O1P") || eq(*b" O2P") || eq(*b" N9 ")
                || eq(*b" N7 ") || eq(*b" O6 ") || eq(*b" N3 ")
            {
                without_h = 1;
            }
            if eq(*b" N1 ") || eq(*b" C2 ") || eq(*b" C8 ") || eq(*b" O2'") {
                with_h = 1;
            }
        }
        b'U' => {
            if eq(*b" O3'") || eq(*b" O4'") || eq(*b" O5'") || eq(*b" O1P") || eq(*b" O2P") || eq(*b" O4 ")
                || eq(*b" O2 ") || eq(*b" N1 ")
            {
                without_h = 1;
            }
            if eq(*b" N3 ") || eq(*b" C5 ") || eq(*b" C6 ") || eq(*b" O2'") {
                with_h = 1;
            }
        }
        b'T' => {
            if eq(*b" O3'") || eq(*b" O4'") || eq(*b" O5'") || eq(*b" O1P") || eq(*b" O2P") || eq(*b" O2 ")
                || eq(*b" O4 ") || eq(*b" N1 ")
            {
                without_h = 1;
            }
            if eq(*b" N3 ") || eq(*b" C5 ") || eq(*b" C6 ") {
                with_h = 1;
            }
        }
        b'C' => {
            if eq(*b" O3'") || eq(*b" O4'") || eq(*b" O5'") || eq(*b" O1P") || eq(*b" O2P") || eq(*b" N1 ")
                || eq(*b" O2 ")
            {
                without_h = 1;
            }
            if eq(*b" N3 ") || eq(*b" N4 ") || eq(*b" C5 ") || eq(*b" C6 ") || eq(*b" O2'") {
                with_h = 1;
            }
        }
        b'P' => {
            if eq(*b" O3'") || eq(*b" O4'") || eq(*b" O5'") || eq(*b" O1P") || eq(*b" O2P") || eq(*b" C5 ")
                || eq(*b" O4 ") || eq(*b" O2 ")
            {
                without_h = 1;
            }
            if eq(*b" N1 ") || eq(*b" N3 ") || eq(*b" C6 ") || eq(*b" O2'") {
                with_h = 1;
            }
        }
        _ => {}
    }

    (without_h, with_h)
}

fn Hbond_pair(
    i: usize,
    j: usize,
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    bseq: &[u8],
    xyz: &[[f64; 4]],
    change: f64,
    c_key: i64,
    bone_key: i64,
) -> Vec<(AtomName4, AtomName4, f64)> {
    let ib = seidx.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0);
    let ie = seidx.get(i).and_then(|v| v.get(2).copied()).unwrap_or(0);
    let jb = seidx.get(j).and_then(|v| v.get(1).copied()).unwrap_or(0);
    let je = seidx.get(j).and_then(|v| v.get(2).copied()).unwrap_or(0);
    if ib == 0 || ie == 0 || jb == 0 || je == 0 {
        return Vec::new();
    }

    let mut out: Vec<(AtomName4, AtomName4, f64)> = Vec::new();

    let bi = to_ascii_upper(*bseq.get(i).unwrap_or(&b'?'));
    let bj = to_ascii_upper(*bseq.get(j).unwrap_or(&b'?'));

    for m in ib..=ie {
        let am = atom_name[m];
        if c_key == 0 && am.0[1] == b'C' {
            continue;
        }
        if bone_key == 0 && (atom_eq(am, *b" O3'") || atom_eq(am, *b" O2P") || atom_eq(am, *b" O5'") || atom_eq(am, *b" O1P"))
        {
            continue;
        }
        if (am.0[1] == b'C' && am.0[3] == b'\'') || am.0[1] == b'P' {
            continue;
        }
        if matches!(bi, b'A' | b'I') {
            if atom_eq(am, *b" C4 ") || atom_eq(am, *b" C5 ") || atom_eq(am, *b" C6 ") {
                continue;
            }
        } else if bi == b'G' {
            if atom_eq(am, *b" C4 ") || atom_eq(am, *b" C5 ") || atom_eq(am, *b" C6 ") || atom_eq(am, *b" C2 ") {
                continue;
            }
        } else if bi == b'P' {
            if atom_eq(am, *b" C4 ") || atom_eq(am, *b" C5 ") {
                continue;
            }
        } else if matches!(bi, b'U' | b'C' | b'T') {
            if atom_eq(am, *b" C4 ") || atom_eq(am, *b" C2 ") {
                continue;
            }
        }

        let (without_h_m, _with_h_m) = H_catalog(i, m, bseq, atom_name);

        for n in jb..=je {
            let an = atom_name[n];
            if c_key == 0 && an.0[1] == b'C' {
                continue;
            }
            if bone_key == 0
                && (atom_eq(an, *b" O3'")
                    || atom_eq(an, *b" O2P")
                    || atom_eq(an, *b" O5'")
                    || atom_eq(an, *b" O1P"))
            {
                continue;
            }
            if (an.0[1] == b'C' && an.0[3] == b'\'') || an.0[1] == b'P' {
                continue;
            }
            if matches!(bj, b'A' | b'I') {
                if atom_eq(an, *b" C4 ") || atom_eq(an, *b" C5 ") || atom_eq(an, *b" C6 ") {
                    continue;
                }
            } else if bj == b'G' {
                if atom_eq(an, *b" C4 ") || atom_eq(an, *b" C5 ") || atom_eq(an, *b" C6 ") || atom_eq(an, *b" C2 ") {
                    continue;
                }
            } else if bj == b'P' {
                if atom_eq(an, *b" C4 ") || atom_eq(an, *b" C5 ") {
                    continue;
                }
            } else if matches!(bj, b'U' | b'C' | b'T') {
                if atom_eq(an, *b" C4 ") || atom_eq(an, *b" C2 ") {
                    continue;
                }
            }

            if am.0[1] == b'C' && an.0[1] == b'C' {
                continue;
            }

            let is_bone_o = |a: AtomName4| {
                atom_eq(a, *b" O3'")
                    || atom_eq(a, *b" O4'")
                    || atom_eq(a, *b" O5'")
                    || atom_eq(a, *b" O1P")
                    || atom_eq(a, *b" O2P")
            };
            if is_bone_o(am) && is_bone_o(an) {
                continue;
            }

            let (without_h_n, _with_h_n) = H_catalog(j, n, bseq, atom_name);
            if without_h_m == 1 && without_h_n == 1 {
                continue;
            }

            let m_base_no = is_in(b"NO", am.0[1]) && am.0[3] != b'\'' && am.0[3] != b'P';
            let n_base_no = is_in(b"NO", an.0[1]) && an.0[3] != b'\'' && an.0[3] != b'P';

            let dist = if m_base_no && n_base_no {
                (3.4 + change).min(4.0)
            } else if (am.0[1] == b'C' && n_base_no) || (an.0[1] == b'C' && m_base_no) {
                (3.6 + change).min(4.0)
            } else if (am.0[1] == b'O' && am.0[3] == b'\'' && n_base_no)
                || (an.0[1] == b'O' && an.0[3] == b'\'' && m_base_no)
            {
                (3.4 + change).min(4.0)
            } else if (am.0[3] == b'P' && an.0[3] != b'\'' && an.0[1] != b'C')
                || (an.0[3] == b'P' && am.0[3] != b'\'' && am.0[1] != b'C')
            {
                (3.2 + change).min(4.0)
            } else {
                (3.1 + change).min(3.8)
            };

            let mut dtmp = [0.0f64; 4];
            for k in 1..=3 {
                dtmp[k] = xyz[m][k] - xyz[n][k];
            }
            let dd = veclen3(&dtmp);
            if dd < dist {
                if out.len() >= BUF512 {
                    break;
                }
                out.push((am, an, dd));
            }
        }
    }

    out
}

pub(crate) fn check_pairs(
    i: usize,
    j: usize,
    bseq: &[u8],
    seidx: &[[usize; 3]],
    xyz: &[[f64; 4]],
    nxyz: &[[f64; 4]],
    orien: &[[f64; 10]],
    org: &[[f64; 4]],
    atom_name: &[AtomName4],
    bprs: &[f64; 7],
    rtn_val: &mut [f64; 21],
    network: i64,
) -> i64 {
    static WC: [[u8; 2]; 8] = [*b"AT", *b"AU", *b"TA", *b"UA", *b"GC", *b"CG", *b"IC", *b"CI"];

    if i == j {
        return 0;
    }

    let mut dorg = [0.0f64; 4];
    let mut dnn = [0.0f64; 4];
    for k in 1..=3 {
        dorg[k] = org[j][k] - org[i][k];
        dnn[k] = nxyz[j][k] - nxyz[i][k];
    }
    rtn_val[1] = veclen3(&dorg);

    let dot_orien = |a: usize, a_off: usize, b: usize, b_off: usize| -> f64 {
        orien[a][a_off + 1] * orien[b][b_off + 1]
            + orien[a][a_off + 2] * orien[b][b_off + 2]
            + orien[a][a_off + 3] * orien[b][b_off + 3]
    };

    let dir_x = dot_orien(i, 0, j, 0);
    let dir_y = dot_orien(i, 3, j, 3);
    let dd = dot_orien(i, 6, j, 6);

    let mut zave = [0.0f64; 4];
    if dd <= 0.0 {
        for k in 1..=3 {
            zave[k] = orien[i][6 + k] - orien[j][6 + k];
        }
    } else {
        for k in 1..=3 {
            zave[k] = orien[i][6 + k] + orien[j][6 + k];
        }
    }
    vec_norm3(&mut zave);

    rtn_val[2] = dot3(&dorg, &zave).abs();
    if rtn_val[2] > bprs[3] {
        return 0;
    }

    rtn_val[3] = 90.0 - (dot2ang(dd) - 90.0).abs();
    if rtn_val[3] > bprs[4] {
        return 0;
    }
    if rtn_val[3] <= 10.0 && rtn_val[2] >= 2.2 {
        return 0;
    }

    rtn_val[4] = veclen3(&dnn);
    if rtn_val[4] < bprs[5] {
        return 0;
    }

    rtn_val[5] = rtn_val[1] + 2.0 * rtn_val[2];

    if network != 0 {
        if rtn_val[2] <= bprs[3] && rtn_val[3] <= bprs[4] && rtn_val[4] >= bprs[5] {
            return 1;
        }
        return 0;
    }

    if j == i + 1 && rtn_val[2] >= 2.0 {
        return 0;
    }

    if rtn_val[1] > bprs[2] {
        return 0;
    }

    let ib = seidx[i][1];
    let ie = seidx[i][2];
    let jb = seidx[j][1];
    let je = seidx[j][2];

    let mut short_contact = false;
    for m in ib..=ie {
        if short_contact {
            break;
        }
        let am = atom_name[m].0;
        if !(is_in(b"ON", am[1]) && am[0] == b' ' && is_ascii_digit(am[2]) && am[3] == b' ') {
            continue;
        }
        for n in jb..=je {
            let an = atom_name[n].0;
            if !(is_in(b"ON", an[1]) && an[0] == b' ' && is_ascii_digit(an[2]) && an[3] == b' ') {
                continue;
            }
            let (without_h_m, _) = H_catalog(i, m, bseq, atom_name);
            let (without_h_n, _) = H_catalog(j, n, bseq, atom_name);
            if without_h_m == 1 && without_h_n == 1 {
                continue;
            }

            let mut tmp = [0.0f64; 4];
            for k in 1..=3 {
                tmp[k] = xyz[n][k] - xyz[m][k];
            }
            if veclen3(&tmp) <= bprs[1] {
                short_contact = true;
                break;
            }
        }
    }
    if !short_contact {
        return 0;
    }

    if base_stack(i, j, bseq, seidx, atom_name, xyz, rtn_val) {
        return 0;
    }

    let bi = to_ascii_upper(*bseq.get(i).unwrap_or(&b'?'));
    let bj = to_ascii_upper(*bseq.get(j).unwrap_or(&b'?'));
    let bpi = [bi, bj];

    let mut bpid: i64 = -1;
    if dir_x > 0.0 && dir_y < 0.0 && dd < 0.0 {
        bpid = 1;
        if rtn_val[1] <= WC_DORG
            && WC.iter().any(|p| p == &bpi)
            && rtn_val[3] <= 40.0
            && rtn_val[2] < 1.5
        {
            bpid = 2;
        }
    }
    if bpid == 2 {
        rtn_val[5] -= 1.5;
    }

    if bpid != 0 {
        // Base-pair origin (avg) and per-base normals, matching legacy check_pairs layout.
        for k in 1..=3 {
            rtn_val[5 + k] = 0.5 * (org[i][k] + org[j][k]);
            rtn_val[8 + k] = orien[i][6 + k];
            rtn_val[11 + k] = orien[j][6 + k];
        }
    }
    bpid
}

fn get_hbond_ij(
    i: usize,
    j: usize,
    hb_upper: f64,
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
) -> Vec<(AtomName4, AtomName4, f64)> {
    let ib = seidx[i][1];
    let ie = seidx[i][2];
    let jb = seidx[j][1];
    let je = seidx[j][2];

    let mut out: Vec<(AtomName4, AtomName4, f64)> = Vec::new();
    for m in ib..=ie {
        let am = atom_name[m].0;
        if am[1] == b'P' || am[3] == b'P' {
            continue;
        }
        for n in jb..=je {
            let an = atom_name[n].0;
            if an[1] == b'P' || an[3] == b'P' {
                continue;
            }
            let mut dtmp = [0.0f64; 4];
            for k in 1..=3 {
                dtmp[k] = xyz[m][k] - xyz[n][k];
            }
            let dd = veclen3(&dtmp);
            if dd < hb_upper {
                if out.len() >= BUF512 {
                    return out;
                }
                out.push((atom_name[m], atom_name[n], dd));
            }
        }
    }
    out
}

fn get_unequility(atoms: &[AtomName4]) -> Vec<AtomName4> {
    let mut uniq: Vec<AtomName4> = Vec::new();
    for &a in atoms {
        if uniq.iter().any(|u| u == &a) {
            continue;
        }
        uniq.push(a);
    }
    uniq
}

fn edge_type(atoms: &[AtomName4], base: u8) -> String {
    let mut watson = 0i64;
    let mut hoogsteen = 0i64;
    let mut sugar = 0i64;

    let eq = |a: AtomName4, s: [u8; 4]| atom_eq(a, s);

    let b = base as char;

    for &a in atoms {
        match b {
            'A' | 'a' => {
                if eq(a, *b" N1 ") || eq(a, *b" C2 ") || eq(a, *b" N6 ") {
                    watson += 1;
                }
                if eq(a, *b" N6 ") || eq(a, *b" C5 ") || eq(a, *b" N7 ") || eq(a, *b" C8 ") {
                    hoogsteen += 1;
                }
                if eq(a, *b" C2 ")
                    || eq(a, *b" N3 ")
                    || eq(a, *b" C4 ")
                    || eq(a, *b" N9 ")
                    || eq(a, *b" C1'")
                    || eq(a, *b" C2'")
                    || eq(a, *b" C3'")
                    || eq(a, *b" O3'")
                    || eq(a, *b" O2'")
                {
                    sugar += 1;
                }
            }
            'I' | 'i' => {
                if eq(a, *b" N1 ") || eq(a, *b" C2 ") || eq(a, *b" O6 ") {
                    watson += 1;
                }
                if eq(a, *b" O6 ") || eq(a, *b" C5 ") || eq(a, *b" N7 ") || eq(a, *b" C8 ") {
                    hoogsteen += 1;
                }
                if eq(a, *b" C2 ")
                    || eq(a, *b" N3 ")
                    || eq(a, *b" C4 ")
                    || eq(a, *b" N9 ")
                    || eq(a, *b" C1'")
                    || eq(a, *b" C2'")
                    || eq(a, *b" C3'")
                    || eq(a, *b" O3'")
                    || eq(a, *b" O2'")
                {
                    sugar += 1;
                }
            }
            'G' | 'g' => {
                if eq(a, *b" N2 ") || eq(a, *b" N1 ") || eq(a, *b" O6 ") {
                    watson += 1;
                }
                if eq(a, *b" O6 ") || eq(a, *b" C5 ") || eq(a, *b" N7 ") || eq(a, *b" C8 ") {
                    hoogsteen += 1;
                }
                if eq(a, *b" N2 ")
                    || eq(a, *b" N3 ")
                    || eq(a, *b" C4 ")
                    || eq(a, *b" N9 ")
                    || eq(a, *b" C1'")
                    || eq(a, *b" C2'")
                    || eq(a, *b" O2'")
                {
                    sugar += 1;
                }
            }
            'C' | 'c' => {
                if eq(a, *b" O2 ") || eq(a, *b" N3 ") || eq(a, *b" N4 ") {
                    watson += 1;
                }
                if eq(a, *b" N4 ") || eq(a, *b" C5 ") || eq(a, *b" C6 ") {
                    hoogsteen += 1;
                }
                if eq(a, *b" O2 ")
                    || eq(a, *b" N1 ")
                    || eq(a, *b" C1'")
                    || eq(a, *b" C2'")
                    || eq(a, *b" C3'")
                    || eq(a, *b" O3'")
                    || eq(a, *b" O2'")
                {
                    sugar += 1;
                }
            }
            'U' | 'u' | 'T' | 't' => {
                if eq(a, *b" O2 ") || eq(a, *b" N3 ") || eq(a, *b" O4 ") {
                    watson += 1;
                }
                if eq(a, *b" O4 ") || eq(a, *b" C5 ") || eq(a, *b" C6 ") {
                    hoogsteen += 1;
                }
                if eq(a, *b" O2 ")
                    || eq(a, *b" N1 ")
                    || eq(a, *b" C1'")
                    || eq(a, *b" C2'")
                    || eq(a, *b" C3'")
                    || eq(a, *b" O3'")
                    || eq(a, *b" O2'")
                {
                    sugar += 1;
                }
            }
            'P' | 'p' => {
                if eq(a, *b" O2 ") || eq(a, *b" N3 ") || eq(a, *b" O4 ") {
                    watson += 1;
                }
                if eq(a, *b" O4 ") || eq(a, *b" N1 ") || eq(a, *b" C6 ") {
                    hoogsteen += 1;
                }
                if eq(a, *b" O2 ")
                    || eq(a, *b" C5 ")
                    || eq(a, *b" C1'")
                    || eq(a, *b" C2'")
                    || eq(a, *b" C3'")
                    || eq(a, *b" O3'")
                    || eq(a, *b" O2'")
                {
                    sugar += 1;
                }
            }
            _ => {}
        }
    }

    let max1 = if watson >= hoogsteen { watson } else { hoogsteen };
    let max2 = if hoogsteen >= sugar { hoogsteen } else { sugar };
    let max = if max1 >= max2 { max1 } else { max2 };

    let mut out = if watson == max {
        "W"
    } else if hoogsteen == max {
        "H"
    } else if sugar == max {
        "S"
    } else {
        "?"
    }
    .to_string();

    if max == 0 {
        out = "?".to_string();
    }
    if max == 1 {
        out = ".".to_string();
    }
    out
}

fn get_pair_type(hbonds: &[(AtomName4, AtomName4, f64)], i: usize, j: usize, bseq: &[u8]) -> String {
    if hbonds.is_empty() {
        return "?/?".to_string();
    }

    let a1: Vec<AtomName4> = hbonds.iter().map(|(a, _, _)| *a).collect();
    let a2: Vec<AtomName4> = hbonds.iter().map(|(_, b, _)| *b).collect();

    let uniq1 = get_unequility(&a1);
    let uniq2 = get_unequility(&a2);

    let t1 = edge_type(&uniq1, bseq[i]);
    let t2 = edge_type(&uniq2, bseq[j]);
    format!("{t1}/{t2}")
}

fn test_orientation(
    i: usize,
    j: usize,
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
) -> (i64, i64) {
    let ib = seidx[i][1];
    let ie = seidx[i][2];
    let jb = seidx[j][1];
    let je = seidx[j][2];

    let mut n1 = 0i64;
    let mut n2 = 0i64;

    for k in ib..=ie {
        if !atom_eq(atom_name[k], *b" O2'") {
            continue;
        }
        for k1 in jb..=je {
            let a = atom_name[k1].0;
            if a[1] != b'N' && a[1] != b'O' {
                continue;
            }
            let mut dx = [0.0f64; 4];
            for m in 1..=3 {
                dx[m] = xyz[k][m] - xyz[k1][m];
            }
            if atom_eq(atom_name[k1], *b" O2'") {
                if veclen3(&dx) < 3.8 {
                    n1 += 1;
                }
            }
            if (a[1] == b'N' || a[1] == b'O') && a[3] != b'\'' {
                if veclen3(&dx) < 3.8 {
                    n2 += 1;
                }
            }
        }
        break;
    }
    (n1, n2)
}

fn get_orientation_SS(
    i: usize,
    j: usize,
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
) -> Option<String> {
    let (n1, n2) = test_orientation(i, j, seidx, atom_name, xyz);
    if n1 > 0 && n2 > 0 {
        return Some("s/S".to_string());
    }
    let (n3, n4) = test_orientation(j, i, seidx, atom_name, xyz);
    if n3 > 0 && n4 > 0 {
        return Some("S/s".to_string());
    }
    None
}

fn get_atom_xyz(
    ib: usize,
    ie: usize,
    aname: AtomName4,
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
) -> [f64; 4] {
    let mut out = [0.0f64; 4];
    if let Some(idx) = find_1st_atom(aname, atom_name, ib, ie) {
        for j in 1..=3 {
            out[j] = xyz[idx][j];
        }
    }
    out
}

fn NC_vector(
    i: usize,
    ib: usize,
    ie: usize,
    atom_name: &[AtomName4],
    bseq: &[u8],
    xyz: &[[f64; 4]],
) -> ([f64; 4], [f64; 4], [f64; 4]) {
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

    let c_xyz = get_atom_xyz(ib, ie, AtomName4::from_bytes(*b" C1'"), atom_name, xyz);

    let base = bseq[i] as char;
    let (natm, n_atom) = if matches!(base, 'A' | 'G' | 'a' | 'g' | 'I' | 'i') {
        (9usize, AtomName4::from_bytes(*b" N9 "))
    } else if matches!(base, 'P' | 'p') {
        (6usize, AtomName4::from_bytes(*b" C5 "))
    } else {
        (6usize, AtomName4::from_bytes(*b" N1 "))
    };

    let n_xyz = get_atom_xyz(ib, ie, n_atom, atom_name, xyz);

    let mut vector_nc = [0.0f64; 4];
    for j in 1..=3 {
        vector_nc[j] = c_xyz[j] - n_xyz[j];
    }

    let mut center = [0.0f64; 4];
    let mut n = 0i64;
    for atom in RING_ATOMS.iter().take(natm) {
        if let Some(k) = find_1st_atom(*atom, atom_name, ib, ie) {
            n += 1;
            for j in 1..=3 {
                center[j] += xyz[k][j];
            }
        }
    }
    if n > 0 {
        for j in 1..=3 {
            center[j] /= n as f64;
        }
    }

    (n_xyz, center, vector_nc)
}

fn cis_or_trans(
    i: usize,
    j: usize,
    bseq: &[u8],
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
) -> String {
    let (n_xyz1, xyz1, vector_nc_1) = {
        let ib = seidx[i][1];
        let ie = seidx[i][2];
        NC_vector(i, ib, ie, atom_name, bseq, xyz)
    };
    let (n_xyz2, xyz2, vector_nc_2) = {
        let ib = seidx[j][1];
        let ie = seidx[j][2];
        NC_vector(j, ib, ie, atom_name, bseq, xyz)
    };

    let mut nn_vec = [0.0f64; 4];
    for k in 1..=3 {
        nn_vec[k] = xyz2[k] - xyz1[k];
    }
    let vc1 = cross3(&vector_nc_1, &nn_vec);
    let vc2 = cross3(&vector_nc_2, &nn_vec);
    let a = dot3(&vc1, &vc2);
    if a > 0.0 {
        " cis ".to_string()
    } else {
        " tran".to_string()
    }
}

fn LW_pair_type(
    i: usize,
    j: usize,
    dist: f64,
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
    bseq: &[u8],
) -> String {
    let debug_pair: Option<(usize, usize)> = std::env::var("RNAVIEW_DEBUG_PAIR")
        .ok()
        .and_then(|s| {
            let s = s.trim();
            if s.is_empty() {
                return None;
            }
            let mut it = s.split(|c| c == '_' || c == '-');
            let a = it.next()?.parse::<usize>().ok()?;
            let b = it.next()?.parse::<usize>().ok()?;
            Some((a.min(b), a.max(b)))
        });

    let mut type_field = String::new();
    let cis_tran = cis_or_trans(i, j, bseq, seidx, atom_name, xyz);

    let hbonds = get_hbond_ij(i, j, dist, seidx, atom_name, xyz);
    if debug_pair == Some((i, j)) {
        eprintln!("RNAVIEW_DEBUG_PAIR {i}_{j}: LW_pair_type dist={dist:.3} nh={}", hbonds.len());
        for (a1, a2, d) in &hbonds {
            eprintln!(
                "RNAVIEW_DEBUG_PAIR {i}_{j}: LW hb {}-{} {:.3}",
                as_chars4(a1.0),
                as_chars4(a2.0),
                d
            );
        }
    }
    type_field = get_pair_type(&hbonds, i, j, bseq);

    if type_field.as_bytes().get(0) == Some(&b'S') && type_field.as_bytes().get(2) == Some(&b'S') {
        if let Some(s) = get_orientation_SS(i, j, seidx, atom_name, xyz) {
            type_field = s;
        }
    }

    type_field.push_str(&cis_tran);

    let tmp = pair_code_from_type_field(&type_field);
    if &tmp == b"SWt" || &tmp == b"WSt" {
        let hbonds = get_hbond_ij(i, j, 5.1, seidx, atom_name, xyz);
        type_field = get_pair_type(&hbonds, i, j, bseq);
        type_field.push_str(&cis_tran);
    }
    type_field
}

fn LW_Saenger_correspond(bs1: u8, bs2: u8, type_field: &str) -> String {
    let bs1 = bs1 as char;
    let bs2 = bs2 as char;
    let type_bytes = type_field.as_bytes();
    let t0 = *type_bytes.get(0).unwrap_or(&b'?') as char;
    let t2 = *type_bytes.get(2).unwrap_or(&b'?') as char;
    let t4 = *type_bytes.get(4).unwrap_or(&b'?') as char;

    let mut base = format!("{bs1}{bs2} {bs2}{bs1}").to_ascii_uppercase();
    let mut base_type = format!("{t0}{t2}{t4}  {t2}{t0}{t4}").to_ascii_uppercase();

    let contains = |hay: &str, needle: &str| hay.contains(needle);

    if contains(&base, "GA") && contains(&base_type, "WWC") {
        "VIII".to_string()
    } else if contains(&base, "CC") && contains(&base_type, "WWC") {
        "n/a".to_string()
    } else if (contains(&base, "GU") || contains(&base, "GT")) && contains(&base_type, "WWC") {
        "XXVIII".to_string()
    } else if (contains(&base, "UC") || contains(&base, "TC")) && contains(&base_type, "WWC") {
        "XVIII".to_string()
    } else if (contains(&base, "UU") || contains(&base, "TT")) && contains(&base_type, "WWC") {
        "XVI".to_string()
    } else if (contains(&base, "AU") || contains(&base, "AT")) && contains(&base_type, "--C") {
        "XX".to_string()
    } else if contains(&base, "GC") && contains(&base_type, "++C") {
        "XIX".to_string()
    } else if (contains(&base, "AU") || contains(&base, "AT")) && contains(&base_type, "WWT") {
        "XXI".to_string()
    } else if contains(&base, "AA") && contains(&base_type, "WWT") {
        "I".to_string()
    } else if contains(&base, "GG") && contains(&base_type, "WWT") {
        "III".to_string()
    } else if contains(&base, "GC") && contains(&base_type, "WWT") {
        "XXII".to_string()
    } else if contains(&base, "AC") && contains(&base_type, "WWT") {
        "XXVI".to_string()
    } else if (contains(&base, "GU") || contains(&base, "GT")) && contains(&base_type, "WWT") {
        "XXVII".to_string()
    } else if (contains(&base, "UC") || contains(&base, "TC")) && contains(&base_type, "WWT") {
        "XVII".to_string()
    } else if contains(&base, "CC") && contains(&base_type, "WWT") {
        "XIV,XV".to_string()
    } else if (contains(&base, "UU") || contains(&base, "TT")) && contains(&base_type, "WWT") {
        "XII,XIII".to_string()
    } else if contains(&base, "GG") && contains(&base_type, "WHC") {
        "VI".to_string()
    } else if (contains(&base, "UA") || contains(&base, "TA")) && contains(&base_type, "WHC") {
        "XXIII".to_string()
    } else if contains(&base, "GA") && contains(&base_type, "WHC") {
        "IX".to_string()
    } else if contains(&base, "AA") && contains(&base_type, "WHT") {
        "V".to_string()
    } else if contains(&base, "GG") && contains(&base_type, "WHT") {
        "VII".to_string()
    } else if (contains(&base, "UA") || contains(&base, "TA")) && contains(&base_type, "WHT") {
        "XXIV".to_string()
    } else if contains(&base, "CA") && contains(&base_type, "WHT") {
        "XXV".to_string()
    } else if contains(&base, "AG") && contains(&base_type, "WSC") {
        "n/a".to_string()
    } else if (contains(&base, "AU") || contains(&base, "AT")) && contains(&base_type, "WSC") {
        "n/a".to_string()
    } else if contains(&base, "AG") && contains(&base_type, "WST") {
        "X".to_string()
    } else if contains(&base, "CG") && contains(&base_type, "WST") {
        "n/a".to_string()
    } else if contains(&base_type, "HHC") {
        "n/a".to_string()
    } else if contains(&base, "AA") && contains(&base_type, "HHT") {
        "II".to_string()
    } else if contains(&base, "AG") && contains(&base_type, "HST") {
        "XI".to_string()
    } else if contains(&base, "AA") && contains(&base_type, "HST") {
        "n/a".to_string()
    } else if (contains(&base, "CU") || contains(&base, "CT")) && contains(&base_type, "HST") {
        "n/a".to_string()
    } else if contains(&base, "GG") && contains(&base_type, "SST") {
        "IV".to_string()
    } else {
        "n/a".to_string()
    }
}

pub(crate) fn syn_or_anti(
    num_residue: usize,
    atom_name: &[AtomName4],
    seidx: &[[usize; 3]],
    xyz: &[[f64; 4]],
    ry: &[i32],
) -> Vec<i64> {
    let mut sugar_syn = vec![0i64; num_residue + 1];
    for i in 1..=num_residue {
        let ib = seidx[i][1];
        let ie = seidx[i][2];

        let chi1 = find_1st_atom(AtomName4::from_bytes(*b" O4'"), atom_name, ib, ie);
        let chi2 = find_1st_atom(AtomName4::from_bytes(*b" C1'"), atom_name, ib, ie);

        let (n_atom, c_atom) = if ry.get(i).copied().unwrap_or(0) == 1 {
            (AtomName4::from_bytes(*b" N9 "), AtomName4::from_bytes(*b" C4 "))
        } else {
            (AtomName4::from_bytes(*b" N1 "), AtomName4::from_bytes(*b" C2 "))
        };
        let chi3 = find_1st_atom(n_atom, atom_name, ib, ie);
        let chi4 = find_1st_atom(c_atom, atom_name, ib, ie);

        let Some(a1) = chi1 else { continue };
        let Some(a2) = chi2 else { continue };
        let Some(a3) = chi3 else { continue };
        let Some(a4) = chi4 else { continue };

        let mut d = [[0.0f64; 4]; 5];
        for k in 1..=3 {
            d[1][k] = xyz[a1][k];
            d[2][k] = xyz[a2][k];
            d[3][k] = xyz[a3][k];
            d[4][k] = xyz[a4][k];
        }

        let chi_angle = torsion(&d);
        if (-90.0..=90.0).contains(&chi_angle) {
            sugar_syn[i] = 1;
        } else {
            sugar_syn[i] = 0;
        }
    }
    sugar_syn
}

pub(crate) fn uncommon_lines(
    num_residue: usize,
    seidx: &[[usize; 3]],
    resname: &[ResName3],
    chain_id: &[u8],
    resseq: &[i32],
    icode: &[u8],
    bseq: &[u8],
    ry: &[i32],
) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for i in 1..=num_residue {
        if ry.get(i).copied().unwrap_or(-1) < 0 {
            continue;
        }
        let ib = seidx[i][1];
        let r = resname.get(ib).copied().unwrap_or(ResName3([b' '; 3]));

        let is_standard = resname_eq(r, *b"  A")
            || resname_eq(r, *b"ADE")
            || resname_eq(r, *b"  G")
            || resname_eq(r, *b"GUA")
            || resname_eq(r, *b"  U")
            || resname_eq(r, *b"URA")
            || resname_eq(r, *b"  C")
            || resname_eq(r, *b"CYT")
            || resname_eq(r, *b"  T")
            || resname_eq(r, *b"THY");
        if is_standard {
            continue;
        }

        let resname_str = as_chars3(r.0);
        let chain = *chain_id.get(ib).unwrap_or(&b' ') as char;
        let rseq = *resseq.get(ib).unwrap_or(&0) as i64;
        let ic = *icode.get(ib).unwrap_or(&b' ') as char;
        let b = *bseq.get(i).unwrap_or(&b'?') as char;
        out.push(format!(
            "uncommon residue {resname_str} {rseq:4}{ic} on chain {chain} [#{i}] assigned to: {b}"
        ));
    }
    out
}

pub(crate) struct AllPairsResult {
    pub(crate) base_pairs: Vec<OutBasePairLine>,
    pub(crate) pair_type_codes: Vec<[u8; 3]>,
    pub(crate) pair_info: Vec<Vec<usize>>,
    pub(crate) num_pair_tot: usize,
}

fn candidate_pairs_by_org(
    num_residue: usize,
    org: &[[f64; 4]],
    orien: &[[f64; 10]],
    cutoff: f64,
    dv_max: f64,
) -> Vec<(usize, usize)> {
    if num_residue <= 1 || !(cutoff > 0.0) || !(dv_max > 0.0) {
        return Vec::new();
    }

    // Mirror rnaview_candidate_pairs_by_org (legacy_ffi) to preserve legacy output while
    // avoiding O(n^2) pair enumeration for large structures.
    let cell_size = cutoff;
    let cutoff2 = cutoff * cutoff;
    let eps = 1e-9;

    let n = num_residue;
    let mut coords: Vec<[f64; 3]> = vec![[0.0; 3]; n + 1];
    let mut norm_z: Vec<[f64; 3]> = vec![[0.0; 3]; n + 1];
    let mut cell: Vec<[i64; 3]> = vec![[0; 3]; n + 1];
    let mut min_cell: [i64; 3] = [i64::MAX; 3];
    let mut max_cell: [i64; 3] = [i64::MIN; 3];

    for i in 1..=n {
        let x = org[i][1];
        let y = org[i][2];
        let z = org[i][3];
        coords[i] = [x, y, z];

        let ix = (x / cell_size).floor() as i64;
        let iy = (y / cell_size).floor() as i64;
        let iz = (z / cell_size).floor() as i64;
        cell[i] = [ix, iy, iz];

        norm_z[i] = [orien[i][7], orien[i][8], orien[i][9]];

        if ix < min_cell[0] {
            min_cell[0] = ix;
        }
        if iy < min_cell[1] {
            min_cell[1] = iy;
        }
        if iz < min_cell[2] {
            min_cell[2] = iz;
        }
        if ix > max_cell[0] {
            max_cell[0] = ix;
        }
        if iy > max_cell[1] {
            max_cell[1] = iy;
        }
        if iz > max_cell[2] {
            max_cell[2] = iz;
        }
    }

    let dx = (max_cell[0] - min_cell[0] + 1) as usize;
    let dy = (max_cell[1] - min_cell[1] + 1) as usize;
    let dz = (max_cell[2] - min_cell[2] + 1) as usize;
    let cell_count = match dx.checked_mul(dy).and_then(|v| v.checked_mul(dz)) {
        Some(v) => v,
        None => return Vec::new(),
    };

    let dv = |i: usize, j: usize| -> f64 {
        let zi = norm_z[i];
        let zj = norm_z[j];
        let dd = zi[0] * zj[0] + zi[1] * zj[1] + zi[2] * zj[2];

        let (mut z0, mut z1, mut z2) = if dd <= 0.0 {
            (zi[0] - zj[0], zi[1] - zj[1], zi[2] - zj[2])
        } else {
            (zi[0] + zj[0], zi[1] + zj[1], zi[2] + zj[2])
        };

        let vlen = (z0 * z0 + z1 * z1 + z2 * z2).sqrt();
        if vlen > XEPS {
            z0 /= vlen;
            z1 /= vlen;
            z2 /= vlen;
        }

        let dx0 = coords[j][0] - coords[i][0];
        let dx1 = coords[j][1] - coords[i][1];
        let dx2 = coords[j][2] - coords[i][2];
        (dx0 * z0 + dx1 * z1 + dx2 * z2).abs()
    };

    const MAX_CELLS: usize = 2_000_000;
    if cell_count == 0 || cell_count > MAX_CELLS {
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        let mut js: Vec<usize> = Vec::new();
        let mut extra: Vec<usize> = Vec::new();
        for i in 1..n {
            js.clear();
            extra.clear();
            for j in (i + 1)..=n {
                let di0 = coords[i][0] - coords[j][0];
                let di1 = coords[i][1] - coords[j][1];
                let di2 = coords[i][2] - coords[j][2];
                let d2 = di0 * di0 + di1 * di1 + di2 * di2;
                if d2 <= cutoff2 + eps {
                    js.push(j);
                }
            }

            js.sort_unstable();
            js.dedup();

            for &j in &js {
                let dvj = dv(i, j);
                if dvj > dv_max && dvj <= dv_max + 0.3 {
                    for k in (i + 1..j).rev() {
                        if dv(i, k) <= dv_max {
                            extra.push(k);
                            break;
                        }
                    }
                }
            }
            js.extend(extra.iter().copied());

            for j in (i + 1..=n).rev() {
                if dv(i, j) <= dv_max {
                    js.push(j);
                    break;
                }
            }

            js.sort_unstable();
            js.dedup();
            for &j in &js {
                pairs.push((i, j));
            }
        }
        return pairs;
    }

    let mut head: Vec<i64> = vec![-1; cell_count];
    let mut next: Vec<i64> = vec![-1; n + 1];

    let idx_of = |c: [i64; 3]| -> usize {
        let ox = (c[0] - min_cell[0]) as usize;
        let oy = (c[1] - min_cell[1]) as usize;
        let oz = (c[2] - min_cell[2]) as usize;
        (ox * dy + oy) * dz + oz
    };

    for i in 1..=n {
        let idx = idx_of(cell[i]);
        next[i] = head[idx];
        head[idx] = i as i64;
    }

    let mut pairs: Vec<(usize, usize)> = Vec::new();
    let mut js: Vec<i64> = Vec::new();
    let mut extra: Vec<i64> = Vec::new();

    for i in 1..n {
        js.clear();
        extra.clear();
        let ci = cell[i];
        for ddx in -1..=1 {
            let cx = ci[0] + ddx;
            if cx < min_cell[0] || cx > max_cell[0] {
                continue;
            }
            for ddy in -1..=1 {
                let cy = ci[1] + ddy;
                if cy < min_cell[1] || cy > max_cell[1] {
                    continue;
                }
                for ddz in -1..=1 {
                    let cz = ci[2] + ddz;
                    if cz < min_cell[2] || cz > max_cell[2] {
                        continue;
                    }
                    let idx = idx_of([cx, cy, cz]);
                    let mut j = head[idx];
                    while j != -1 {
                        let ju = j as usize;
                        if ju > i {
                            let dx0 = coords[i][0] - coords[ju][0];
                            let dx1 = coords[i][1] - coords[ju][1];
                            let dx2 = coords[i][2] - coords[ju][2];
                            let d2 = dx0 * dx0 + dx1 * dx1 + dx2 * dx2;
                            if d2 <= cutoff2 + eps {
                                js.push(j);
                            }
                        }
                        j = next[ju];
                    }
                }
            }
        }

        js.sort_unstable();
        js.dedup();

        for &j in &js {
            let ju = j as usize;
            let dvj = dv(i, ju);
            if dvj > dv_max && dvj <= dv_max + 0.3 {
                for k in (i + 1..ju).rev() {
                    if dv(i, k) <= dv_max {
                        extra.push(k as i64);
                        break;
                    }
                }
            }
        }
        js.extend(extra.iter().copied());

        for j in (i + 1..=n).rev() {
            if dv(i, j) <= dv_max {
                js.push(j as i64);
                break;
            }
        }

        js.sort_unstable();
        js.dedup();
        for &j in &js {
            pairs.push((i, j as usize));
        }
    }

    pairs
}

pub(crate) fn all_pairs(
    num_residue: usize,
    ry: &[i32],
    nxyz: &[[f64; 4]],
    orien: &[[f64; 10]],
    org: &[[f64; 4]],
    bprs: &[f64; 7],
    seidx: &[[usize; 3]],
    xyz: &[[f64; 4]],
    atom_name: &[AtomName4],
    chain_id: &[u8],
    resseq: &[i32],
    bseq: &[u8],
) -> Result<AllPairsResult, String> {
    let sugar_syn = syn_or_anti(num_residue, atom_name, seidx, xyz, ry);

    let cand_pairs = candidate_pairs_by_org(num_residue, org, orien, bprs[2], bprs[3]);

    let debug_pair: Option<(usize, usize)> = std::env::var("RNAVIEW_DEBUG_PAIR")
        .ok()
        .and_then(|s| {
            let s = s.trim();
            if s.is_empty() {
                return None;
            }
            let mut it = s.split(|c| c == '_' || c == '-');
            let a = it.next()?.parse::<usize>().ok()?;
            let b = it.next()?.parse::<usize>().ok()?;
            Some((a.min(b), a.max(b)))
        });

    let mut rtn_val = [0.0f64; 21];
    let mut last_angle = 0.0f64;
    let mut have_last_angle = false;
    let mut tmp_str = String::new();

    let mut num_bp: usize = 0;
    let mut pair_type_codes: Vec<[u8; 3]> = Vec::new();
    let mut pair_info: Vec<Vec<usize>> = vec![Vec::new(); num_residue + 1];

    let mut main_lines: Vec<OutBasePairLine> = Vec::new();
    let mut base_single: Vec<OutBasePairLine> = Vec::new();
    let mut base_sugar: Vec<OutBasePairLine> = Vec::new();
    let mut base_p: Vec<OutBasePairLine> = Vec::new();
    let mut other: Vec<OutBasePairLine> = Vec::new();

    for (i, j) in cand_pairs {
        let syn_i = sugar_syn.get(i).copied().unwrap_or(0) == 1;
        let syn_j = sugar_syn.get(j).copied().unwrap_or(0) == 1;

        let bpid = check_pairs(
            i, j, bseq, seidx, xyz, nxyz, orien, org, atom_name, bprs, &mut rtn_val, 0,
        );
        if debug_pair == Some((i, j)) {
            eprintln!(
                "RNAVIEW_DEBUG_PAIR {i}_{j}: bpid={bpid} dorg={:.3} dv={:.3} angle={:.3} dNN={:.3}",
                rtn_val[1], rtn_val[2], rtn_val[3], rtn_val[4]
            );
        }
        if rtn_val[2] <= bprs[3] {
            last_angle = rtn_val[3];
            have_last_angle = true;
        }

        if bpid != 0 {
            let hbonds = Hbond_pair(i, j, seidx, atom_name, bseq, xyz, 0.35, 1, 0);
            let mut nnh = 0i64;
            for (a1, a2, _) in &hbonds {
                if a1.0[3] != b' ' && a2.0[3] != b' ' {
                    continue;
                }
                nnh += 1;
            }
            if debug_pair == Some((i, j)) {
                eprintln!("RNAVIEW_DEBUG_PAIR {i}_{j}: base-base nh={} nnh={nnh}", hbonds.len());
                for (a1, a2, d) in &hbonds {
                    eprintln!(
                        "RNAVIEW_DEBUG_PAIR {i}_{j}: bb hb {}-{} {:.3}",
                        as_chars4(a1.0),
                        as_chars4(a2.0),
                        d
                    );
                }
            }

            if nnh == 1 {
                let type_field = LW_pair_type(i, j, 5.2, seidx, atom_name, xyz, bseq);
                base_single.push(OutBasePairLine {
                    i: i as i32,
                    j: j as i32,
                    chain_i: *chain_id.get(seidx[i][1]).unwrap_or(&b' ') as char,
                    resseq_i: *resseq.get(seidx[i][1]).unwrap_or(&0),
                    base_i: *bseq.get(i).unwrap_or(&b'?') as char,
                    base_j: *bseq.get(j).unwrap_or(&b'?') as char,
                    resseq_j: *resseq.get(seidx[j][1]).unwrap_or(&0),
                    chain_j: *chain_id.get(seidx[j][1]).unwrap_or(&b' ') as char,
                    kind: OutPairKind::Pair,
                    type_field: Some(type_field),
                    syn_i,
                    syn_j,
                    note: Some("!1H(b_b)".to_string()),
                });
                continue;
            }

            let hb_atom1: Vec<AtomName4> = hbonds.iter().map(|(a, _, _)| *a).collect();
            let hb_atom2: Vec<AtomName4> = hbonds.iter().map(|(_, b, _)| *b).collect();
            let nh1 = get_unequility(&hb_atom1).len() as i64;
            let nh2 = get_unequility(&hb_atom2).len() as i64;
            if debug_pair == Some((i, j)) {
                eprintln!("RNAVIEW_DEBUG_PAIR {i}_{j}: nh1={nh1} nh2={nh2}");
            }

            let mut type_field: String = String::new();

            if (nh1 == 1 && nh2 > 1) || (nh2 == 1 && nh1 > 1) {
                if rtn_val[2] > bprs[3] - 0.6 || rtn_val[3] > bprs[4] - 20.0 {
                    let type_field = LW_pair_type(i, j, 5.0, seidx, atom_name, xyz, bseq);
                    base_single.push(OutBasePairLine {
                        i: i as i32,
                        j: j as i32,
                        chain_i: *chain_id.get(seidx[i][1]).unwrap_or(&b' ') as char,
                        resseq_i: *resseq.get(seidx[i][1]).unwrap_or(&0),
                        base_i: *bseq.get(i).unwrap_or(&b'?') as char,
                        base_j: *bseq.get(j).unwrap_or(&b'?') as char,
                        resseq_j: *resseq.get(seidx[j][1]).unwrap_or(&0),
                        chain_j: *chain_id.get(seidx[j][1]).unwrap_or(&b' ') as char,
                        kind: OutPairKind::Pair,
                        type_field: Some(type_field),
                        syn_i,
                        syn_j,
                        note: Some("!1H(b_b).".to_string()),
                    });
                    continue;
                }

                type_field = LW_pair_type(i, j, 4.329, seidx, atom_name, xyz, bseq);
                tmp_str = String::from_utf8_lossy(&pair_code_from_type_field(&type_field)).to_string();

                if tmp_str.contains('.') || tmp_str.contains('?') || rtn_val[2] > bprs[3] - 0.4 || rtn_val[3] > bprs[4] - 15.0 {
                    let type_field = LW_pair_type(i, j, 5.0, seidx, atom_name, xyz, bseq);
                    base_single.push(OutBasePairLine {
                        i: i as i32,
                        j: j as i32,
                        chain_i: *chain_id.get(seidx[i][1]).unwrap_or(&b' ') as char,
                        resseq_i: *resseq.get(seidx[i][1]).unwrap_or(&0),
                        base_i: *bseq.get(i).unwrap_or(&b'?') as char,
                        base_j: *bseq.get(j).unwrap_or(&b'?') as char,
                        resseq_j: *resseq.get(seidx[j][1]).unwrap_or(&0),
                        chain_j: *chain_id.get(seidx[j][1]).unwrap_or(&b' ') as char,
                        kind: OutPairKind::Pair,
                        type_field: Some(type_field),
                        syn_i,
                        syn_j,
                        note: Some("!1H(b_b).".to_string()),
                    });
                    continue;
                }
            } else {
                if bpid == 2 {
                    let bi = to_ascii_upper(bseq[i]);
                    let bj = to_ascii_upper(bseq[j]);
                    if (bi == b'A' && (bj == b'U' || bj == b'T')) || (bj == b'A' && (bi == b'U' || bi == b'T')) {
                        type_field = "-/- cis ".to_string();
                    } else if (bi == b'G' && bj == b'C') || (bj == b'G' && bi == b'C') {
                        type_field = "+/+ cis ".to_string();
                    } else {
                        type_field = "X/X cis ".to_string();
                    }
                } else {
                    type_field = LW_pair_type(i, j, 4.1, seidx, atom_name, xyz, bseq);
                    tmp_str = String::from_utf8_lossy(&pair_code_from_type_field(&type_field)).to_string();
                    if tmp_str.contains('.') || tmp_str.contains('?') {
                        type_field = LW_pair_type(i, j, 4.3, seidx, atom_name, xyz, bseq);
                    }
                }
            }

            let corresp = LW_Saenger_correspond(bseq[i], bseq[j], &type_field);
            main_lines.push(OutBasePairLine {
                i: i as i32,
                j: j as i32,
                chain_i: *chain_id.get(seidx[i][1]).unwrap_or(&b' ') as char,
                resseq_i: *resseq.get(seidx[i][1]).unwrap_or(&0),
                base_i: *bseq.get(i).unwrap_or(&b'?') as char,
                base_j: *bseq.get(j).unwrap_or(&b'?') as char,
                resseq_j: *resseq.get(seidx[j][1]).unwrap_or(&0),
                chain_j: *chain_id.get(seidx[j][1]).unwrap_or(&b' ') as char,
                kind: OutPairKind::Pair,
                type_field: Some(type_field.clone()),
                syn_i,
                syn_j,
                note: Some(corresp),
            });

            num_bp += 1;
            tmp_str = String::from_utf8_lossy(&pair_code_from_type_field(&type_field)).to_string();
            pair_type_codes.push(pair_code_from_type_field(&type_field));

            pair_info[i].push(j);
            pair_info[j].push(i);
            continue;
        }

        // bpid == 0 (non base-base): possible stacking or tertiary interactions.
        if rtn_val[1] > bprs[2] {
            continue;
        }
        if rtn_val[2] > bprs[3] + 0.3 {
            continue;
        }
        let mut angle = rtn_val[3];
        if rtn_val[2] > bprs[3] && have_last_angle {
            angle = last_angle;
        }
        if angle > bprs[4] {
            continue;
        }

        let stack_key = base_stack(i, j, bseq, seidx, atom_name, xyz, &rtn_val);
        if debug_pair == Some((i, j)) {
            eprintln!("RNAVIEW_DEBUG_PAIR {i}_{j}: stack_key={stack_key}");
        }
        if stack_key {
            if rtn_val[3] < 40.0 {
                main_lines.push(OutBasePairLine {
                    i: i as i32,
                    j: j as i32,
                    chain_i: *chain_id.get(seidx[i][1]).unwrap_or(&b' ') as char,
                    resseq_i: *resseq.get(seidx[i][1]).unwrap_or(&0),
                    base_i: *bseq.get(i).unwrap_or(&b'?') as char,
                    base_j: *bseq.get(j).unwrap_or(&b'?') as char,
                    resseq_j: *resseq.get(seidx[j][1]).unwrap_or(&0),
                    chain_j: *chain_id.get(seidx[j][1]).unwrap_or(&b' ') as char,
                    kind: OutPairKind::Stacked,
                    type_field: None,
                    syn_i,
                    syn_j,
                    note: None,
                });
            }
            continue;
        }

        if rtn_val[2] > bprs[3] {
            if debug_pair == Some((i, j)) {
                eprintln!(
                    "RNAVIEW_DEBUG_PAIR {i}_{j}: dv {:.3} > {:.3} (skip tertiary)",
                    rtn_val[2], bprs[3]
                );
            }
            continue;
        }

        let hbonds = Hbond_pair(i, j, seidx, atom_name, bseq, xyz, 0.0, 0, 1);
        if debug_pair == Some((i, j)) {
            eprintln!("RNAVIEW_DEBUG_PAIR {i}_{j}: tertiary nh={}", hbonds.len());
            for (a1, a2, d) in &hbonds {
                eprintln!(
                    "RNAVIEW_DEBUG_PAIR {i}_{j}: hb {}-{} {:.3}",
                    as_chars4(a1.0),
                    as_chars4(a2.0),
                    d
                );
            }
        }
        if hbonds.is_empty() {
            continue;
        }

        let mut nc1 = 0i64;
        let mut nc2 = 0i64;
        for (a1, a2, _) in &hbonds {
            let a1b = a1.0;
            let a2b = a2.0;

            let a1_o2o4 = atom_eq(*a1, *b" O2'") || atom_eq(*a1, *b" O4'");
            let a2_o2o4 = atom_eq(*a2, *b" O2'") || atom_eq(*a2, *b" O4'");

            let a1_is_base_no = is_in(b"NO", a1b[1]) && a1b[3] != b'\'' && a1b[3] != b'P';
            let a2_is_base_no = is_in(b"NO", a2b[1]) && a2b[3] != b'\'' && a2b[3] != b'P';

            if (a1_o2o4 && a2_is_base_no) || (a2_o2o4 && a1_is_base_no) {
                nc1 += 1;
                continue;
            }

            let a1_o1o2p = atom_eq(*a1, *b" O1P") || atom_eq(*a1, *b" O2P");
            let a2_o1o2p = atom_eq(*a2, *b" O1P") || atom_eq(*a2, *b" O2P");

            if (a1_o1o2p && a2b[3] != b'\'' && a2b[1] != b'C') || (a2_o1o2p && a1b[3] != b'\'' && a1b[1] != b'C')
            {
                nc2 += 1;
            }
        }

        let nh = hbonds.len() as i64;
        if nh == nc1 {
            let mut type_field = LW_pair_type(i, j, 4.8, seidx, atom_name, xyz, bseq);
            if tmp_str.contains('.') || tmp_str.contains('?') {
                type_field = LW_pair_type(i, j, 5.8, seidx, atom_name, xyz, bseq);
            }
            base_sugar.push(OutBasePairLine {
                i: i as i32,
                j: j as i32,
                chain_i: *chain_id.get(seidx[i][1]).unwrap_or(&b' ') as char,
                resseq_i: *resseq.get(seidx[i][1]).unwrap_or(&0),
                base_i: *bseq.get(i).unwrap_or(&b'?') as char,
                base_j: *bseq.get(j).unwrap_or(&b'?') as char,
                resseq_j: *resseq.get(seidx[j][1]).unwrap_or(&0),
                chain_j: *chain_id.get(seidx[j][1]).unwrap_or(&b' ') as char,
                kind: OutPairKind::Pair,
                type_field: Some(type_field),
                syn_i,
                syn_j,
                note: Some("!(b_s)".to_string()),
            });
        } else if nh == nc2 {
            let mut type_field = LW_pair_type(i, j, 4.8, seidx, atom_name, xyz, bseq);
            if tmp_str.contains('.') || tmp_str.contains('?') {
                type_field = LW_pair_type(i, j, 5.8, seidx, atom_name, xyz, bseq);
            }
            base_p.push(OutBasePairLine {
                i: i as i32,
                j: j as i32,
                chain_i: *chain_id.get(seidx[i][1]).unwrap_or(&b' ') as char,
                resseq_i: *resseq.get(seidx[i][1]).unwrap_or(&0),
                base_i: *bseq.get(i).unwrap_or(&b'?') as char,
                base_j: *bseq.get(j).unwrap_or(&b'?') as char,
                resseq_j: *resseq.get(seidx[j][1]).unwrap_or(&0),
                chain_j: *chain_id.get(seidx[j][1]).unwrap_or(&b' ') as char,
                kind: OutPairKind::Pair,
                type_field: Some(type_field),
                syn_i,
                syn_j,
                note: Some("!b_(O1P,O2P)".to_string()),
            });
        } else {
            let mut type_field = LW_pair_type(i, j, 5.2, seidx, atom_name, xyz, bseq);
            if tmp_str.contains('.') || tmp_str.contains('?') {
                type_field = LW_pair_type(i, j, 6.0, seidx, atom_name, xyz, bseq);
            }
            other.push(OutBasePairLine {
                i: i as i32,
                j: j as i32,
                chain_i: *chain_id.get(seidx[i][1]).unwrap_or(&b' ') as char,
                resseq_i: *resseq.get(seidx[i][1]).unwrap_or(&0),
                base_i: *bseq.get(i).unwrap_or(&b'?') as char,
                base_j: *bseq.get(j).unwrap_or(&b'?') as char,
                resseq_j: *resseq.get(seidx[j][1]).unwrap_or(&0),
                chain_j: *chain_id.get(seidx[j][1]).unwrap_or(&b' ') as char,
                kind: OutPairKind::Pair,
                type_field: Some(type_field),
                syn_i,
                syn_j,
                note: Some("!(s_s)".to_string()),
            });
        }
    }

    let mut base_pairs = main_lines;
    base_pairs.extend(base_single);
    base_pairs.extend(base_sugar);
    base_pairs.extend(base_p);
    base_pairs.extend(other);

    Ok(AllPairsResult {
        base_pairs,
        pair_type_codes,
        pair_info,
        num_pair_tot: num_bp,
    })
}

fn check_pair_network(
    i: usize,
    j: usize,
    seidx: &[[usize; 3]],
    xyz: &[[f64; 4]],
    nxyz: &[[f64; 4]],
    orien: &[[f64; 10]],
    org: &[[f64; 4]],
    atom_name: &[AtomName4],
    bseq: &[u8],
    bprs: &[f64; 7],
) -> (bool, f64) {
    let mut rtn_val = [0.0f64; 21];
    let bpid = check_pairs(
        i, j, bseq, seidx, xyz, nxyz, orien, org, atom_name, bprs, &mut rtn_val, 1,
    );
    (bpid != 0, rtn_val[2])
}

fn round_mfactor(v: f64) -> i64 {
    (v * MFACTOR).round() as i64
}

fn lsort_indices(a: &mut [i64]) -> Vec<usize> {
    if a.len() <= 2 {
        return (0..a.len()).collect();
    }
    let n = a.len() - 1; // 1-based
    let mut idx: Vec<usize> = (0..=n).collect();

    let mut inc: usize = 1;
    loop {
        inc = inc * 3 + 1;
        if inc > n {
            break;
        }
    }

    while inc > 1 {
        inc /= 3;
        for i in (inc + 1)..=n {
            let v = a[i];
            let iv = idx[i];
            let mut j = i;
            while a[j - inc] > v {
                a[j] = a[j - inc];
                idx[j] = idx[j - inc];
                j -= inc;
                if j <= inc {
                    break;
                }
            }
            a[j] = v;
            idx[j] = iv;
        }
    }

    idx
}

pub(crate) fn bp_network_multiplets(
    num_residue: usize,
    ry: &[i32],
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    chain_id: &[u8],
    resseq: &[i32],
    xyz: &[[f64; 4]],
    bseq: &[u8],
    pair_info_in: &[Vec<usize>],
    nxyz: &[[f64; 4]],
    orien: &[[f64; 10]],
    org: &[[f64; 4]],
    bprs: &[f64; 7],
) -> Result<Option<Vec<String>>, String> {
    let debug_net: Option<usize> = std::env::var("RNAVIEW_DEBUG_NET")
        .ok()
        .and_then(|s| s.trim().parse::<usize>().ok());

    let mut pair_info: Vec<Vec<isize>> = vec![Vec::new(); num_residue + 1];
    for i in 1..=num_residue {
        pair_info[i] = pair_info_in
            .get(i)
            .map(|v| v.iter().map(|&x| x as isize).collect())
            .unwrap_or_default();
    }

    let mut num_ple = 0i64;
    let mut max_ple: i64 = -1;

    for i in 1..=num_residue {
        if ry.get(i).copied().unwrap_or(-1) < 0 {
            continue;
        }

        if debug_net == Some(i) {
            let direct = pair_info_in.get(i).cloned().unwrap_or_default();
            eprintln!("RNAVIEW_DEBUG_NET {i}: direct_pairs={direct:?}");
        }

        let mut ivec: Vec<isize> = vec![0; num_residue + 1];
        ivec[1] = i as isize;
        let mut inum_base: usize = 1;

        let mut m = 1usize;
        while m <= num_residue && ivec[m] != 0 {
            let ir = ivec[m] as usize;
            m += 1;
            for &p in pair_info[ir].iter() {
                let mut seen = false;
                for k in 1..=inum_base {
                    if ivec[k] == p {
                        seen = true;
                        break;
                    }
                }
                if !seen {
                    inum_base += 1;
                    ivec[inum_base] = p;
                }
            }
        }

        if debug_net == Some(i) {
            eprintln!("RNAVIEW_DEBUG_NET {i}: bfs={:?}", &ivec[1..=inum_base]);
        }

        for j in 1..=(inum_base.saturating_sub(1)) {
            if ivec[j] < 0 {
                continue;
            }
            let mut idx1: Vec<i64> = vec![0; inum_base - j + 1];

            let mut mcount = 0usize;
            for k in (j + 1)..=inum_base {
                mcount += 1;
                if ivec[k] < 0 {
                    idx1[mcount] = 1_000_000 - ivec[k] as i64;
                    continue;
                }
                let (ok, dv) = check_pair_network(
                    ivec[j] as usize,
                    ivec[k] as usize,
                    seidx,
                    xyz,
                    nxyz,
                    orien,
                    org,
                    atom_name,
                    bseq,
                    bprs,
                );
                if !ok {
                    idx1[mcount] = 1_000_000 + ivec[k] as i64;
                    ivec[k] = -ivec[k];
                } else {
                    idx1[mcount] = round_mfactor(dv);
                }
            }

            if mcount > 1 {
                let mut tmp: Vec<isize> = vec![0; mcount + 1];
                let order = lsort_indices(&mut idx1);
                for k in 1..=mcount {
                    tmp[k] = ivec[j + order[k]];
                }
                for k in 1..=mcount {
                    ivec[j + k] = tmp[k];
                }
            }
        }

        let mut keep: Vec<isize> = Vec::new();
        for j in 2..=inum_base {
            if ivec[j] > 0 {
                keep.push(ivec[j]);
                if keep.len() >= NP - 1 {
                    break;
                }
            }
        }
        pair_info[i] = keep;

        let k = pair_info[i].len();
        if debug_net == Some(i) {
            let kept: Vec<usize> = pair_info[i].iter().map(|&v| v as usize).collect();
            eprintln!("RNAVIEW_DEBUG_NET {i}: kept_pairs={kept:?} (count={})", kept.len());
        }
        let mut kk = k as i64;
        if kk > 1 {
            num_ple += 1;
            kk += 1;
            if kk > max_ple {
                max_ple = kk;
            }
        }
    }

    if num_ple == 0 {
        return Ok(None);
    }

    let mut unique: Vec<Vec<usize>> = Vec::new();
    let mut lines: Vec<String> = Vec::new();

    for i in 1..=num_residue {
        if pair_info[i].len() <= 1 {
            continue;
        }
        let mut group: Vec<usize> = Vec::with_capacity(pair_info[i].len() + 1);
        group.push(i);
        for &p in &pair_info[i] {
            group.push(p as usize);
        }
        group.sort_unstable();

        if unique.iter().any(|u| u == &group) {
            continue;
        }
        unique.push(group.clone());

        let n_unique = unique.len();
        let inum_base = group.len();
        let mut pairnum = String::new();
        let mut pairstr = String::new();

        for (idx, &jr) in group.iter().enumerate() {
            pairnum.push_str(&format!("{jr}_"));
            let k = seidx[jr][1];
            let chain = *chain_id.get(k).unwrap_or(&b' ') as char;
            let rseq = *resseq.get(k).unwrap_or(&0) as i64;
            let base = *bseq.get(jr).unwrap_or(&b'?') as char;
            if idx + 1 == inum_base {
                pairstr.push_str(&format!("{chain}: {rseq} {base}"));
            } else {
                pairstr.push_str(&format!("{chain}: {rseq} {base}  +  "));
            }
        }

        lines.push(format!("{pairnum}| [{n_unique} {inum_base}]  {pairstr}"));
    }

    Ok(Some(lines))
}

pub(crate) fn pair_type_statistics(pair_type_codes: &[[u8; 3]]) -> ([i64; 7], [i64; 6]) {
    let mut stat = [0i64; 16];
    for t in pair_type_codes {
        if t == b"--c" || t == b"++c" {
            stat[1] += 1;
        } else if t == b"WWc" {
            stat[2] += 1;
        } else if t == b"WWt" {
            stat[3] += 1;
        } else if t == b"HHc" {
            stat[4] += 1;
        } else if t == b"HHt" {
            stat[5] += 1;
        } else if t == b"SSc" || t == b"Ssc" || t == b"sSc" {
            stat[6] += 1;
        } else if t == b"SSt" || t == b"Sst" || t == b"sSt" {
            stat[7] += 1;
        } else if t == b"WHc" || t == b"HWc" {
            stat[8] += 1;
        } else if t == b"WHt" || t == b"HWt" {
            stat[9] += 1;
        } else if t == b"WSc" || t == b"SWc" {
            stat[10] += 1;
        } else if t == b"WSt" || t == b"SWt" {
            stat[11] += 1;
        } else if t == b"HSc" || t == b"SHc" {
            stat[12] += 1;
        } else if t == b"HSt" || t == b"SHt" {
            stat[13] += 1;
        }
    }

    (
        [stat[1], stat[2], stat[3], stat[4], stat[5], stat[6], stat[7]],
        [stat[8], stat[9], stat[10], stat[11], stat[12], stat[13]],
    )
}
