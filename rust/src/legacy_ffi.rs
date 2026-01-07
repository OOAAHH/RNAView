#![allow(non_snake_case)]

use libc::{c_char, c_double, c_int, c_long, c_longlong};

const WC_DORG: c_double = 2.5;
const BUF512: c_long = 512;
const XEPS: c_double = 1.0e-7;

extern "C" {
    fn veclen(va: *mut c_double) -> c_double;
    fn dot(va: *mut c_double, vb: *mut c_double) -> c_double;
    fn vec_norm(va: *mut c_double);
    fn dot2ang(dotval: c_double) -> c_double;

    fn base_stack(
        i: c_long,
        j: c_long,
        bseq: *mut c_char,
        seidx: *mut *mut c_long,
        AtomName: *mut *mut c_char,
        xyz: *mut *mut c_double,
        rtn_val: *mut c_double,
        yes: *mut c_long,
    );

    fn num_strmatch(str: *mut c_char, strmat: *mut *mut c_char, nb: c_long, ne: c_long) -> c_long;

    fn bpstep_par(
        rot1: *mut *mut c_double,
        org1: *mut c_double,
        rot2: *mut *mut c_double,
        org2: *mut c_double,
        pars: *mut c_double,
        mst_orien: *mut *mut c_double,
        mst_org: *mut c_double,
    );

    fn dmatrix(nrl: c_long, nrh: c_long, ncl: c_long, nch: c_long) -> *mut *mut c_double;
    fn free_dmatrix(m: *mut *mut c_double, nrl: c_long, nrh: c_long, ncl: c_long, nch: c_long);

    fn H_catalog(
        i: c_long,
        m: c_long,
        bseq: *mut c_char,
        AtomName: *mut *mut c_char,
        without_H: *mut c_long,
        with_H: *mut c_long,
    );

    fn nrerror(error_text: *mut c_char);

    fn get_hbond_ij(
        i: c_long,
        j: c_long,
        HB_UPPER: *mut c_double,
        seidx: *mut *mut c_long,
        AtomName: *mut *mut c_char,
        HB_ATOM: *mut c_char,
        xyz: *mut *mut c_double,
        nh: *mut c_long,
        hb_atom1: *mut *mut c_char,
        hb_atom2: *mut *mut c_char,
        hb_dist: *mut c_double,
    );

    fn get_pair_type(
        num_hbonds: c_long,
        hb_atom1: *mut *mut c_char,
        hb_atom2: *mut *mut c_char,
        i: c_long,
        j: c_long,
        bseq: *mut c_char,
        type_: *mut c_char,
    );

    fn get_orientation_SS(
        i: c_long,
        j: c_long,
        seidx: *mut *mut c_long,
        AtomName: *mut *mut c_char,
        xyz: *mut *mut c_double,
        type_: *mut c_char,
    );

    fn cis_or_trans(
        i: c_long,
        j: c_long,
        bseq: *mut c_char,
        seidx: *mut *mut c_long,
        AtomName: *mut *mut c_char,
        xyz: *mut *mut c_double,
        cis_tran: *mut c_char,
    );

    fn rnaview_profile_is_enabled() -> c_int;
    fn rnaview_profile_now_ns() -> c_longlong;
    fn rnaview_profile_add_all_pairs_hbond_pair_h_catalog(delta_ns: c_longlong);
    fn rnaview_profile_add_all_pairs_lw_get_hbond_ij(delta_ns: c_longlong);
}

fn to_ascii_upper(c: c_char) -> c_char {
    let b = c as u8;
    if (b'a'..=b'z').contains(&b) {
        (b - 32) as c_char
    } else {
        c
    }
}

fn is_ascii_digit(c: c_char) -> bool {
    (c as u8).is_ascii_digit()
}

fn is_in(chars: &[u8], c: c_char) -> bool {
    chars.contains(&(c as u8))
}

unsafe fn mat_row<'a, T>(m: *mut *mut T, idx: c_long) -> *mut T {
    *m.offset(idx as isize)
}

unsafe fn cstr_eq(a: *const c_char, b: &[u8]) -> bool {
    libc::strcmp(a, b.as_ptr() as *const c_char) == 0
}

#[no_mangle]
pub unsafe extern "C" fn check_pairs(
    i: c_long,
    j: c_long,
    bseq: *mut c_char,
    seidx: *mut *mut c_long,
    xyz: *mut *mut c_double,
    Nxyz: *mut *mut c_double,
    orien: *mut *mut c_double,
    org: *mut *mut c_double,
    AtomName: *mut *mut c_char,
    BPRS: *mut c_double,
    rtn_val: *mut c_double,
    bpid: *mut c_long,
    network: c_long,
) {
    static WC_BYTES: [&[u8]; 9] = [
        b"XX\0",
        b"AT\0",
        b"AU\0",
        b"TA\0",
        b"UA\0",
        b"GC\0",
        b"CG\0",
        b"IC\0",
        b"CI\0",
    ];

    let mut wc_ptrs: [*mut c_char; 9] = [std::ptr::null_mut(); 9];
    for (idx, s) in WC_BYTES.iter().enumerate() {
        wc_ptrs[idx] = s.as_ptr() as *mut c_char;
    }

    *bpid = 0;
    if i == j {
        return;
    }

    let org_i = mat_row(org, i);
    let org_j = mat_row(org, j);
    let nxyz_i = mat_row(Nxyz, i);
    let nxyz_j = mat_row(Nxyz, j);

    let mut dorg: [c_double; 4] = [0.0; 4];
    let mut dnn_vec: [c_double; 4] = [0.0; 4];
    for k in 1..=3 {
        dorg[k] = *org_j.add(k as usize) - *org_i.add(k as usize);
        dnn_vec[k] = *nxyz_j.add(k as usize) - *nxyz_i.add(k as usize);
    }

    *rtn_val.add(1) = veclen(dorg.as_mut_ptr());

    let orien_i = mat_row(orien, i);
    let orien_j = mat_row(orien, j);

    let dir_x = dot(orien_i, orien_j);
    let dir_y = dot(orien_i.add(3), orien_j.add(3));
    let dd = dot(orien_i.add(6), orien_j.add(6));

    let mut zave: [c_double; 4] = [0.0; 4];
    if dd <= 0.0 {
        for k in 1..=3 {
            zave[k] = *orien_i.add((6 + k) as usize) - *orien_j.add((6 + k) as usize);
        }
    } else {
        for k in 1..=3 {
            zave[k] = *orien_i.add((6 + k) as usize) + *orien_j.add((6 + k) as usize);
        }
    }
    vec_norm(zave.as_mut_ptr());

    *rtn_val.add(2) = dot(dorg.as_mut_ptr(), zave.as_mut_ptr()).abs();
    if *rtn_val.add(2) > *BPRS.add(3) {
        return;
    }

    *rtn_val.add(3) = 90.0 - (dot2ang(dd) - 90.0).abs();
    if *rtn_val.add(3) > *BPRS.add(4) {
        return;
    }

    if *rtn_val.add(3) <= 10.0 && *rtn_val.add(2) >= 2.2 {
        return;
    }

    *rtn_val.add(4) = veclen(dnn_vec.as_mut_ptr());
    if *rtn_val.add(4) < *BPRS.add(5) {
        return;
    }

    if j == i + 1 && *rtn_val.add(2) >= 2.0 {
        return;
    }

    *rtn_val.add(5) = *rtn_val.add(1) + 2.0 * *rtn_val.add(2);

    if network != 0 {
        if *rtn_val.add(2) <= *BPRS.add(3) && *rtn_val.add(3) <= *BPRS.add(4) && *rtn_val.add(4) >= *BPRS.add(5)
        {
            *bpid = 1;
        }
        return;
    }

    if *rtn_val.add(1) > *BPRS.add(2) {
        return;
    }

    let seidx_i = mat_row(seidx, i);
    let seidx_j = mat_row(seidx, j);
    let i_start = *seidx_i.add(1) as c_long;
    let i_end = *seidx_i.add(2) as c_long;
    let j_start = *seidx_j.add(1) as c_long;
    let j_end = *seidx_j.add(2) as c_long;

    let mut short_contact: c_long = 0;
    let mut without_h_m: c_long = 0;
    let mut with_h_m: c_long = 0;
    let mut without_h_n: c_long = 0;
    let mut with_h_n: c_long = 0;

    for m in i_start..=i_end {
        if short_contact != 0 {
            break;
        }
        let atom_m = mat_row(AtomName, m);
        let atom_m_1 = *atom_m.add(1);
        let atom_m_0 = *atom_m.add(0);
        let atom_m_2 = *atom_m.add(2);
        let atom_m_3 = *atom_m.add(3);

        for n in j_start..=j_end {
            if short_contact != 0 {
                break;
            }
            let atom_n = mat_row(AtomName, n);
            let atom_n_1 = *atom_n.add(1);
            let atom_n_0 = *atom_n.add(0);
            let atom_n_2 = *atom_n.add(2);
            let atom_n_3 = *atom_n.add(3);

            if !(is_in(b"ON", atom_m_1) && is_in(b"ON", atom_n_1)) {
                continue;
            }
            if !(atom_m_0 == b' ' as c_char && is_ascii_digit(atom_m_2) && atom_m_3 == b' ' as c_char) {
                continue;
            }
            if !(atom_n_0 == b' ' as c_char && is_ascii_digit(atom_n_2) && atom_n_3 == b' ' as c_char) {
                continue;
            }

            H_catalog(i, m, bseq, AtomName, &mut without_h_m, &mut with_h_m);
            H_catalog(j, n, bseq, AtomName, &mut without_h_n, &mut with_h_n);
            if without_h_m == 1 && without_h_n == 1 {
                continue;
            }

            let xyz_m = mat_row(xyz, m);
            let xyz_n = mat_row(xyz, n);
            for k in 1..=3 {
                dorg[k] = *xyz_n.add(k as usize) - *xyz_m.add(k as usize);
            }
            if veclen(dorg.as_mut_ptr()) <= *BPRS.add(1) {
                short_contact = 1;
            }
        }
    }

    if short_contact == 0 {
        return;
    }

    let mut stack_key: c_long = 0;
    base_stack(i, j, bseq, seidx, AtomName, xyz, rtn_val, &mut stack_key);
    if stack_key > 0 {
        return;
    }

    let bi = to_ascii_upper(*bseq.add(i as usize));
    let bj = to_ascii_upper(*bseq.add(j as usize));
    let mut bpi: [c_char; 3] = [0; 3];
    bpi[0] = bi;
    bpi[1] = bj;
    bpi[2] = 0;

    *bpid = -1;
    if dir_x > 0.0 && dir_y < 0.0 && dd < 0.0 {
        *bpid = 1;
        if *rtn_val.add(1) <= WC_DORG
            && num_strmatch(bpi.as_mut_ptr(), wc_ptrs.as_mut_ptr(), 1, 8) != 0
            && *rtn_val.add(3) <= 40.0
            && *rtn_val.add(2) < 1.5
        {
            *bpid = 2;
        }
    }
    if *bpid == 2 {
        *rtn_val.add(5) -= 1.5;
    }

    let r1 = dmatrix(1, 3, 1, 3);
    let r2 = dmatrix(1, 3, 1, 3);
    let mst = dmatrix(1, 3, 1, 3);

    for k in 1..=3 {
        *rtn_val.add((k + 8) as usize) = *orien_i.add((6 + k) as usize);
        *rtn_val.add((k + 11) as usize) = *orien_j.add((6 + k) as usize);

        let koffset = (k - 1) * 3;
        for l in 1..=3 {
            *(*r1.offset(l as isize)).add(k as usize) = *orien_i.add((koffset + l) as usize);
            *(*r2.offset(l as isize)).add(k as usize) = if k == 1 {
                *orien_j.add((koffset + l) as usize)
            } else {
                -*orien_j.add((koffset + l) as usize)
            };
        }
    }

    let mut pars: [c_double; 7] = [0.0; 7];
    bpstep_par(r2, org_j, r1, org_i, pars.as_mut_ptr(), mst, rtn_val.add(5));
    for k in 1..=6 {
        *rtn_val.add((14 + k) as usize) = pars[k as usize];
    }

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
}

#[no_mangle]
pub unsafe extern "C" fn Hbond_pair(
    i: c_long,
    j: c_long,
    seidx: *mut *mut c_long,
    AtomName: *mut *mut c_char,
    bseq: *mut c_char,
    xyz: *mut *mut c_double,
    change: c_double,
    nh: *mut c_long,
    hb_atom1: *mut *mut c_char,
    hb_atom2: *mut *mut c_char,
    hb_dist: *mut c_double,
    c_key: c_long,
    bone_key: c_long,
) {
    const O3P: &[u8] = b" O3'\0";
    const O2P: &[u8] = b" O2P\0";
    const O5P: &[u8] = b" O5'\0";
    const O1P: &[u8] = b" O1P\0";
    const O4P: &[u8] = b" O4'\0";
    const C4: &[u8] = b" C4 \0";
    const C5: &[u8] = b" C5 \0";
    const C6: &[u8] = b" C6 \0";
    const C2: &[u8] = b" C2 \0";

    let seidx_i = mat_row(seidx, i);
    let seidx_j = mat_row(seidx, j);
    let i_start = *seidx_i.add(1) as c_long;
    let i_end = *seidx_i.add(2) as c_long;
    let j_start = *seidx_j.add(1) as c_long;
    let j_end = *seidx_j.add(2) as c_long;

    let mut num_hbonds: c_long = 0;
    let mut without_h_m: c_long = 0;
    let mut with_h_m: c_long = 0;
    let mut without_h_n: c_long = 0;
    let mut with_h_n: c_long = 0;
    let prof = rnaview_profile_is_enabled() != 0;

    for m in i_start..=i_end {
        let atom_m = mat_row(AtomName, m);

        if c_key == 0 && *atom_m.add(1) == b'C' as c_char {
            continue;
        }

        if bone_key == 0
            && (cstr_eq(atom_m, O3P) || cstr_eq(atom_m, O2P) || cstr_eq(atom_m, O5P) || cstr_eq(atom_m, O1P))
        {
            continue;
        }

        if (*atom_m.add(1) == b'C' as c_char && *atom_m.add(3) == b'\'' as c_char) || *atom_m.add(1) == b'P' as c_char
        {
            continue;
        }

        let bi = to_ascii_upper(*bseq.add(i as usize));
        if bi == b'A' as c_char || bi == b'I' as c_char {
            if cstr_eq(atom_m, C4) || cstr_eq(atom_m, C5) || cstr_eq(atom_m, C6) {
                continue;
            }
        } else if bi == b'G' as c_char {
            if cstr_eq(atom_m, C4) || cstr_eq(atom_m, C5) || cstr_eq(atom_m, C6) || cstr_eq(atom_m, C2) {
                continue;
            }
        } else if bi == b'P' as c_char {
            if cstr_eq(atom_m, C4) || cstr_eq(atom_m, C5) {
                continue;
            }
        } else if bi == b'U' as c_char || bi == b'C' as c_char || bi == b'T' as c_char {
            if cstr_eq(atom_m, C4) || cstr_eq(atom_m, C2) {
                continue;
            }
        }

        if prof {
            let t0 = rnaview_profile_now_ns();
            H_catalog(i, m, bseq, AtomName, &mut without_h_m, &mut with_h_m);
            rnaview_profile_add_all_pairs_hbond_pair_h_catalog(rnaview_profile_now_ns() - t0);
        } else {
            H_catalog(i, m, bseq, AtomName, &mut without_h_m, &mut with_h_m);
        }

        for n in j_start..=j_end {
            let atom_n = mat_row(AtomName, n);

            if c_key == 0 && *atom_n.add(1) == b'C' as c_char {
                continue;
            }

            if bone_key == 0
                && (cstr_eq(atom_n, O3P) || cstr_eq(atom_n, O2P) || cstr_eq(atom_n, O5P) || cstr_eq(atom_n, O1P))
            {
                continue;
            }

            if (*atom_n.add(1) == b'C' as c_char && *atom_n.add(3) == b'\'' as c_char) || *atom_n.add(1) == b'P' as c_char
            {
                continue;
            }

            let bj = to_ascii_upper(*bseq.add(j as usize));
            if bj == b'A' as c_char || bj == b'I' as c_char {
                if cstr_eq(atom_n, C4) || cstr_eq(atom_n, C5) || cstr_eq(atom_n, C6) {
                    continue;
                }
            } else if bj == b'G' as c_char {
                if cstr_eq(atom_n, C4) || cstr_eq(atom_n, C5) || cstr_eq(atom_n, C6) || cstr_eq(atom_n, C2) {
                    continue;
                }
            } else if bj == b'P' as c_char {
                if cstr_eq(atom_n, C4) || cstr_eq(atom_n, C5) {
                    continue;
                }
            } else if bj == b'U' as c_char || bj == b'C' as c_char || bj == b'T' as c_char {
                if cstr_eq(atom_n, C4) || cstr_eq(atom_n, C2) {
                    continue;
                }
            }

            if *atom_m.add(1) == b'C' as c_char && *atom_n.add(1) == b'C' as c_char {
                continue;
            }

            if (cstr_eq(atom_m, O3P) || cstr_eq(atom_m, O4P) || cstr_eq(atom_m, O5P) || cstr_eq(atom_m, O1P) || cstr_eq(atom_m, O2P))
                && (cstr_eq(atom_n, O3P) || cstr_eq(atom_n, O4P) || cstr_eq(atom_n, O5P) || cstr_eq(atom_n, O1P) || cstr_eq(atom_n, O2P))
            {
                continue;
            }

            if prof {
                let t0 = rnaview_profile_now_ns();
                H_catalog(j, n, bseq, AtomName, &mut without_h_n, &mut with_h_n);
                rnaview_profile_add_all_pairs_hbond_pair_h_catalog(rnaview_profile_now_ns() - t0);
            } else {
                H_catalog(j, n, bseq, AtomName, &mut without_h_n, &mut with_h_n);
            }
            if without_h_m == 1 && without_h_n == 1 {
                continue;
            }

            let dist = if is_in(b"NO", *atom_m.add(1)) && *atom_m.add(3) != b'\'' as c_char && *atom_m.add(3) != b'P' as c_char
                && is_in(b"NO", *atom_n.add(1)) && *atom_n.add(3) != b'\'' as c_char && *atom_n.add(3) != b'P' as c_char
            {
                let mut d = 3.4 + change;
                if d >= 4.0 {
                    d = 4.0;
                }
                d
            } else if (*atom_m.add(1) == b'C' as c_char
                && *atom_n.add(3) != b'\'' as c_char
                && *atom_n.add(3) != b'P' as c_char
                && is_in(b"NO", *atom_n.add(1)))
                || (*atom_n.add(1) == b'C' as c_char
                    && *atom_m.add(3) != b'\'' as c_char
                    && *atom_m.add(3) != b'P' as c_char
                    && is_in(b"NO", *atom_m.add(1)))
            {
                let mut d = 3.6 + change;
                if d >= 4.0 {
                    d = 4.0;
                }
                d
            } else if (*atom_m.add(1) == b'O' as c_char
                && *atom_m.add(3) == b'\'' as c_char
                && is_in(b"NO", *atom_n.add(1))
                && *atom_n.add(3) != b'\'' as c_char
                && *atom_n.add(3) != b'P' as c_char)
                || (*atom_n.add(1) == b'O' as c_char
                    && *atom_n.add(3) == b'\'' as c_char
                    && is_in(b"NO", *atom_m.add(1))
                    && *atom_m.add(3) != b'\'' as c_char
                    && *atom_m.add(3) != b'P' as c_char)
            {
                let mut d = 3.4 + change;
                if d >= 4.0 {
                    d = 4.0;
                }
                d
            } else if (*atom_m.add(3) == b'P' as c_char && *atom_n.add(3) != b'\'' as c_char && *atom_n.add(1) != b'C' as c_char)
                || (*atom_n.add(3) == b'P' as c_char && *atom_m.add(3) != b'\'' as c_char && *atom_m.add(1) != b'C' as c_char)
            {
                let mut d = 3.2 + change;
                if d >= 4.0 {
                    d = 4.0;
                }
                d
            } else {
                let mut d = 3.1 + change;
                if d >= 3.8 {
                    d = 3.8;
                }
                d
            };

            let xyz_m = mat_row(xyz, m);
            let xyz_n = mat_row(xyz, n);
            let mut dtmp: [c_double; 4] = [0.0; 4];
            for k in 1..=3 {
                dtmp[k] = *xyz_m.add(k as usize) - *xyz_n.add(k as usize);
            }
            let dd = veclen(dtmp.as_mut_ptr());
            if dd < dist {
                num_hbonds += 1;
                if num_hbonds > BUF512 {
                    let msg = b"Too many possible H-bonds between two bases\0";
                    nrerror(msg.as_ptr() as *mut c_char);
                    *nh = BUF512;
                    return;
                }
                libc::strcpy(mat_row(hb_atom1, num_hbonds), atom_m);
                libc::strcpy(mat_row(hb_atom2, num_hbonds), atom_n);
                *hb_dist.add(num_hbonds as usize) = dd;
            }
        }
    }

    *nh = num_hbonds;
}

#[no_mangle]
pub unsafe extern "C" fn LW_pair_type(
    i: c_long,
    j: c_long,
    dist: c_double,
    seidx: *mut *mut c_long,
    AtomName: *mut *mut c_char,
    HB_ATOM: *mut c_char,
    xyz: *mut *mut c_double,
    bseq: *mut c_char,
    hb_atom1: *mut *mut c_char,
    hb_atom2: *mut *mut c_char,
    hb_dist: *mut c_double,
    type_: *mut c_char,
) {
    *type_ = 0;
    let prof = rnaview_profile_is_enabled() != 0;
    let mut cis_tran: [c_char; 10] = [0; 10];

    let mut hb_upper_new: [c_double; 3] = [0.0; 3];
    hb_upper_new[0] = dist;

    let mut nh: c_long = 0;
    if prof {
        let t0 = rnaview_profile_now_ns();
        get_hbond_ij(
            i,
            j,
            hb_upper_new.as_mut_ptr(),
            seidx,
            AtomName,
            HB_ATOM,
            xyz,
            &mut nh,
            hb_atom1,
            hb_atom2,
            hb_dist,
        );
        rnaview_profile_add_all_pairs_lw_get_hbond_ij(rnaview_profile_now_ns() - t0);
    } else {
        get_hbond_ij(
            i,
            j,
            hb_upper_new.as_mut_ptr(),
            seidx,
            AtomName,
            HB_ATOM,
            xyz,
            &mut nh,
            hb_atom1,
            hb_atom2,
            hb_dist,
        );
    }

    get_pair_type(nh, hb_atom1, hb_atom2, i, j, bseq, type_);

    if *type_.add(0) == b'S' as c_char && *type_.add(2) == b'S' as c_char {
        get_orientation_SS(i, j, seidx, AtomName, xyz, type_);
    }

    cis_or_trans(i, j, bseq, seidx, AtomName, xyz, cis_tran.as_mut_ptr());
    libc::strcat(type_, cis_tran.as_ptr());

    let mut tmp: [u8; 4] = [0; 4];
    tmp[0] = *type_.add(0) as u8;
    tmp[1] = *type_.add(2) as u8;
    tmp[2] = *type_.add(4) as u8;
    tmp[3] = 0;

    if tmp == *b"SWt\0" || tmp == *b"WSt\0" {
        *type_ = 0;
        hb_upper_new[0] = 5.1;
        if prof {
            let t0 = rnaview_profile_now_ns();
            get_hbond_ij(
                i,
                j,
                hb_upper_new.as_mut_ptr(),
                seidx,
                AtomName,
                HB_ATOM,
                xyz,
                &mut nh,
                hb_atom1,
                hb_atom2,
                hb_dist,
            );
            rnaview_profile_add_all_pairs_lw_get_hbond_ij(rnaview_profile_now_ns() - t0);
        } else {
            get_hbond_ij(
                i,
                j,
                hb_upper_new.as_mut_ptr(),
                seidx,
                AtomName,
                HB_ATOM,
                xyz,
                &mut nh,
                hb_atom1,
                hb_atom2,
                hb_dist,
            );
        }
        get_pair_type(nh, hb_atom1, hb_atom2, i, j, bseq, type_);
        libc::strcat(type_, cis_tran.as_ptr());
    }
}

#[no_mangle]
pub unsafe extern "C" fn rnaview_candidate_pairs_by_org(
    num_residue: c_long,
    org: *mut *mut c_double,
    orien: *mut *mut c_double,
    cutoff: c_double,
    dv_max: c_double,
    out_len_pairs: *mut c_long,
) -> *mut c_long {
    if out_len_pairs.is_null() {
        return std::ptr::null_mut();
    }
    *out_len_pairs = 0;

    if num_residue <= 1 || org.is_null() || orien.is_null() || !(cutoff > 0.0) || !(dv_max > 0.0) {
        return std::ptr::null_mut();
    }

    let n = num_residue as usize;
    let cell_size = cutoff;
    let cutoff2 = cutoff * cutoff;
    let eps = 1e-9;

    let mut coords: Vec<[c_double; 3]> = vec![[0.0; 3]; n + 1];
    let mut norm_z: Vec<[c_double; 3]> = vec![[0.0; 3]; n + 1];
    let mut cell: Vec<[i64; 3]> = vec![[0; 3]; n + 1];
    let mut min_cell: [i64; 3] = [i64::MAX; 3];
    let mut max_cell: [i64; 3] = [i64::MIN; 3];

    for i in 1..=n {
        let row = mat_row(org, i as c_long);
        let x = *row.add(1);
        let y = *row.add(2);
        let z = *row.add(3);
        coords[i] = [x, y, z];

        let ix = (x / cell_size).floor() as i64;
        let iy = (y / cell_size).floor() as i64;
        let iz = (z / cell_size).floor() as i64;
        cell[i] = [ix, iy, iz];

        let o = mat_row(orien, i as c_long);
        norm_z[i] = [*o.add(7), *o.add(8), *o.add(9)];

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

    let Some(cell_count) = dx.checked_mul(dy).and_then(|v| v.checked_mul(dz)) else {
        return std::ptr::null_mut();
    };

    let dv = |i: usize, j: usize| -> c_double {
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
        let mut pairs: Vec<c_long> = Vec::new();
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
                pairs.push(i as c_long);
                pairs.push(j as c_long);
            }
        }
        let pair_count = (pairs.len() / 2) as c_long;
        if pair_count == 0 {
            return std::ptr::null_mut();
        }
        let bytes = pairs.len() * std::mem::size_of::<c_long>();
        let ptr = libc::malloc(bytes) as *mut c_long;
        if ptr.is_null() {
            return std::ptr::null_mut();
        }
        std::ptr::copy_nonoverlapping(pairs.as_ptr(), ptr, pairs.len());
        *out_len_pairs = pair_count;
        return ptr;
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

    let mut pairs: Vec<c_long> = Vec::new();
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
            pairs.push(i as c_long);
            pairs.push(j as c_long);
        }
    }

    let pair_count = (pairs.len() / 2) as c_long;
    if pair_count == 0 {
        return std::ptr::null_mut();
    }

    let bytes = pairs.len() * std::mem::size_of::<c_long>();
    let ptr = libc::malloc(bytes) as *mut c_long;
    if ptr.is_null() {
        return std::ptr::null_mut();
    }
    std::ptr::copy_nonoverlapping(pairs.as_ptr(), ptr, pairs.len());
    *out_len_pairs = pair_count;
    ptr
}
