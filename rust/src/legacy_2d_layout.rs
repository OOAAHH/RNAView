// Legacy 2D layout (best-pair selection + helix partitioning + XY coordinates).
//
// Ported from legacy C:
// - Best-pair selection + re-ordering: src/rnaview.c, src/fpair.c, src/fpair_sub.c
// - 2D layout: src/ps-xy.c, src/ps-xy-sub.c
//
// This module intentionally uses 1-based indexing (like legacy) to maximize output parity.

use crate::legacy_alg::{compute_base_info, cross, dot, dot2ang, find_1st_atom, vec_norm, veclen, AtomName4};
use crate::legacy_pairing::check_pairs;
use crate::PairsJson;

const BUF512: usize = 512;
const EMPTY_NUMBER: f64 = -9999.99;
const XBIG: f64 = 1.0e18;
const HLXANG: f64 = 0.26;
const MFACTOR: f64 = 10000.0;
const XEPS: f64 = 1.0e-7;

#[derive(Debug, Clone)]
pub(crate) struct Layout2d {
    pub(crate) base_xy: Vec<[f64; 4]>, // 1-based: [_, x, y, _]
    pub(crate) xml_nh: usize,
    pub(crate) xml_helix: Vec<[i64; 3]>, // 1-based rows: [_, left, right] in working residue indices
    pub(crate) xml_helix_len: Vec<i64>,  // 1-based: length per helix
    pub(crate) xml_ns: usize,
    pub(crate) xml_bases: Vec<i64>, // 1-based: residue indices
    pub(crate) num_loop: usize,
    pub(crate) loop_: Vec<[i64; 3]>, // 1-based rows: [_, start, end]
}

pub(crate) fn compute_layout_2d(
    num_residue: usize,
    bseq: &[u8],
    sugar_class: &[u8],
    seidx: &[[usize; 3]],
    ry: &[i32],
    atom_name: &[AtomName4],
    chain_id: &[u8],
    xyz: &mut [[f64; 4]],
    pairs_all: &PairsJson,
) -> Result<Layout2d, String> {
    if num_residue == 0 {
        return Err("compute_layout_2d: num_residue=0".to_string());
    }

    let base_info = compute_base_info(num_residue, bseq, seidx, ry, atom_name, xyz)?;

    // Candidate list derived from pairs.json core pairs; this corresponds to legacy `pair_info`.
    let mut cand: Vec<Vec<usize>> = vec![Vec::new(); num_residue + 1];
    for bp in &pairs_all.core.base_pairs {
        if bp.kind != "pair" {
            continue;
        }
        let i = bp.i;
        let j = bp.j;
        if i <= 0 || j <= 0 {
            continue;
        }
        let iu = i as usize;
        let ju = j as usize;
        if iu == 0 || ju == 0 || iu > num_residue || ju > num_residue {
            continue;
        }
        cand[iu].push(ju);
        cand[ju].push(iu);
    }
    for v in &mut cand {
        v.sort_unstable();
        v.dedup();
    }

    // Best-pair selection (mutual best match).
    let mut matched_idx: Vec<i64> = vec![0; num_residue + 1];
    let mut base_pairs: Vec<[i64; 18]> = vec![[0; 18]]; // 1-based

    for i in 1..=num_residue {
        let Some(pi) = best_pair(
            i,
            &cand,
            ry,
            &matched_idx,
            bseq,
            sugar_class,
            seidx,
            xyz,
            &base_info.nxyz,
            &base_info.orien,
            &base_info.org,
            atom_name,
            &base_info.bprs,
        )? else {
            continue;
        };
        let j = pi.partner;
        let Some(pj) = best_pair(
            j,
            &cand,
            ry,
            &matched_idx,
            bseq,
            sugar_class,
            seidx,
            xyz,
            &base_info.nxyz,
            &base_info.orien,
            &base_info.org,
            atom_name,
            &base_info.bprs,
        )? else {
            continue;
        };
        if pj.partner != i {
            continue;
        }

        matched_idx[i] = 1;
        matched_idx[j] = 1;

        let mut row = pi.as_base_pairs_row();
        row[1] = i as i64;
        base_pairs.push(row);
    }

    let num_bp = base_pairs.len().saturating_sub(1);
    if num_bp == 0 {
        return Err("compute_layout_2d: no best pairs".to_string());
    }

    // Re-order base-pairs into helix regions (legacy `re_ordering`).
    let o3_p = build_o3_p(num_residue, seidx, atom_name, xyz, ry)?;
    let mut bp_idx: Vec<i64> = vec![0; num_bp + 1];
    let mut helix_marker: Vec<i64> = vec![0; num_bp + 1];
    let mut helix_idx: Vec<[i64; 8]> = vec![[0; 8]; num_bp + 1]; // [][1..7]
    let mut num_helix: i64 = 1;
    re_ordering(
        num_bp as i64,
        &mut base_pairs,
        &mut bp_idx,
        &mut helix_marker,
        &mut helix_idx,
        &base_info.bprs,
        &mut num_helix,
        &o3_p,
    )?;
    if num_helix <= 0 {
        return Err("compute_layout_2d: no helix".to_string());
    }

    // 2D layout (legacy `process_2d_fig`).
    let mut base_xy: Vec<[f64; 4]> = vec![[0.0; 4]; num_residue + 1];
    let mut loop_: Vec<[i64; 3]> = vec![[0; 3]; (num_helix as usize) * 2 + 2]; // 1-based rows: [_, start, end]
    let mut xml_helix: Vec<[i64; 3]> = vec![[0; 3]; num_residue + 1];
    let mut xml_helix_len: Vec<i64> = vec![0; num_residue + 1];
    let mut xml_bases: Vec<i64> = vec![0; num_residue + 1];
    let mut num_loop: i64 = 0;
    let mut xml_nh: i64 = 0;
    let mut xml_ns: i64 = 0;

    process_2d_fig(
        num_residue as i64,
        seidx,
        ry,
        atom_name,
        chain_id,
        xyz,
        num_helix,
        &helix_idx,
        &bp_idx,
        &base_pairs,
        &mut base_xy,
        &mut num_loop,
        &mut loop_,
        &mut xml_nh,
        &mut xml_helix_len,
        &mut xml_helix,
        &mut xml_ns,
        &mut xml_bases,
        bseq,
    )?;

    if xml_nh <= 0 {
        return Err("compute_layout_2d: xml_nh=0 (no anti-parallel helix)".to_string());
    }

    Ok(Layout2d {
        base_xy,
        xml_nh: xml_nh as usize,
        xml_helix,
        xml_helix_len,
        xml_ns: xml_ns as usize,
        xml_bases,
        num_loop: num_loop as usize,
        loop_,
    })
}

#[derive(Debug, Clone)]
struct PairStat {
    partner: usize,
    bpid: i64,
    rounded: [i64; 14], // get_round(MFACTOR*rtn_val[1..14])
}

impl PairStat {
    fn as_base_pairs_row(&self) -> [i64; 18] {
        let mut row = [0i64; 18];
        row[2] = self.partner as i64;
        row[3] = self.bpid;
        for k in 0..14 {
            row[4 + k] = self.rounded[k];
        }
        row
    }
}

fn get_round(d: f64) -> i64 {
    if d > 0.0 {
        (d + 0.5) as i64
    } else {
        (d - 0.5) as i64
    }
}

fn best_pair(
    i: usize,
    cand: &[Vec<usize>],
    ry: &[i32],
    matched_idx: &[i64],
    bseq: &[u8],
    sugar_class: &[u8],
    seidx: &[[usize; 3]],
    xyz: &[[f64; 4]],
    nxyz: &[[f64; 4]],
    orien: &[[f64; 10]],
    org: &[[f64; 4]],
    atom_name: &[AtomName4],
    bprs: &[f64; 7],
) -> Result<Option<PairStat>, String> {
    let mut ddmin = XBIG;
    let mut best: Option<PairStat> = None;

    let mut rtn_val = [0.0f64; 21];
    let Some(cands) = cand.get(i) else {
        return Ok(None);
    };
    for &j in cands {
        if j == i {
            continue;
        }
        if ry.get(j).copied().unwrap_or(-1) < 0 {
            continue;
        }
        if matched_idx.get(j).copied().unwrap_or(0) != 0 {
            continue;
        }

        let bpid = check_pairs(
            i,
            j,
            bseq,
            sugar_class,
            seidx,
            xyz,
            nxyz,
            orien,
            org,
            atom_name,
            bprs,
            &mut rtn_val,
            0,
            None,
        );
        if bpid == 0 {
            continue;
        }
        if rtn_val[5] >= ddmin {
            continue;
        }
        ddmin = rtn_val[5];

        let mut rounded = [0i64; 14];
        for k in 1..=14 {
            rounded[k - 1] = get_round(MFACTOR * rtn_val[k]);
        }
        best = Some(PairStat {
            partner: j,
            bpid,
            rounded,
        });
    }
    Ok(best)
}

fn build_o3_p(
    num_residue: usize,
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
    ry: &[i32],
) -> Result<Vec<[f64; 9]>, String> {
    let mut out: Vec<[f64; 9]> = vec![[0.0; 9]; num_residue + 1];
    for i in 1..=num_residue {
        if ry.get(i).copied().unwrap_or(-1) < 0 {
            out[i][4] = -1.0;
            out[i][8] = -1.0;
            continue;
        }
        let ib = seidx.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0);
        let ie = seidx.get(i).and_then(|v| v.get(2).copied()).unwrap_or(0);
        if ib == 0 || ie == 0 || ib > ie {
            out[i][4] = -1.0;
            out[i][8] = -1.0;
            continue;
        }
        o3_p_xyz(
            ib,
            ie,
            AtomName4::from_bytes(*b" O3'"),
            atom_name,
            xyz,
            &mut out[i],
            4,
        );
        o3_p_xyz(
            ib,
            ie,
            AtomName4::from_bytes(*b" P  "),
            atom_name,
            xyz,
            &mut out[i],
            8,
        );
    }
    Ok(out)
}

fn o3_p_xyz(
    ib: usize,
    ie: usize,
    aname: AtomName4,
    atom_name: &[AtomName4],
    xyz: &[[f64; 4]],
    o3_or_p: &mut [f64; 9],
    idx: usize,
) {
    let found = find_1st_atom(aname, atom_name, ib, ie);
    if let Some(i) = found {
        for j in 1..=3 {
            o3_or_p[idx - 4 + j] = xyz[i][j];
        }
        o3_or_p[idx] = 1.0;
    } else {
        o3_or_p[idx] = -1.0;
    }
}

fn re_ordering(
    num_bp: i64,
    base_pairs: &mut [[i64; 18]],
    bp_idx: &mut [i64],
    helix_marker: &mut [i64],
    helix_idx: &mut [[i64; 8]],
    bprs: &[f64; 7],
    num_helix: &mut i64,
    o3_p: &[[f64; 9]],
) -> Result<(), String> {
    if num_bp <= 0 {
        *num_helix = 0;
        return Ok(());
    }

    let num_bp_usize = num_bp as usize;
    let mut bp_xyz: Vec<[f64; 10]> = vec![[0.0; 10]; num_bp_usize + 1]; // 1-based [][1..9]
    for i in 1..=num_bp_usize {
        for j in 1..=9 {
            bp_xyz[i][j] = (base_pairs[i][j + 8] as f64) / MFACTOR;
        }
    }

    let mut bp_order: Vec<[i64; 4]> = vec![[0; 4]; num_bp_usize + 1];
    let mut end_list: Vec<[i64; 4]> = vec![[0; 4]; num_bp_usize + 1];
    let mut num_ends: i64 = 0;

    bp_context1(
        num_bp,
        base_pairs,
        bprs.get(6).copied().unwrap_or(0.0),
        &bp_xyz,
        &mut bp_order,
        &mut end_list,
        &mut num_ends,
    );

    locate_helix1(
        num_bp,
        helix_idx,
        num_ends,
        num_helix,
        &end_list,
        &bp_order,
        bp_idx,
        helix_marker,
    );

    five2three(
        num_bp,
        num_helix,
        helix_idx,
        bp_idx,
        &mut bp_xyz,
        base_pairs,
        o3_p,
    );

    check_zdna(num_helix, helix_idx, bp_idx, &bp_xyz);

    Ok(())
}

fn bp_context1(
    num_bp: i64,
    _base_pairs: &[[i64; 18]],
    helix_chg: f64,
    bp_xyz: &[[f64; 10]],
    bp_order: &mut [[i64; 4]],
    end_list: &mut [[i64; 4]],
    num_ends: &mut i64,
) {
    let mut temp = 0.0f64;
    let mut overlap = 0i64;
    let cnum = 8usize;

    let num_bp_usize = num_bp as usize;
    for i in 1..=num_bp_usize {
        let mut ddmin = [0.0f64; 9];
        let mut ddidx = [0i64; 9];
        for j in 1..=cnum {
            ddmin[j] = XBIG;
            ddidx[j] = 0;
        }

        for j in 1..=num_bp_usize {
            if j == i {
                continue;
            }
            let mut txyz = [0.0f64; 4];
            for k in 1..=3 {
                txyz[k] = bp_xyz[i][k] - bp_xyz[j][k];
            }
            temp = veclen(&txyz);
            for k in 1..=cnum {
                if temp < ddmin[k] {
                    for m in (k + 1..=cnum).rev() {
                        let n = m - 1;
                        if ddidx[n] != 0 {
                            ddmin[m] = ddmin[n];
                            ddidx[m] = ddidx[n];
                        }
                    }
                    ddmin[k] = temp;
                    ddidx[k] = j as i64;
                    break;
                }
            }
        }

        if ddidx[1] != 0 && ddidx[2] != 0 {
            if ddmin[1] > helix_chg {
                *num_ends += 1;
                end_list[*num_ends as usize][1] = i as i64;
            } else {
                if ddmin[1] < 1.25 {
                    overlap += 1;
                }
                let mut txyz = [0.0f64; 4];
                for j in 1..=3 {
                    txyz[j] = bp_xyz[i][j] - bp_xyz[ddidx[1] as usize][j];
                }
                vec_norm(&mut txyz);
                let mut n: i64 = 0;

                for j in 2..=cnum {
                    if ddidx[j] == 0 {
                        continue;
                    }
                    let mut txyz2 = [0.0f64; 4];
                    for k in 1..=3 {
                        txyz2[k] = bp_xyz[i][k] - bp_xyz[ddidx[j] as usize][k];
                    }
                    vec_norm(&mut txyz2);
                    if dot(&txyz, &txyz2) < HLXANG {
                        if ddmin[j] <= helix_chg {
                            n = j as i64;
                            bp_order[i][1] = -1;
                            bp_order[i][2] = ddidx[1];
                            bp_order[i][3] = ddidx[j];
                        } else {
                            n = 9999;
                        }
                        break;
                    }
                }

                if n == 0 || n == 9999 {
                    let n = 2;
                    *num_ends += 1;
                    end_list[*num_ends as usize][1] = i as i64;
                    end_list[*num_ends as usize][2] = ddidx[1];
                    bp_order[i][2] = ddidx[1];
                    let mut txyz2 = [0.0f64; 4];
                    for j in 1..=3 {
                        txyz2[j] = bp_xyz[ddidx[2] as usize][j] - bp_xyz[ddidx[1] as usize][j];
                    }
                    if dot(&txyz, &txyz2) < 0.0 && veclen(&txyz2) <= helix_chg {
                        end_list[*num_ends as usize][3] = ddidx[n];
                        bp_order[i][3] = ddidx[n];
                    }
                }
            }
        }
    }

    if *num_ends == 0 {
        *num_ends = 1;
        end_list[1][1] = 1;
        if num_bp == 2 {
            if temp <= helix_chg {
                end_list[1][2] = 2;
                *num_ends = 2;
                end_list[2][1] = 2;
                end_list[2][2] = 1;
            } else {
                *num_ends = 2;
                end_list[2][1] = 2;
            }
        }
    }

    if overlap != 0 {
        // legacy prints a warning; ignore for deterministic output.
    }
}

fn locate_helix1(
    num_bp: i64,
    helix_idx: &mut [[i64; 8]],
    num_ends: i64,
    num_helix: &mut i64,
    end_list: &[[i64; 4]],
    bp_order: &[[i64; 4]],
    bp_idx: &mut [i64],
    helix_marker: &mut [i64],
) {
    let num_bp_usize = num_bp as usize;
    let mut ip: i64 = 0;

    helix_idx[*num_helix as usize][1] = 1;
    let mut matched_idx: Vec<i64> = vec![0; num_bp_usize + 1];

    for i in 1..=num_ends {
        if ip >= num_bp {
            break;
        }
        let mut k_sum = 0i64;
        let mut k0 = 0i64;
        for j in 1..=3 {
            let k = end_list[i as usize][j];
            if k != 0 {
                k_sum += matched_idx[k as usize];
                k0 += 1;
            }
        }
        if k_sum == k0 {
            continue;
        }
        for j in 1..=3 {
            if ip >= num_bp {
                break;
            }
            let k = end_list[i as usize][j];
            if k != 0 && matched_idx[k as usize] == 0 {
                ip += 1;
                bp_idx[ip as usize] = k;
                matched_idx[k as usize] = 1;
            }
        }

        for _j in 1..=num_bp {
            let k = bp_idx[ip as usize] as usize;
            let k2 = bp_order[k][2];
            let k3 = bp_order[k][3];
            if bp_order[k][1] == 0 {
                if k2 != 0 && matched_idx[k2 as usize] == 0 && k3 == 0 {
                    ip += 1;
                    bp_idx[ip as usize] = k2;
                    matched_idx[k2 as usize] = 1;
                }
                break;
            }
            let m = matched_idx[k2 as usize] + matched_idx[k3 as usize];
            if m == 2 || m == 0 {
                break;
            }
            let prev = bp_idx[(ip - 1) as usize];
            if k2 == prev {
                ip += 1;
                bp_idx[ip as usize] = k3;
                matched_idx[k3 as usize] = 1;
            } else if k3 == prev {
                ip += 1;
                bp_idx[ip as usize] = k2;
                matched_idx[k2 as usize] = 1;
            } else {
                break;
            }
        }

        helix_idx[*num_helix as usize][2] = ip;
        helix_marker[ip as usize] = 1;
        if ip < num_bp {
            *num_helix += 1;
            helix_idx[*num_helix as usize][1] = ip + 1;
        }
    }

    if ip < num_bp {
        helix_idx[*num_helix as usize][7] = 1;
        helix_idx[*num_helix as usize][2] = num_bp;
        helix_marker[num_bp as usize] = 1;
        for j in 1..=num_bp_usize {
            if matched_idx[j] == 0 {
                ip += 1;
                bp_idx[ip as usize] = j as i64;
            }
        }
    }
}

fn distance_ab(o3_p: &[[f64; 9]], ia: i64, ib: i64, ipa: usize, ipb: usize) -> f64 {
    let iau = ia as usize;
    let ibu = ib as usize;
    if iau >= o3_p.len() || ibu >= o3_p.len() {
        return -1.0;
    }
    if o3_p[iau][ipa] > 0.0 && o3_p[ibu][ipb] > 0.0 {
        let mut txyz = [0.0f64; 4];
        for i in 1..=3 {
            txyz[i] = o3_p[iau][ipa - 4 + i] - o3_p[ibu][ipb - 4 + i];
        }
        veclen(&txyz)
    } else {
        -1.0
    }
}

fn get_ij(m: i64, swapped: &[i64], base_pairs: &[[i64; 18]]) -> (i64, i64) {
    let mu = m as usize;
    if swapped.get(mu).copied().unwrap_or(0) != 0 {
        (base_pairs[mu][2], base_pairs[mu][1])
    } else {
        (base_pairs[mu][1], base_pairs[mu][2])
    }
}

fn get_d1_d2(m: i64, n: i64, swapped: &[i64], bp_xyz: &[[f64; 10]]) -> (f64, f64) {
    let mut dm = [0.0f64; 4];
    let mut dn = [0.0f64; 4];
    let mut zave = [0.0f64; 4];
    let mut dorg = [0.0f64; 4];

    let mu = m as usize;
    let nu = n as usize;

    for i in 1..=3 {
        dorg[i] = bp_xyz[nu][i] - bp_xyz[mu][i];
        let (idx1, idx2) = if swapped.get(mu).copied().unwrap_or(0) != 0 { (6, 3) } else { (3, 6) };
        dm[i] = bp_xyz[mu][i + idx1] - bp_xyz[mu][i + idx2];
        let (idx1, idx2) = if swapped.get(nu).copied().unwrap_or(0) != 0 { (6, 3) } else { (3, 6) };
        dn[i] = bp_xyz[nu][i + idx1] - bp_xyz[nu][i + idx2];
    }
    vec_norm(&mut dm);
    vec_norm(&mut dn);

    let mut d1 = dot(&dm, &dn);
    let d2: f64;
    if d1 < 0.0 && d1 > -HLXANG {
        d1 = 1.0;
        d2 = 3.4;
    } else {
        for i in 1..=3 {
            zave[i] = dm[i] + if d1 > 0.0 { dn[i] } else { -dn[i] };
        }
        vec_norm(&mut zave);
        d2 = dot(&zave, &dorg);
    }
    (d1, d2)
}

fn deg2rad(ang_deg: f64) -> f64 {
    ang_deg * std::f64::consts::PI / 180.0
}

fn magang(mut a: [f64; 4], mut b: [f64; 4]) -> f64 {
    if veclen(&a) < XEPS || veclen(&b) < XEPS {
        return 0.0;
    }
    vec_norm(&mut a);
    vec_norm(&mut b);
    dot2ang(dot(&a, &b))
}

fn identity_matrix_3() -> [[f64; 4]; 4] {
    let mut out = [[0.0f64; 4]; 4];
    for i in 1..=3 {
        out[i][i] = 1.0;
    }
    out
}

fn arb_rotation(mut va: [f64; 4], ang_deg: f64) -> [[f64; 4]; 4] {
    let vlen = veclen(&va);
    if vlen < XEPS {
        return identity_matrix_3();
    }
    for i in 1..=3 {
        va[i] /= vlen;
    }
    let ang = deg2rad(ang_deg);
    let c = ang.cos();
    let s = ang.sin();
    let dc = 1.0 - c;

    let mut rot = [[0.0f64; 4]; 4];
    rot[1][1] = c + dc * va[1] * va[1];
    rot[1][2] = va[1] * va[2] * dc - va[3] * s;
    rot[1][3] = va[1] * va[3] * dc + va[2] * s;
    rot[2][1] = va[1] * va[2] * dc + va[3] * s;
    rot[2][2] = c + dc * va[2] * va[2];
    rot[2][3] = va[2] * va[3] * dc - va[1] * s;
    rot[3][1] = va[1] * va[3] * dc - va[2] * s;
    rot[3][2] = va[2] * va[3] * dc + va[1] * s;
    rot[3][3] = c + dc * va[3] * va[3];
    rot
}

fn rot_2_lsplane(num_atoms: usize, atom_name: &[AtomName4], xyz: &mut [[f64; 4]]) -> Result<(), String> {
    if num_atoms == 0 {
        return Ok(());
    }

    let atom_p = AtomName4::from_bytes(*b" P  ");
    let atom_o5 = AtomName4::from_bytes(*b" O5'");
    let atom_c5 = AtomName4::from_bytes(*b" C5'");
    let atom_c4 = AtomName4::from_bytes(*b" C4'");
    let atom_c3 = AtomName4::from_bytes(*b" C3'");
    let atom_o3 = AtomName4::from_bytes(*b" O3'");

    let mut plane_points: Vec<[f64; 4]> = vec![[0.0; 4]];
    for j in 1..=num_atoms {
        let an = atom_name.get(j).copied().unwrap_or(AtomName4([0; 4]));
        if an != atom_p && an != atom_o5 && an != atom_c5 && an != atom_c4 && an != atom_c3 && an != atom_o3 {
            continue;
        }
        let p = xyz.get(j).copied().unwrap_or([0.0; 4]);
        plane_points.push(p);
    }
    let n = plane_points.len().saturating_sub(1);
    if n < 3 {
        return Ok(());
    }

    let mut adist: Vec<f64> = vec![0.0; n + 2];
    let mut z = [0.0f64; 4];
    let mut ppos = [0.0f64; 4];
    let mut odist = 0.0f64;
    ls_plane(
        &plane_points,
        n as i64,
        &mut z,
        &mut ppos,
        &mut odist,
        &mut adist,
    )?;

    let zphy = [0.0f64, 0.0, 0.0, 1.0];
    let hinge = cross(&z, &zphy);
    let ang_deg = magang(zphy, z);
    let rotmat = arb_rotation(hinge, ang_deg);

    for j in 1..=num_atoms {
        if j >= xyz.len() {
            break;
        }
        let old = xyz[j];
        let mut rotated = [0.0f64; 4];
        for k in 1..=3 {
            rotated[k] = dot(&old, &rotmat[k]);
        }
        xyz[j][1] = rotated[1];
        xyz[j][2] = rotated[2];
        xyz[j][3] = rotated[3];
    }
    Ok(())
}

fn lreverse(ia: i64, n: i64, lvec: &mut [i64]) {
    let n_usize = n as usize;
    let ia_usize = ia as usize;
    let mut tmp = vec![0i64; n_usize + 1];
    for i in 1..=n_usize {
        tmp[i] = lvec[n_usize + ia_usize - i];
    }
    for i in 1..=n_usize {
        lvec[ia_usize + i - 1] = tmp[i];
    }
}

fn five2three(
    num_bp: i64,
    num_helix: &mut i64,
    helix_idx: &mut [[i64; 8]],
    bp_idx: &mut [i64],
    bp_xyz: &mut [[f64; 10]],
    base_pairs: &mut [[i64; 18]],
    o3_p: &[[f64; 9]],
) {
    let o3p = 5.0f64;
    let num_bp_usize = num_bp as usize;

    // Check wrong O3'[i] to P[i] linkage.
    for i in 1..=*num_helix as usize {
        for j in helix_idx[i][1]..=helix_idx[i][2] {
            let i1 = base_pairs[bp_idx[j as usize] as usize][1];
            let di1_i2 = distance_ab(o3_p, i1, i1, 4, 8);
            if di1_i2 > 0.0 && di1_i2 <= 2.5 {
                return;
            }
        }
    }

    let mut swapped: Vec<i64> = vec![0; num_bp_usize + 1];
    for i in 1..=*num_helix as usize {
        helix_idx[i][3] = helix_idx[i][2] - helix_idx[i][1] + 1;
        let mut nwc: i64 = 0;
        for j in helix_idx[i][1]..helix_idx[i][2] {
            let m = bp_idx[j as usize];
            let n = bp_idx[(j + 1) as usize];
            if base_pairs[m as usize][3] > 0 && base_pairs[n as usize][3] > 0 {
                nwc += 1;
                let (d1, d2) = get_d1_d2(m, n, &swapped, bp_xyz);
                if d1 < 0.0 {
                    swapped[n as usize] = 1 - swapped[n as usize];
                    if d2 < 0.0 {
                        swapped[m as usize] = 1 - swapped[m as usize];
                    }
                } else {
                    if d2 < 0.0 {
                        swapped[n as usize] = 1 - swapped[n as usize];
                        if nwc == 1 {
                            swapped[m as usize] = 1 - swapped[m as usize];
                        }
                    }
                }
                let (d1, d2) = get_d1_d2(m, n, &swapped, bp_xyz);
                if d1 < 0.0 {
                    swapped[n as usize] = 1 - swapped[n as usize];
                }
                if d2 < 0.0 {
                    // legacy prints weird warning; ignore.
                }
            } else {
                nwc = 0;
                let (i1, j1) = get_ij(m, &swapped, base_pairs);
                let (i2, j2) = get_ij(n, &swapped, base_pairs);
                let di1_i2 = distance_ab(o3_p, i1, i2, 4, 4);
                let di1_j2 = distance_ab(o3_p, i1, j2, 4, 4);
                let dj1_i2 = distance_ab(o3_p, j1, i2, 4, 4);
                let dj1_j2 = distance_ab(o3_p, j1, j2, 4, 4);
                if (di1_i2 > 0.0 && di1_j2 > 0.0 && di1_i2 > di1_j2)
                    && (dj1_i2 > 0.0 && dj1_j2 > 0.0 && dj1_j2 > dj1_i2)
                {
                    swapped[n as usize] = 1 - swapped[n as usize];
                }
            }
        }

        // Check if strand I is in 5'->3' direction.
        let mut direction = [0i64; 9];
        for j in helix_idx[i][1]..helix_idx[i][2] {
            let m = bp_idx[j as usize];
            let (i1, j1) = get_ij(m, &swapped, base_pairs);
            let n = bp_idx[(j + 1) as usize];
            let (i2, j2) = get_ij(n, &swapped, base_pairs);

            let di1_i2 = distance_ab(o3_p, i1, i2, 4, 8);
            let dj2_j1 = distance_ab(o3_p, j2, j1, 4, 8);
            if di1_i2 > o3p {
                direction[1] += 1;
            } else if di1_i2 > 0.0 {
                direction[2] += 1;
            }
            if dj2_j1 > o3p {
                direction[3] += 1;
            } else if dj2_j1 > 0.0 {
                direction[4] += 1;
            }

            let di2_i1 = distance_ab(o3_p, i2, i1, 4, 8);
            let dj1_j2 = distance_ab(o3_p, j1, j2, 4, 8);
            if di2_i1 > o3p {
                direction[5] += 1;
            } else if di2_i1 > 0.0 {
                direction[6] += 1;
            }
            if dj1_j2 > o3p {
                direction[7] += 1;
            } else if dj1_j2 > 0.0 {
                direction[8] += 1;
            }
        }

        if (direction[1] - direction[2]) * (direction[5] - direction[6]) > 0
            || (direction[3] - direction[4]) * (direction[7] - direction[8]) > 0
        {
            helix_idx[i][7] = 1;
        } else {
            if (direction[1] != 0 && direction[2] != 0)
                || (direction[3] != 0 && direction[4] != 0)
                || (direction[5] != 0 && direction[6] != 0)
                || (direction[7] != 0 && direction[8] != 0)
            {
                helix_idx[i][5] = 1;
            }
            if direction[1] > direction[2] || direction[5] < direction[6] {
                if direction[3] > direction[4] || direction[7] < direction[8] {
                    let m = bp_idx[helix_idx[i][1] as usize];
                    let (i1, j1) = get_ij(m, &swapped, base_pairs);
                    let n = bp_idx[helix_idx[i][2] as usize];
                    let (i2, _j2) = get_ij(n, &swapped, base_pairs);
                    if i2 < j1 {
                        lreverse(helix_idx[i][1], helix_idx[i][3], bp_idx);
                    } else {
                        for j in helix_idx[i][1]..=helix_idx[i][2] {
                            let idx = bp_idx[j as usize] as usize;
                            swapped[idx] = 1 - swapped[idx];
                        }
                    }
                } else {
                    helix_idx[i][6] = 1;
                    lreverse(helix_idx[i][1], helix_idx[i][3], bp_idx);
                }
            } else {
                if direction[3] > direction[4] || direction[7] < direction[8] {
                    helix_idx[i][6] = 1;
                }
            }
        }

        let mut nswap = 0i64;
        for j in helix_idx[i][1]..=helix_idx[i][2] {
            if swapped[bp_idx[j as usize] as usize] != 0 {
                nswap += 1;
            }
        }
        if nswap != 0 {
            for j in helix_idx[i][1]..=helix_idx[i][2] {
                let m = bp_idx[j as usize] as usize;
                if swapped[m] != 0 {
                    // swap the two residue indices
                    let t = base_pairs[m][1];
                    base_pairs[m][1] = base_pairs[m][2];
                    base_pairs[m][2] = t;
                    // swap base normals
                    for k in 1..=3 {
                        let t = bp_xyz[m][k + 3];
                        bp_xyz[m][k + 3] = bp_xyz[m][k + 6];
                        bp_xyz[m][k + 6] = t;

                        let t = base_pairs[m][k + 11];
                        base_pairs[m][k + 11] = base_pairs[m][k + 14];
                        base_pairs[m][k + 14] = t;
                    }
                }
            }
        }
    }
}

fn check_zdna(num_helix: &mut i64, helix_idx: &mut [[i64; 8]], bp_idx: &[i64], bp_xyz: &[[f64; 10]]) {
    let mut nwired = 0i64;
    let mut mixed_rl = 0i64;

    for i in 1..=*num_helix as usize {
        if helix_idx[i][5] != 0 || helix_idx[i][6] != 0 || helix_idx[i][7] != 0 || helix_idx[i][3] <= 1 {
            nwired += 1;
            continue;
        }
        let mut nrev = 0i64;
        let mut txyz = [0.0f64; 4];
        for j in helix_idx[i][1]..=helix_idx[i][2] {
            let m = bp_idx[j as usize] as usize;
            if helix_idx[i][3] == 1 {
                continue;
            }
            if j < helix_idx[i][2] {
                let n = bp_idx[(j + 1) as usize] as usize;
                for k in 1..=3 {
                    txyz[k] = bp_xyz[n][k] - bp_xyz[m][k];
                }
            }
            if dot(&txyz, &[0.0, bp_xyz[m][4], bp_xyz[m][5], bp_xyz[m][6]]) < 0.0 {
                nrev += 1;
            } else {
                break;
            }
        }
        if nrev == helix_idx[i][3] {
            helix_idx[i][4] = 1;
            mixed_rl += 1;
        }
    }

    if nwired == 0 && mixed_rl != 0 && mixed_rl != *num_helix {
        // legacy prints a warning; ignore.
    }
}

// === ps-xy.c + ps-xy-sub.c port ===

fn process_2d_fig(
    num_residue: i64,
    seidx: &[[usize; 3]],
    ry: &[i32],
    atom_name: &[AtomName4],
    chain_id: &[u8],
    xyz: &mut [[f64; 4]],
    num_helixs: i64,
    helix_idx: &[[i64; 8]],
    bp_idx: &[i64],
    base_pairs: &[[i64; 18]],
    xy_bs: &mut [[f64; 4]],
    num_loop1: &mut i64,
    loop_: &mut [[i64; 3]],
    xmlnh: &mut i64,
    xml_helix_len: &mut [i64],
    xml_helix: &mut [[i64; 3]],
    xml_ns: &mut i64,
    xml_bases: &mut [i64],
    bseq: &[u8],
) -> Result<(), String> {
    let num_residue_usize = num_residue as usize;
    if num_residue <= 0 {
        *xmlnh = 0;
        return Ok(());
    }

    let mut bs_sub: Vec<i64> = vec![0; num_residue_usize + 2];
    let mut helix_len: Vec<f64> = vec![0.0; (num_helixs as usize) + 1];
    let mut dxyz: Vec<[f64; 4]> = vec![[0.0; 4]; (num_helixs as usize) + 1];
    let mut slope: Vec<f64> = vec![0.0; (num_helixs as usize) + 1];

    let mut xy: Vec<[f64; 4]> = vec![[0.0; 4]; (num_helixs as usize) * 2 + 2];
    let mut nxy: Vec<[f64; 4]> = vec![[0.0; 4]; (num_helixs as usize) * 2 + 2];

    let mut bs_1: Vec<Vec<i64>> = vec![vec![0; num_residue_usize + 2]; (num_helixs as usize) + 1];
    let mut bs_2: Vec<Vec<i64>> = vec![vec![0; num_residue_usize + 2]; (num_helixs as usize) + 1];
    let mut bs1_lk: Vec<Vec<i64>> = vec![vec![0; num_residue_usize + 2]; (num_helixs as usize) + 1];
    let mut bs2_lk: Vec<Vec<i64>> = vec![vec![0; num_residue_usize + 2]; (num_helixs as usize) + 1];
    let mut link: Vec<i64> = vec![0; num_residue_usize + 2];
    let mut xy_bso: Vec<[f64; 4]> = vec![[0.0; 4]; num_residue_usize + 1];
    let mut bsp_1: Vec<Vec<i64>> = vec![vec![0; num_residue_usize + 2]; (num_helixs as usize) + 1];
    let mut bsp_2: Vec<Vec<i64>> = vec![vec![0; num_residue_usize + 2]; (num_helixs as usize) + 1];
    let mut bs_1_add: Vec<i64> = vec![0; num_residue_usize + 2];
    let mut bs_2_add: Vec<i64> = vec![0; num_residue_usize + 2];
    let mut add: Vec<i64> = vec![0; num_residue_usize + 2];
    let mut npair_per_helix: Vec<i64> = vec![0; (num_helixs as usize) + 1];
    let mut num_per_helix: Vec<i64> = vec![0; (num_helixs as usize) + 1];
    let mut nregion: Vec<i64> = vec![0; num_residue_usize + 2];
    let mut sub_helix1: Vec<Vec<i64>> = vec![vec![0; num_residue_usize + 2]; (num_helixs as usize) + 1];
    let mut sub_helix2: Vec<Vec<i64>> = vec![vec![0; num_residue_usize + 2]; (num_helixs as usize) + 1];

    // get rid of the isolated pair
    let mut nhelix: i64 = 0;
    helix_regions(
        num_helixs,
        helix_idx,
        bp_idx,
        base_pairs,
        &mut nhelix,
        &mut npair_per_helix,
        &mut bs_1,
        &mut bs_2,
    );

    // only count the continuing anti-parallel helix regions in the longer helix.
    let mut k_out = 0i64;
    for i in 1..=nhelix as usize {
        head_to_tail(
            i as i64,
            &mut npair_per_helix,
            &mut bs_1,
            &mut bs_2,
            &mut nregion,
            &mut sub_helix1,
            &mut sub_helix2,
        );
        let n = npair_per_helix[i];
        if n == 0 || nregion[i] == 0 {
            continue;
        }
        k_out += 1;
        for j in 1..=n as usize {
            bs_1[k_out as usize][j] = bs_1[i][j];
            bs_2[k_out as usize][j] = bs_2[i][j];
        }
        npair_per_helix[k_out as usize] = npair_per_helix[i];
        nregion[k_out as usize] = nregion[i];
        for j in 1..=nregion[k_out as usize] as usize {
            sub_helix1[k_out as usize][j] = sub_helix1[i][j];
            sub_helix2[k_out as usize][j] = sub_helix2[i][j];
        }
    }

    nhelix = k_out;
    if nhelix <= 0 {
        *xmlnh = 0;
        return Ok(());
    }

    for i in 1..=nhelix as usize {
        // make sure the first increases and second decreases within each region
        for j in 1..=nregion[i] as usize {
            let m = sub_helix1[i][j] as usize; // first element in region
            let n = sub_helix2[i][j] as usize; // last element in region
            if bs_1[i][n] < bs_1[i][m] {
                // descending order: swap strands for all pairs in region
                for k in sub_helix1[i][j]..=sub_helix2[i][j] {
                    let ku = k as usize;
                    let t = bs_1[i][ku];
                    bs_1[i][ku] = bs_2[i][ku];
                    bs_2[i][ku] = t;
                }
            }
        }
    }

    // NOTE: for rnaml where 5' and 3' are defined
    let mut xml_nh = 0i64;
    for k in 1..=nhelix as usize {
        for j in 1..=nregion[k] as usize {
            xml_nh += 1;
            xml_helix_len[xml_nh as usize] = sub_helix2[k][j] - sub_helix1[k][j] + 1;
            let tmp = sub_helix1[k][j] as usize;
            xml_helix[xml_nh as usize][1] = bs_1[k][tmp];
            xml_helix[xml_nh as usize][2] = bs_2[k][tmp];
        }
    }
    *xmlnh = xml_nh;

    for k in 1..=xml_nh as usize {
        if xml_helix[k][1] > xml_helix[k][2] {
            let i = xml_helix[k][1] + xml_helix_len[k] - 1;
            let j = xml_helix[k][2] - xml_helix_len[k] + 1;
            if j < 0 || i < 0 {
                return Err("process_2d_fig: error in defining helix".to_string());
            }
            xml_helix[k][1] = j;
            xml_helix[k][2] = i;
        }
    }

    // get the rest of bases by subtracting the helix pairs from the num_residue
    let mut nsub: i64 = 0;
    rest_bases(
        num_residue,
        nhelix,
        &npair_per_helix,
        &bs_1,
        &bs_2,
        &mut nsub,
        &mut bs_sub,
    );

    *xml_ns = nsub;
    for i in 1..=nsub as usize {
        xml_bases[i] = bs_sub[i];
    }

    // loop check at both ends
    let mut loop_key: Vec<[i64; 3]> = vec![[0; 3]; (num_helixs as usize) + 1];
    let mut num_loop: i64 = 0;
    for i in 1..=nhelix as usize {
        loop_key[i][1] = 0;
        loop_key[i][2] = 0;

        let n = npair_per_helix[i];
        let mut loop_tmp = [0i64; 3];
        let mut yes = 0i64;

        helix_head(
            i as i64,
            1,
            &bs_1,
            &bs_2,
            &mut nsub,
            &mut bs_sub,
            chain_id,
            seidx,
            &mut loop_tmp,
            &mut yes,
        );
        if yes > 0 {
            loop_key[i][1] = 1;
            num_loop += 1;
            loop_[num_loop as usize][1] = loop_tmp[1];
            loop_[num_loop as usize][2] = loop_tmp[2];
        }

        helix_tail(
            i as i64,
            n,
            &bs_1,
            &bs_2,
            &mut nsub,
            &mut bs_sub,
            chain_id,
            seidx,
            &mut loop_tmp,
            &mut yes,
        );
        if yes > 0 {
            loop_key[i][2] = 1;
            num_loop += 1;
            loop_[num_loop as usize][1] = loop_tmp[1];
            loop_[num_loop as usize][2] = loop_tmp[2];
        }
    }
    *num_loop1 = num_loop;

    // Complete each longer helix and make the chain residues continue.
    for i in 1..=nhelix as usize {
        let n = npair_per_helix[i];
        if n == 0 || nregion[i] <= 0 {
            continue;
        }

        let mut m = 0i64;
        link[i] = 0;
        let mut nn = 0i64;

        if nregion[i] > 1 {
            for j in 1..=(nregion[i] - 1) as usize {
                let n1 = sub_helix2[i][j];
                let k = j + 1;
                let n2 = sub_helix1[i][k];
                if (k as i64) > nregion[i] {
                    break;
                }

                add[j] = 0;
                let mut yes = 0i64;
                check_link(i as i64, n1, n2, nsub, &bs_sub, &bs_1, &bs_2, &mut yes);

                if yes <= 0 {
                    // the two segment of real helix is not linked
                    add[j] = 1;
                    bs_1_add[add[j] as usize] = 0;
                    bs_2_add[add[j] as usize] = 0;

                    nn += 1;
                    bs1_lk[i][nn as usize] = bs_1[i][n1 as usize];
                    bs2_lk[i][nn as usize] = bs_2[i][n1 as usize];
                    nn += 1;
                    bs1_lk[i][nn as usize] = bs_1[i][n2 as usize];
                    bs2_lk[i][nn as usize] = bs_2[i][n2 as usize];
                    link[i] = nn;
                }

                if yes > 0 {
                    add_bs_2helix(
                        i as i64,
                        j as i64,
                        n1,
                        n2,
                        num_residue,
                        &bs_1,
                        &bs_2,
                        &mut nsub,
                        &mut bs_sub,
                        &mut add,
                        &mut bs_1_add,
                        &mut bs_2_add,
                    );
                }

                for k2 in sub_helix1[i][j]..=sub_helix2[i][j] {
                    m += 1;
                    let ku = k2 as usize;
                    bsp_1[i][m as usize] = bs_1[i][ku];
                    bsp_2[i][m as usize] = bs_2[i][ku];
                }
                for k2 in 1..=add[j] {
                    m += 1;
                    let ku = k2 as usize;
                    bsp_1[i][m as usize] = bs_1_add[ku];
                    bsp_2[i][m as usize] = bs_2_add[ku];
                }
            }
        }

        let j = nregion[i] as usize;
        for k2 in sub_helix1[i][j]..=sub_helix2[i][j] {
            m += 1;
            let ku = k2 as usize;
            bsp_1[i][m as usize] = bs_1[i][ku];
            bsp_2[i][m as usize] = bs_2[i][ku];
        }
        num_per_helix[i] = m;
    }

    // determine the longest helix axis length (pre-rotation)
    let mut longest: i64 = 1;
    let mut maxdm = -XBIG;
    for i in 1..=nhelix as usize {
        let n = npair_per_helix[i];
        let mut xy1 = [0.0f64; 4];
        let mut xy2 = [0.0f64; 4];
        let mut a = 0.0f64;
        HelixAxis(
            i as i64,
            n,
            &bs_1,
            &bs_2,
            num_residue,
            seidx,
            ry,
            xyz,
            atom_name,
            bseq,
            &mut a,
            &mut xy1,
            &mut xy2,
            &mut helix_len,
            &mut dxyz,
        )?;
        if helix_len[i] >= maxdm {
            maxdm = helix_len[i];
            longest = i as i64;
        }
    }

    // rotate the molecule to the least-squares plane (legacy behavior)
    let mut atmnum: usize = 0;
    for i in 1..=num_residue_usize {
        let ib = seidx.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0);
        let ie = seidx.get(i).and_then(|v| v.get(2).copied()).unwrap_or(0);
        if ie >= ib {
            atmnum += ie - ib + 1;
        }
    }
    rot_2_lsplane(atmnum, atom_name, xyz)?;

    // determine the slope and front/tail of each helix (post-rotation)
    let mut mm: i64 = 0;
    let mut helix_2d_vec = [0.0f64; 4];
    for i in 1..=nhelix as usize {
        let n = npair_per_helix[i];
        let mut xy1 = [0.0f64; 4];
        let mut xy2 = [0.0f64; 4];
        let mut a = 0.0f64;
        HelixAxis(
            i as i64,
            n,
            &bs_1,
            &bs_2,
            num_residue,
            seidx,
            ry,
            xyz,
            atom_name,
            bseq,
            &mut a,
            &mut xy1,
            &mut xy2,
            &mut helix_len,
            &mut dxyz,
        )?;

        if a.abs() <= 0.00000001 {
            if a < 0.0 {
                a = -0.00000001;
            } else {
                a = 0.00000001;
            }
        }
        slope[i] = a;
        mm += 1;
        xy[mm as usize][1] = xy1[1];
        xy[mm as usize][2] = xy1[2];
        mm += 1;
        xy[mm as usize][1] = xy2[1];
        xy[mm as usize][2] = xy2[2];
        if (i as i64) == longest {
            helix_2d_vec[1] = xy2[1] - xy1[1];
            helix_2d_vec[2] = xy2[2] - xy1[2];
            helix_2d_vec[3] = 0.0;
        }
    }
    if mm != 2 * nhelix {
        // legacy prints warning
    }

    // rescale the helix axis points to PS format (1st time)
    xy4ps(1, &mut xy, mm, &mut nxy);
    let n1 = (longest * 2 - 1) as usize;
    let n2 = (longest * 2) as usize;
    let x1 = nxy[n1][1] - nxy[n2][1];
    let y1 = nxy[n1][2] - nxy[n2][2];
    let d1 = (x1 * x1 + y1 * y1).sqrt() / ((num_per_helix[longest as usize] - 1) as f64);

    let mut nh: i64 = 0;
    for i in 1..=nhelix as usize {
        let n = num_per_helix[i];
        nh += 1;
        let mut xy1 = [0.0f64; 4];
        xy1[1] = nxy[nh as usize][1];
        xy1[2] = nxy[nh as usize][2];
        nh += 1;
        let mut xy2 = [0.0f64; 4];
        xy2[1] = nxy[nh as usize][1];
        xy2[2] = nxy[nh as usize][2];

        let a = slope[i];
        let d2 = 0.5 * d1 * ((n - 1) as f64);
        new_xy(a, d2, &mut xy1, &mut xy2);
        gen_xy_cood(
            i as i64,
            num_residue,
            n,
            a,
            &xy1,
            &xy2,
            &bsp_1,
            &bsp_2,
            &mut nsub,
            &mut bs_sub,
            &loop_key,
            chain_id,
            seidx,
            &link,
            &bs1_lk,
            &bs2_lk,
            &mut xy_bso,
        );
    }

    rot_2D_to_Yaxis(num_residue, &mut helix_2d_vec, &mut xy_bso);

    // rescale all xy to PS format (2nd time)
    xy4ps(2, &mut xy_bso, num_residue, xy_bs);

    Ok(())
}

fn helix_regions(
    num_helixs: i64,
    helix_idx: &[[i64; 8]],
    bp_idx: &[i64],
    base_pairs: &[[i64; 18]],
    nhelix: &mut i64,
    npair_per_helix: &mut [i64],
    bs_1: &mut [Vec<i64>],
    bs_2: &mut [Vec<i64>],
) {
    let mut nh = 0i64;
    for i in 1..=num_helixs as usize {
        if helix_idx[i][2] <= helix_idx[i][1] {
            continue;
        }
        nh += 1;
        let mut n = 0i64;
        for j in helix_idx[i][1]..=helix_idx[i][2] {
            n += 1;
            let k = bp_idx[j as usize] as usize;
            bs_1[nh as usize][n as usize] = base_pairs[k][1];
            bs_2[nh as usize][n as usize] = base_pairs[k][2];
        }
        npair_per_helix[nh as usize] = n;
    }
    *nhelix = nh;
}

fn head_to_tail(
    j: i64,
    npair_per_helix: &mut [i64],
    bs_1: &mut [Vec<i64>],
    bs_2: &mut [Vec<i64>],
    nregions: &mut [i64],
    sub_helix1: &mut [Vec<i64>],
    sub_helix2: &mut [Vec<i64>],
) {
    let n = npair_per_helix[j as usize];
    if n <= 0 {
        nregions[j as usize] = 0;
        return;
    }
    let mut helix: Vec<[i64; 3]> = vec![[0; 3]; (n as usize) + 1];

    let mut nregion = 1i64;
    let mut i = 1i64;
    while i <= n {
        helix[nregion as usize][1] = i;
        i += 1;
        while i <= n
            && (((bs_1[j as usize][i as usize] == bs_1[j as usize][(i - 1) as usize] + 1)
                && (bs_2[j as usize][i as usize] == bs_2[j as usize][(i - 1) as usize] - 1))
                || ((bs_2[j as usize][i as usize] == bs_2[j as usize][(i - 1) as usize] + 1)
                    && (bs_1[j as usize][i as usize] == bs_1[j as usize][(i - 1) as usize] - 1)))
        {
            helix[nregion as usize][2] = i;
            i += 1;
        }
        i -= 1;
        if helix[nregion as usize][2] - helix[nregion as usize][1] >= 1 {
            nregion += 1;
        }
        i += 1;
    }

    nregions[j as usize] = nregion - 1;
    let mut nh = 1i64;
    for m in 1..=nregions[j as usize] as usize {
        sub_helix1[j as usize][m] = nh;
        for k in helix[m][1]..=helix[m][2] {
            let ku = k as usize;
            bs_1[j as usize][nh as usize] = bs_1[j as usize][ku];
            bs_2[j as usize][nh as usize] = bs_2[j as usize][ku];
            nh += 1;
        }
        sub_helix2[j as usize][m] = nh - 1;
    }
    npair_per_helix[j as usize] = nh - 1;
}

fn rest_bases(
    num_residue: i64,
    nhelix: i64,
    npair_per_helix: &[i64],
    bs_1: &[Vec<i64>],
    bs_2: &[Vec<i64>],
    nsub: &mut i64,
    bs_sub: &mut [i64],
) {
    let mut ns = 0i64;
    for k in 1..=num_residue {
        let mut ntest = 0i64;
        for i in 1..=nhelix as usize {
            ntest = 0;
            let n = npair_per_helix[i];
            for j in 1..=n as usize {
                if bs_1[i][j] == k || bs_2[i][j] == k {
                    ntest = 1;
                    break;
                }
            }
            if ntest == 1 {
                break;
            }
        }
        if ntest == 0 {
            ns += 1;
            bs_sub[ns as usize] = k;
        }
    }
    *nsub = ns;
}

fn helix_head(
    k: i64,
    _n: i64,
    bs_1: &[Vec<i64>],
    bs_2: &[Vec<i64>],
    nsub: &mut i64,
    bs_sub: &mut [i64],
    chain_id: &[u8],
    seidx: &[[usize; 3]],
    loop_out: &mut [i64; 3],
    yes: &mut i64,
) {
    *yes = 0;
    let (mut n1, mut n2) = (0i64, 0i64);
    if (bs_1[k as usize][1] > bs_1[k as usize][2] && bs_2[k as usize][1] < bs_2[k as usize][2])
        && (bs_2[k as usize][1] > bs_1[k as usize][1])
    {
        n1 = bs_1[k as usize][1];
        n2 = bs_2[k as usize][1];
        loop_proc(k, n1, n2, nsub, bs_sub, chain_id, seidx, yes);
    } else if (bs_2[k as usize][1] > bs_2[k as usize][2] && bs_1[k as usize][1] < bs_1[k as usize][2])
        && (bs_1[k as usize][1] > bs_2[k as usize][1])
    {
        n1 = bs_2[k as usize][1];
        n2 = bs_1[k as usize][1];
        loop_proc(k, n1, n2, nsub, bs_sub, chain_id, seidx, yes);
    }
    if *yes > 0 {
        if n1 > n2 {
            loop_out[2] = n1 - 1;
            loop_out[1] = n2 + 1;
        } else {
            loop_out[1] = n1 + 1;
            loop_out[2] = n2 - 1;
        }
    }
}

fn helix_tail(
    k: i64,
    n: i64,
    bs_1: &[Vec<i64>],
    bs_2: &[Vec<i64>],
    nsub: &mut i64,
    bs_sub: &mut [i64],
    chain_id: &[u8],
    seidx: &[[usize; 3]],
    loop_out: &mut [i64; 3],
    yes: &mut i64,
) {
    *yes = 0;
    let (mut n1, mut n2) = (0i64, 0i64);
    let nu = n as usize;

    if (bs_1[k as usize][nu] > bs_1[k as usize][nu - 1] && bs_2[k as usize][nu] < bs_2[k as usize][nu - 1])
        && (bs_2[k as usize][nu] > bs_1[k as usize][nu])
    {
        n1 = bs_1[k as usize][nu];
        n2 = bs_2[k as usize][nu];
        loop_proc(k, n1, n2, nsub, bs_sub, chain_id, seidx, yes);
    } else if (bs_2[k as usize][nu] > bs_2[k as usize][nu - 1] && bs_1[k as usize][nu] < bs_1[k as usize][nu - 1])
        && (bs_1[k as usize][nu] > bs_2[k as usize][nu])
    {
        n1 = bs_2[k as usize][nu];
        n2 = bs_1[k as usize][nu];
        loop_proc(k, n1, n2, nsub, bs_sub, chain_id, seidx, yes);
    }

    if *yes > 0 {
        if n1 > n2 {
            loop_out[2] = n1 - 1;
            loop_out[1] = n2 + 1;
        } else {
            loop_out[1] = n1 + 1;
            loop_out[2] = n2 - 1;
        }
    }
}

fn loop_proc(
    _k: i64,
    mut n1: i64,
    n2: i64,
    nsub: &mut i64,
    bs_sub: &mut [i64],
    chain_id: &[u8],
    seidx: &[[usize; 3]],
    yes: &mut i64,
) {
    *yes = 0;
    let k1 = seidx[n1 as usize][1];
    let k2 = seidx[n2 as usize][1];

    let mut n = 0i64;
    for j in (n1 + 1)..=(n2 - 1) {
        if n2 - n1 >= 20 || chain_id[k1] != chain_id[k2] {
            *yes = 0;
            return;
        }
        for i in 1..=*nsub as usize {
            if j == bs_sub[i] {
                n += 1;
            }
        }
    }

    if n == n2 - n1 - 1 && chain_id[k1] == chain_id[k2] {
        *yes = 1;
        loop {
            n1 += 1;
            if n1 > n2 {
                break;
            }
            let mut ns = 0i64;
            for i in 1..=*nsub as usize {
                if n1 != bs_sub[i] {
                    ns += 1;
                    bs_sub[ns as usize] = bs_sub[i];
                }
            }
            *nsub = ns;
            if n1 == n2 {
                break;
            }
            if n1 > n2 {
                break;
            }
        }
    } else {
        *yes = 0;
    }
}

fn check(m: i64, nsub: &mut i64, bs_sub: &mut [i64], yes: &mut i64) {
    let mut ns = 0i64;
    *yes = 0;
    for i in 1..=*nsub as usize {
        if m == bs_sub[i] {
            *yes = 1;
        } else {
            ns += 1;
            bs_sub[ns as usize] = bs_sub[i];
        }
    }
    *nsub = ns;
}

fn add_bs_2helix(
    i: i64,
    j: i64,
    n1: i64,
    n2: i64,
    num_residue: i64,
    bs_1: &[Vec<i64>],
    bs_2: &[Vec<i64>],
    nsub: &mut i64,
    bs_sub: &mut [i64],
    add: &mut [i64],
    bs_1_add: &mut [i64],
    bs_2_add: &mut [i64],
) {
    let mut nn1 = bs_1[i as usize][n1 as usize];
    let mut nn2 = bs_2[i as usize][n1 as usize];
    let mut tmp1 = 0i64;
    let mut tmp2 = 0i64;

    if bs_1[i as usize][n1 as usize] > bs_1[i as usize][(n1 - 1) as usize] {
        // increasing
        for _m in 1..=num_residue {
            nn1 += 1;
            nn2 -= 1;

            let mut yes = 0i64;
            check(nn1, nsub, bs_sub, &mut yes);
            if yes > 0 {
                add[j as usize] += 1;
                bs_1_add[add[j as usize] as usize] = nn1;
                bs_2_add[add[j as usize] as usize] = 0;
            } else {
                nn1 = -num_residue;
                tmp1 = 1;
            }

            yes = 0;
            check(nn2, nsub, bs_sub, &mut yes);
            if yes > 0 {
                add[j as usize] += 1;
                bs_1_add[add[j as usize] as usize] = 0;
                bs_2_add[add[j as usize] as usize] = nn2;
            } else {
                tmp2 = 1;
                nn2 = -num_residue;
            }

            if nn1 > bs_1[i as usize][n2 as usize] {
                tmp1 = 1;
            }
            if nn2 < bs_2[i as usize][n2 as usize] {
                tmp2 = 1;
            }
            if tmp1 == 1 && tmp2 == 1 {
                break;
            }
        }
    } else {
        // decreasing
        for _m in 1..=num_residue {
            nn1 -= 1;
            nn2 += 1;

            let mut yes = 0i64;
            check(nn1, nsub, bs_sub, &mut yes);
            if yes > 0 {
                add[j as usize] += 1;
                bs_1_add[add[j as usize] as usize] = nn1;
                bs_2_add[add[j as usize] as usize] = 0;
            } else {
                tmp1 = 1;
                nn1 = -num_residue;
            }

            yes = 0;
            check(nn2, nsub, bs_sub, &mut yes);
            if yes > 0 {
                add[j as usize] += 1;
                bs_1_add[add[j as usize] as usize] = 0;
                bs_2_add[add[j as usize] as usize] = nn2;
            } else {
                tmp2 = 1;
                nn2 = -num_residue;
            }

            if nn1 < bs_1[i as usize][n2 as usize] {
                tmp1 = 1;
            }
            if nn2 > bs_2[i as usize][n2 as usize] {
                tmp2 = 1;
            }
            if tmp1 == 1 && tmp2 == 1 {
                break;
            }
        }
    }
}

fn chck_lk(diff: i64, m1: i64, m2: i64, nsub: i64, bs_sub: &[i64]) -> i64 {
    let mut n = 0i64;
    for m in m1..=m2 {
        for i in 1..=nsub as usize {
            if m == bs_sub[i] {
                n += 1;
            }
        }
    }
    if n < diff {
        0
    } else {
        1
    }
}

fn check_link(i: i64, n1: i64, n2: i64, nsub: i64, bs_sub: &[i64], bs_1: &[Vec<i64>], bs_2: &[Vec<i64>], yes: &mut i64) {
    *yes = 1;
    let (mut m1, mut m2) = if bs_1[i as usize][n1 as usize] > bs_1[i as usize][n2 as usize] {
        (bs_1[i as usize][n2 as usize], bs_1[i as usize][n1 as usize])
    } else {
        (bs_1[i as usize][n1 as usize], bs_1[i as usize][n2 as usize])
    };
    m1 += 1;
    m2 -= 1;
    let diff = m2 - m1 + 1;
    let yes1 = chck_lk(diff, m1, m2, nsub, bs_sub);

    let (mut m1b, mut m2b) = if bs_2[i as usize][n1 as usize] > bs_2[i as usize][n2 as usize] {
        (bs_2[i as usize][n2 as usize], bs_2[i as usize][n1 as usize])
    } else {
        (bs_2[i as usize][n1 as usize], bs_2[i as usize][n2 as usize])
    };
    m1b += 1;
    m2b -= 1;
    let diff2 = m2b - m1b + 1;
    let yes2 = chck_lk(diff2, m1b, m2b, nsub, bs_sub);

    if yes1 == 0 || yes2 == 0 {
        *yes = 0;
    }
}

fn xy_base(j: i64, n: i64, n1: i64, n01: i64, a: f64, d1: f64, xy1: &[f64; 4], xy2: &[f64; 4], xy_bs: &mut [[f64; 4]]) {
    let a0 = d1 / (1.0 + a * a).sqrt();
    let b0 = a * d1 / (1.0 + a * a).sqrt();
    let x0i = xy_bs[n01 as usize][1];
    let y0i = xy_bs[n01 as usize][2];

    let xi = xy1[1];
    let yi = xy1[2];
    let xf = xy2[1];
    let yf = xy2[2];

    let (mut x0, mut y0);
    if n == 1 {
        x0 = x0i - a0 * (j as f64);
        y0 = y0i - b0 * (j as f64);
        let c1 = (y0 - y0i) * (yi - yf) + (x0 - x0i) * (xi - xf);
        let sign = c1 / c1.abs();
        x0 = x0i - sign * a0 * (j as f64);
        y0 = y0i - sign * b0 * (j as f64);
    } else {
        x0 = x0i + a0 * (j as f64);
        y0 = y0i + b0 * (j as f64);
        let c1 = (y0 - y0i) * (yf - yi) + (x0 - x0i) * (xf - xi);
        let sign = c1 / c1.abs();
        x0 = x0i + sign * a0 * (j as f64);
        y0 = y0i + sign * b0 * (j as f64);
    }

    xy_bs[n1 as usize][1] = x0;
    xy_bs[n1 as usize][2] = y0;
}

fn link_helix(
    n: i64,
    mut n1: i64,
    mut n2: i64,
    num_residue: i64,
    mut d1: f64,
    _d2: f64,
    a: f64,
    chain_id: &[u8],
    seidx: &[[usize; 3]],
    xy1: &[f64; 4],
    xy2: &[f64; 4],
    nsub: &mut i64,
    bs_sub: &mut [i64],
    xy_bs: &mut [[f64; 4]],
) {
    d1 = 0.6 * d1;
    let n01 = n1;
    let n02 = n2;
    let mut nn1 = n1;
    let mut nn2 = n2;
    if n1 == 0 || n2 == 0 {
        return;
    }

    let mut yes = 0i64;
    let mut add = 0i64;
    for _m in 1..=50 {
        n1 += 1;
        if n1 >= num_residue {
            break;
        }
        let k1 = seidx[n01 as usize][1];
        let k2 = seidx[n1 as usize][1];
        if chain_id[k1] == chain_id[k2] {
            check(n1, nsub, bs_sub, &mut yes);
            if yes > 0 {
                add += 1;
                xy_base(add, n, n1, n01, a, d1, xy1, xy2, xy_bs);
            } else {
                break;
            }
        }
    }
    yes = 0;
    add = 0;
    for _m in 1..=50 {
        nn1 -= 1;
        if nn1 <= 0 {
            break;
        }
        let k1 = seidx[n01 as usize][1];
        let k2 = seidx[nn1 as usize][1];
        if chain_id[k1] == chain_id[k2] {
            check(nn1, nsub, bs_sub, &mut yes);
            if yes > 0 {
                add += 1;
                xy_base(add, n, nn1, n01, a, d1, xy1, xy2, xy_bs);
            } else {
                break;
            }
        }
    }
    yes = 0;
    add = 0;
    for _m in 1..=50 {
        n2 -= 1;
        if n2 <= 0 {
            break;
        }
        let k1 = seidx[n02 as usize][1];
        let k2 = seidx[n2 as usize][1];
        if chain_id[k1] == chain_id[k2] {
            check(n2, nsub, bs_sub, &mut yes);
            if yes > 0 {
                add += 1;
                xy_base(add, n, n2, n02, a, d1, xy1, xy2, xy_bs);
            } else {
                break;
            }
        }
    }
    yes = 0;
    add = 0;
    for _m in 1..=50 {
        nn2 += 1;
        if nn2 > num_residue {
            break;
        }
        let k1 = seidx[n02 as usize][1];
        let k2 = seidx[nn2 as usize][1];
        if chain_id[k1] == chain_id[k2] {
            check(nn2, nsub, bs_sub, &mut yes);
            if yes > 0 {
                add += 1;
                xy_base(add, n, nn2, n02, a, d1, xy1, xy2, xy_bs);
            } else {
                break;
            }
        }
    }
}

fn gen_xy_cood(
    i: i64,
    num_residue: i64,
    n: i64,
    a: f64,
    xy1: &[f64; 4],
    xy2: &[f64; 4],
    bs_1: &[Vec<i64>],
    bs_2: &[Vec<i64>],
    nsub: &mut i64,
    bs_sub: &mut [i64],
    loop_key: &[[i64; 3]],
    chain_id: &[u8],
    seidx: &[[usize; 3]],
    link: &[i64],
    bs1_lk: &[Vec<i64>],
    bs2_lk: &[Vec<i64>],
    xy_bs: &mut [[f64; 4]],
) {
    let xi = xy1[1];
    let yi = xy1[2];
    let xf = xy2[1];
    let yf = xy2[2];
    let at = -1.0 / a;
    let len = ((xf - xi) * (xf - xi) + (yf - yi) * (yf - yi)).sqrt();
    let d1 = len / ((n - 1) as f64);
    let d2 = 1.2 * d1;

    let a0 = d1 / (1.0 + a * a).sqrt();
    let b0 = a * d1 / (1.0 + a * a).sqrt();
    let a1 = d2 / (1.0 + at * at).sqrt();
    let b1 = at * d2 / (1.0 + at * at).sqrt();

    for j in 1..=n {
        let mut x12 = [0.0f64; 3];
        let mut y12 = [0.0f64; 3];
        xy_at_base12(a0, b0, a1, b1, j, xi, yi, xf, yf, &mut x12, &mut y12);

        let k = bs_1[i as usize][j as usize] as usize;
        xy_bs[k][1] = x12[1];
        xy_bs[k][2] = y12[1];

        let k = bs_2[i as usize][j as usize] as usize;
        xy_bs[k][1] = x12[2];
        xy_bs[k][2] = y12[2];
    }

    if loop_key[i as usize][1] > 0 {
        // for the head of helix
        loop_xy(i, 1, bs_1, bs_2, xy1, xy2, a, xy_bs);
    } else {
        let n1 = bs_1[i as usize][1];
        let n2 = bs_2[i as usize][1];
        link_helix(1, n1, n2, num_residue, d1, d2, a, chain_id, seidx, xy1, xy2, nsub, bs_sub, xy_bs);
    }

    if loop_key[i as usize][2] > 0 {
        // for the tail of helix
        loop_xy(i, n, bs_1, bs_2, xy2, xy1, a, xy_bs);
    } else {
        let n1 = bs_1[i as usize][n as usize];
        let n2 = bs_2[i as usize][n as usize];
        link_helix(2, n1, n2, num_residue, d1, d2, a, chain_id, seidx, xy1, xy2, nsub, bs_sub, xy_bs);
    }

    let mut num = 0i64;
    for k in 1..=link[i as usize] {
        link_xy(i, k, d2, nsub, bs_sub, bs1_lk, bs2_lk, xy_bs, &mut num);
    }
}

fn link_xy(i: i64, j: i64, d: f64, nsub: &mut i64, bs_sub: &mut [i64], bs_1: &[Vec<i64>], bs_2: &[Vec<i64>], xy_bs: &mut [[f64; 4]], num: &mut i64) {
    let k01 = bs_1[i as usize][j as usize];
    let k02 = bs_2[i as usize][j as usize];

    let mut yes = 0i64;
    let mut k1 = k01;
    let mut k2 = k02;
    let mut m = 0i64;

    // increasing k1
    k1 = k01;
    k2 = k02;
    m = 0;
    loop {
        k1 += 1;
        check(k1, nsub, bs_sub, &mut yes);
        if yes > 0 {
            m += 1;
            link_xy_proc(m, d, k01, k02, k1, xy_bs);
        } else {
            break;
        }
        if k1 >= k01 + 20 {
            break;
        }
    }

    // decreasing k1
    k1 = k01;
    k2 = k02;
    m = 0;
    loop {
        k1 -= 1;
        check(k1, nsub, bs_sub, &mut yes);
        if yes > 0 {
            m += 1;
            link_xy_proc(m, d, k01, k02, k1, xy_bs);
        } else {
            break;
        }
        if k1 <= 0 {
            break;
        }
    }

    // increasing k2
    k1 = k01;
    k2 = k02;
    m = 0;
    loop {
        k2 += 1;
        check(k2, nsub, bs_sub, &mut yes);
        if yes > 0 {
            m += 1;
            link_xy_proc(m, d, k02, k01, k2, xy_bs);
        } else {
            break;
        }
        if k2 >= k02 + 20 {
            break;
        }
    }

    // decreasing k2
    k1 = k01;
    k2 = k02;
    m = 0;
    loop {
        k2 -= 1;
        check(k2, nsub, bs_sub, &mut yes);
        if yes > 0 {
            m += 1;
            link_xy_proc(m, d, k02, k01, k2, xy_bs);
        } else {
            k2 += 1;
            break;
        }
        if k2 <= 0 {
            break;
        }
    }

    // legacy signature includes this output but doesn't assign it.
    let _ = (k1, k2, bs_1, bs_2);
    *num = 0;
}

fn link_xy_proc(m: i64, d: f64, k01: i64, k02: i64, k1: i64, xy_bs: &mut [[f64; 4]]) {
    let d1 = 0.5 * d;
    let xf = xy_bs[k01 as usize][1];
    let yf = xy_bs[k01 as usize][2];
    let xi = xy_bs[k02 as usize][1];
    let yi = xy_bs[k02 as usize][2];
    let mut dn = xf - xi;
    if dn.abs() <= 0.0000001 {
        if dn < 0.0 {
            dn = -0.0000001;
        } else {
            dn = 0.0000001;
        }
    }

    let a = (yf - yi) / dn;
    let a0 = d1 / (1.0 + a * a).sqrt();
    let b0 = a * d1 / (1.0 + a * a).sqrt();
    let mut x0 = xf + a0 * (m as f64);
    let mut y0 = yf + b0 * (m as f64);
    let sign = (y0 - yf) * (yf - yi) + (x0 - xf) * (xf - xi);
    if sign >= 0.0 {
        x0 = xf + a0 * (m as f64);
        y0 = yf + b0 * (m as f64);
    } else {
        x0 = xf - a0 * (m as f64);
        y0 = yf - b0 * (m as f64);
    }
    xy_bs[k1 as usize][1] = x0;
    xy_bs[k1 as usize][2] = y0;
}

fn loop_xy(i: i64, n: i64, bs_1: &[Vec<i64>], bs_2: &[Vec<i64>], xy1: &[f64; 4], xy2: &[f64; 4], a: f64, xy_bs: &mut [[f64; 4]]) {
    let x01 = xy1[1];
    let y01 = xy1[2];
    let x02 = xy2[1];
    let y02 = xy2[2];

    let n1 = bs_1[i as usize][n as usize];
    let n2 = bs_2[i as usize][n as usize];
    let x1 = xy_bs[n1 as usize][1];
    let y1 = xy_bs[n1 as usize][2];

    let m = (n1 - n2).abs() - 1;
    let mut ang = 90.0 - (180.0 / std::f64::consts::PI) * a.atan();
    let d = ((y1 - y01) * (y1 - y01) + (x1 - x01) * (x1 - x01)).sqrt();
    let h = 0.15 * (m as f64) * d;
    let r = (d * d + h * h).sqrt();
    let a0 = h / (1.0 + a * a).sqrt();
    let b0 = a * h / (1.0 + a * a).sqrt();

    let mut x0 = x01 + a0;
    let mut y0 = y01 + b0;
    let sign = (y0 - y01) * (y02 - y01) + (x0 - x01) * (x02 - x01);
    if sign >= 0.0 {
        x0 = x01 - a0;
        y0 = y01 - b0;
    } else {
        x0 = x01 + a0;
        y0 = y01 + b0;
    }

    let ap = (180.0 / std::f64::consts::PI) * (h / r).asin();
    let alfa = (180.0 + 2.0 * ap) / ((1 + m) as f64);

    let c01 = x0 - x01;
    let c02 = y0 - y01;
    let c11 = x1 - x01;
    let c12 = y1 - y01;
    if c01 < 0.0 {
        ang = 180.0 + ang;
    }

    let sign = c11 * c02 - c12 * c01;
    if sign >= 0.0 {
        loop_xy_proc(i, n, m, alfa, ang, ap, bs_1, r, x0, y0, xy_bs);
    } else {
        loop_xy_proc(i, n, m, alfa, ang, ap, bs_2, r, x0, y0, xy_bs);
    }
}

fn loop_xy_proc(i: i64, n: i64, m: i64, alfa: f64, ang: f64, ap: f64, bs: &[Vec<i64>], r: f64, x0: f64, y0: f64, xy_bs: &mut [[f64; 4]]) {
    if n > 1 {
        // at the end
        if bs[i as usize][n as usize] < bs[i as usize][(n - 1) as usize] {
            // decreasing
            for j in 1..=m {
                let alf = (j as f64) * alfa;
                let a2 = alf - ap;
                let angl = (ang - a2) * std::f64::consts::PI / 180.0;
                let k = bs[i as usize][n as usize] - j;
                xy_bs[k as usize][1] = x0 + r * angl.cos();
                xy_bs[k as usize][2] = y0 - r * angl.sin();
            }
        } else {
            for j in 1..=m {
                let alf = (j as f64) * alfa;
                let a2 = alf - ap;
                let angl = (ang - a2) * std::f64::consts::PI / 180.0;
                let k = bs[i as usize][n as usize] + j;
                xy_bs[k as usize][1] = x0 + r * angl.cos();
                xy_bs[k as usize][2] = y0 - r * angl.sin();
            }
        }
    } else {
        // at the beginning
        if bs[i as usize][n as usize] < bs[i as usize][(n + 1) as usize] {
            // decreasing
            for j in 1..=m {
                let alf = (j as f64) * alfa;
                let a2 = alf - ap;
                let angl = (ang - a2) * std::f64::consts::PI / 180.0;
                let k = bs[i as usize][n as usize] - j;
                xy_bs[k as usize][1] = x0 + r * angl.cos();
                xy_bs[k as usize][2] = y0 - r * angl.sin();
            }
        } else {
            for j in 1..=m {
                let alf = (j as f64) * alfa;
                let a2 = alf - ap;
                let angl = (ang - a2) * std::f64::consts::PI / 180.0;
                let k = bs[i as usize][n as usize] + j;
                xy_bs[k as usize][1] = x0 + r * angl.cos();
                xy_bs[k as usize][2] = y0 - r * angl.sin();
            }
        }
    }
}

fn xy_at_base12(a0: f64, b0: f64, a1: f64, b1: f64, j: i64, xi: f64, yi: f64, xf: f64, yf: f64, x12: &mut [f64; 3], y12: &mut [f64; 3]) {
    let mut x0 = xi + a0 * ((j - 1) as f64);
    let mut y0 = yi + b0 * ((j - 1) as f64);
    let sign = ((yf - y0) * (y0 - yi) + (xf - x0) * (x0 - xi)) as i64;
    if sign >= 0 {
        x0 = xi + a0 * ((j - 1) as f64);
        y0 = yi + b0 * ((j - 1) as f64);
    } else {
        x0 = xi - a0 * ((j - 1) as f64);
        y0 = yi - b0 * ((j - 1) as f64);
    }
    x12[1] = x0 + a1;
    y12[1] = y0 + b1;
    x12[2] = x0 - a1;
    y12[2] = y0 - b1;
}

fn new_xy(a: f64, d: f64, xy1: &mut [f64; 4], xy2: &mut [f64; 4]) {
    let x0 = 0.5 * (xy1[1] + xy2[1]);
    let y0 = 0.5 * (xy1[2] + xy2[2]);

    let a0 = d / (1.0 + a * a).sqrt();
    let b0 = a * d / (1.0 + a * a).sqrt();

    let xi = xy1[1];
    let yi = xy1[2];
    let xf = xy2[1];
    let yf = xy2[2];

    let mut x00 = x0 - a0;
    let mut y00 = y0 - b0;
    let c1 = (y00 - y0) * (yf - yi) + (x00 - x0) * (xf - xi);
    if c1 >= 0.0 {
        xy2[1] = x00;
        xy2[2] = y00;
    } else {
        xy1[1] = x00;
        xy1[2] = y00;
    }

    x00 = x0 + a0;
    y00 = y0 + b0;
    let c1 = (y00 - y0) * (yf - yi) + (x00 - x0) * (xf - xi);
    if c1 >= 0.0 {
        xy2[1] = x00;
        xy2[2] = y00;
    } else {
        xy1[1] = x00;
        xy1[2] = y00;
    }
}

fn max_dmatrix(oxy: &[[f64; 4]], num: i64, nc: usize, max_xy: &mut [f64; 3]) {
    max_xy[1] = -XBIG;
    max_xy[2] = -XBIG;
    for i in 1..=num as usize {
        for j in 1..=nc {
            if oxy[i][j] > max_xy[j] {
                max_xy[j] = oxy[i][j];
            }
        }
    }
}

fn min_dmatrix(oxy: &[[f64; 4]], num: i64, nc: usize, min_xy: &mut [f64; 3]) {
    min_xy[1] = XBIG;
    min_xy[2] = XBIG;
    for i in 1..=num as usize {
        for j in 1..=nc {
            if oxy[i][j] < min_xy[j] {
                min_xy[j] = oxy[i][j];
            }
        }
    }
}

fn move_position(oxy: &mut [[f64; 4]], num: i64, nc: usize, min_xy: &[f64; 3]) {
    for i in 1..=num as usize {
        for j in 1..=nc {
            oxy[i][j] -= min_xy[j];
        }
    }
}

fn dmax(a: f64, b: f64) -> f64 {
    if a > b { a } else { b }
}

fn xy4ps(n: i64, oxy: &mut [[f64; 4]], num: i64, nxy: &mut [[f64; 4]]) {
    let default_size = 550.0f64;
    let mut max_xy = [0.0f64; 3];
    let mut min_xy = [0.0f64; 3];

    max_dmatrix(oxy, num, 2, &mut max_xy);
    min_dmatrix(oxy, num, 2, &mut min_xy);
    move_position(oxy, num, 2, &min_xy);

    let temp = dmax(max_xy[1] - min_xy[1], max_xy[2] - min_xy[2]);
    let scale_factor = default_size / temp.abs();
    for i in 1..=num as usize {
        for j in 1..=2 {
            nxy[i][j] = oxy[i][j] * scale_factor;
        }
    }

    if n > 1 {
        max_dmatrix(nxy, num, 2, &mut max_xy);
    }
}

fn rot_2D_to_Yaxis(num_residue: i64, z: &mut [f64; 4], xy_bso: &mut [[f64; 4]]) {
    let mut yphy = [0.0f64; 4];
    yphy[0] = EMPTY_NUMBER;
    yphy[1] = 0.0;
    yphy[2] = 1.0;
    yphy[3] = 0.0;

    let hinge = cross(z, &yphy);
    let rotmat = arb_rotation(hinge, magang(yphy, *z));

    let mut xyz_tmp: Vec<[f64; 4]> = vec![[0.0; 4]; (num_residue as usize) + 1];
    for i in 1..=num_residue as usize {
        xyz_tmp[i][1] = xy_bso[i][1];
        xyz_tmp[i][2] = xy_bso[i][2];
        xyz_tmp[i][3] = 0.0;
    }

    for i in 1..=num_residue as usize {
        let mut out = [0.0f64; 4];
        for k in 1..=3 {
            out[k] = dot(&xyz_tmp[i], &rotmat[k]);
        }
        xy_bso[i][1] = out[1];
        xy_bso[i][2] = out[2];
    }
}

fn HelixAxis(
    ii: i64,
    n: i64,
    bs_1: &[Vec<i64>],
    bs_2: &[Vec<i64>],
    num_residue: i64,
    seidx: &[[usize; 3]],
    ry: &[i32],
    xyz: &[[f64; 4]],
    atom_name: &[AtomName4],
    bseq: &[u8],
    a: &mut f64,
    xy1: &mut [f64; 4],
    xy2: &mut [f64; 4],
    helix_len: &mut [f64],
    dxyz: &mut [[f64; 4]],
) -> Result<(), String> {
    let mut nxyz: Vec<[f64; 4]> = vec![[0.0; 4]; (num_residue as usize) * 100 + 2];
    let mut chi: Vec<Vec<usize>> = vec![vec![0usize; (n as usize) * 4 + 2]; 4]; // [1..3][1..n*4]
    let mut residx: Vec<[usize; 3]> = vec![[0usize; 3]; (num_residue as usize) + 2];
    let mut atmnam: Vec<AtomName4> = vec![AtomName4([0; 4]); (num_residue as usize) * 100 + 2];

    let mut atmnum: usize = 1;
    for i in 1..=n as usize {
        let n1 = bs_1[ii as usize][i] as usize;
        residx[n1][1] = atmnum;
        for j in seidx[n1][1]..=seidx[n1][2] {
            atmnam[atmnum] = atom_name[j];
            for k in 1..=3 {
                nxyz[atmnum][k] = xyz[j][k];
            }
            atmnum += 1;
        }
        residx[n1][2] = atmnum - 1;

        let n2 = bs_2[ii as usize][i] as usize;
        residx[n2][1] = atmnum;
        for j in seidx[n2][1]..=seidx[n2][2] {
            atmnam[atmnum] = atom_name[j];
            for k in 1..=3 {
                nxyz[atmnum][k] = xyz[j][k];
            }
            atmnum += 1;
        }
        residx[n2][2] = atmnum - 1;
    }
    let atmnum_total = atmnum - 1;

    get_chi(1, ii, n, bs_1, &residx, &atmnam, bseq, ry, &mut chi)?;
    get_chi(2, ii, n, bs_2, &residx, &atmnam, bseq, ry, &mut chi)?;

    let mut hstart = [0.0f64; 4];
    let mut hend = [0.0f64; 4];
    axis_start_end(n, atmnum_total as i64, &chi, &nxyz, &mut hstart, &mut hend)?;

    let dx1 = hstart[1] - hend[1];
    let dx2 = hstart[2] - hend[2];
    let dx3 = hstart[3] - hend[3];
    helix_len[ii as usize] = (dx1 * dx1 + dx2 * dx2 + dx3 * dx3).sqrt();
    dxyz[ii as usize][1] = -dx1 / helix_len[ii as usize];
    dxyz[ii as usize][2] = -dx2 / helix_len[ii as usize];
    dxyz[ii as usize][3] = -dx3 / helix_len[ii as usize];

    xy1[1] = hstart[1];
    xy1[2] = hstart[2];
    xy2[1] = hend[1];
    xy2[2] = hend[2];

    let mut div = hend[1] - hstart[1];
    if div.abs() <= 0.00000001 {
        if div <= 0.0 {
            div = -0.00000001;
        } else {
            div = 0.00000001;
        }
    }
    *a = (hend[2] - hstart[2]) / div;
    Ok(())
}

fn get_chi(
    which: i64,
    ii: i64,
    n: i64,
    bs: &[Vec<i64>],
    residx: &[[usize; 3]],
    atmnam: &[AtomName4],
    bseq: &[u8],
    ry: &[i32],
    chi: &mut [Vec<usize>],
) -> Result<(), String> {
    let o4 = AtomName4::from_bytes(*b" O4'");
    let c1 = AtomName4::from_bytes(*b" C1'");
    let n9 = AtomName4::from_bytes(*b" N9 ");
    let c4 = AtomName4::from_bytes(*b" C4 ");
    let n1 = AtomName4::from_bytes(*b" N1 ");
    let c2 = AtomName4::from_bytes(*b" C2 ");
    let c5 = AtomName4::from_bytes(*b" C5 ");

    let mut offset = 0usize;
    for j in 1..=n as usize {
        let idx = bs[ii as usize][j] as usize;
        let ib = residx[idx][1];
        let ie = residx[idx][2];
        let o4_idx = find_1st_atom(o4, atmnam, ib, ie);
        let c1_idx = find_1st_atom(c1, atmnam, ib, ie);

        let mut n_atom = n1;
        let mut c_atom = c2;
        if ry.get(idx).copied().unwrap_or(0) == 1 {
            n_atom = n9;
            c_atom = c4;
        } else if matches!(bseq.get(idx).copied().unwrap_or(b'?'), b'P' | b'p') {
            n_atom = c5;
            c_atom = c4;
        }
        let n1n9_idx = find_1st_atom(n_atom, atmnam, ib, ie);
        let c2c4_idx = find_1st_atom(c_atom, atmnam, ib, ie);

        offset = (j - 1) * 4;
        chi[which as usize][offset + 1] = o4_idx.unwrap_or(0);
        chi[which as usize][offset + 2] = c1_idx.unwrap_or(0);
        chi[which as usize][offset + 3] = n1n9_idx.unwrap_or(0);
        chi[which as usize][offset + 4] = c2c4_idx.unwrap_or(0);
    }
    Ok(())
}

fn axis_start_end(
    num_bp: i64,
    num: i64,
    chi: &[Vec<usize>],
    xyz: &[[f64; 4]],
    hstart: &mut [f64; 4],
    hend: &mut [f64; 4],
) -> Result<(), String> {
    let nbpm1 = num_bp - 1;
    if nbpm1 <= 0 {
        return Ok(());
    }

    // collect vectors
    let mut idx: Vec<[usize; 3]> = vec![[0usize; 3]; (4 * nbpm1) as usize + 2];
    let mut nvec = 0usize;
    for i in 1..=2usize {
        for j in 1..=nbpm1 as usize {
            let ioffset = (j - 1) * 4;
            let joffset = j * 4;
            for k in 2..=3usize {
                let ia = chi[i][ioffset + k];
                let ib = chi[i][joffset + k];
                if ia != 0 && ib != 0 {
                    nvec += 1;
                    idx[nvec][1] = ia;
                    idx[nvec][2] = ib;
                }
            }
        }
    }
    if nvec < 3 {
        return Ok(());
    }

    let mut vxyz: Vec<[f64; 4]> = vec![[0.0; 4]; nvec + 1];
    let mut drise: Vec<f64> = vec![0.0; nvec + 1];
    for i in 1..=nvec {
        for j in 1..=3 {
            vxyz[i][j] = xyz[idx[i][2]][j] - xyz[idx[i][1]][j];
        }
    }
    let mut haxis = [0.0f64; 4];
    let mut hinge = [0.0f64; 4];
    let mut hrise = 0.0f64;
    ls_plane(&vxyz, nvec as i64, &mut haxis, &mut hinge, &mut hrise, &mut drise)?;
    if hrise < 0.0 {
        hrise = -hrise;
        for i in 1..=3 {
            haxis[i] = -haxis[i];
        }
    }

    // align haxis to global z-axis
    let mut z = [0.0f64; 4];
    z[0] = EMPTY_NUMBER;
    z[3] = 1.0;
    let hinge2 = cross(&haxis, &z);
    let ang_deg = magang(haxis, z);
    let rotmat = arb_rotation(hinge2, ang_deg);

    // xyzH = xyz * rotmatT (legacy uses transpose + multi_matrix)
    let mut xyzH: Vec<[f64; 4]> = vec![[0.0; 4]; (num as usize) + 1];
    for i in 1..=num as usize {
        for k in 1..=3 {
            xyzH[i][k] = dot(&xyz[i], &rotmat[k]);
        }
    }

    // locate xy-coordinate the helix passes through
    let mut dxy: Vec<[f64; 3]> = vec![[0.0; 3]; nvec + 1];
    let mut g: Vec<f64> = vec![0.0; nvec + 1];
    for i in 1..=nvec {
        g[i] = 0.0;
        for j in 1..=2 {
            let tb = xyzH[idx[i][1]][j];
            let te = xyzH[idx[i][2]][j];
            dxy[i][j] = 2.0 * (te - tb);
            g[i] += te * te - tb * tb;
        }
    }

    // t2 = g * dxy (1xnvec * nvecx2) -> 2
    let mut t2 = [0.0f64; 3];
    for j in 1..=2 {
        let mut s = 0.0;
        for i in 1..=nvec {
            s += g[i] * dxy[i][j];
        }
        t2[j] = s;
    }

    // dd = dxyT * dxy (2xnvec * nvecx2) -> 2x2
    let mut dd = [[0.0f64; 3]; 3];
    for i in 1..=2 {
        for j in 1..=2 {
            let mut s = 0.0;
            for k in 1..=nvec {
                s += dxy[k][i] * dxy[k][j];
            }
            dd[i][j] = s;
        }
    }
    let inv_dd = dinverse_2(&dd)?;
    let mut org_xyz = [0.0f64; 4];
    for i in 1..=2 {
        // Legacy multi_vec_matrix: o[i] = sum_j a[j] * b[j][i]
        org_xyz[i] = t2[1] * inv_dd[1][i] + t2[2] * inv_dd[2][i];
    }

    org_xyz[3] = 0.5 * (xyzH[chi[1][2]][3] + xyzH[chi[2][2]][3]);
    multi_vec_matrix(&org_xyz, &rotmat, hstart);

    let ioffset = (nbpm1 as usize) * 4;
    org_xyz[3] = 0.5 * (xyzH[chi[1][ioffset + 2]][3] + xyzH[chi[2][ioffset + 2]][3]);
    multi_vec_matrix(&org_xyz, &rotmat, hend);

    Ok(())
}

fn dinverse_2(dd: &[[f64; 3]; 3]) -> Result<[[f64; 3]; 3], String> {
    let det = dd[1][1] * dd[2][2] - dd[1][2] * dd[2][1];
    if det.abs() < XEPS {
        return Err("dinverse: singular 2x2".to_string());
    }
    let mut inv = [[0.0f64; 3]; 3];
    inv[1][1] = dd[2][2] / det;
    inv[1][2] = -dd[1][2] / det;
    inv[2][1] = -dd[2][1] / det;
    inv[2][2] = dd[1][1] / det;
    Ok(inv)
}

fn multi_vec_matrix(v: &[f64; 4], m: &[[f64; 4]; 4], out: &mut [f64; 4]) {
    // Legacy multi_vec_matrix: o[i] = sum_j a[j] * b[j][i]
    for i in 1..=3 {
        let mut s = 0.0f64;
        for j in 1..=3 {
            s += v[j] * m[j][i];
        }
        out[i] = s;
    }
}

fn cov_matrix(points: &[[f64; 4]], n: i64, cov: &mut [[f64; 4]; 4]) {
    let mut ave = [0.0f64; 4];
    for i in 1..=n as usize {
        for k in 1..=3 {
            ave[k] += points[i][k];
        }
    }
    for k in 1..=3 {
        ave[k] /= n as f64;
    }

    for i in 1..=3 {
        for j in 1..=3 {
            let mut s = 0.0;
            for k in 1..=n as usize {
                s += (points[k][i] - ave[i]) * (points[k][j] - ave[j]);
            }
            cov[i][j] = s / (n as f64);
        }
    }
}

fn jacobi(a: &mut [[f64; 4]; 4], n: usize, d: &mut [f64; 4], v: &mut [[f64; 4]; 4]) {
    // Minimal Jacobi for 3x3, ported from legacy_alg.rs style.
    let mut b = [0.0f64; 4];
    let mut z = [0.0f64; 4];
    for ip in 1..=n {
        for iq in 1..=n {
            v[ip][iq] = if ip == iq { 1.0 } else { 0.0 };
        }
        b[ip] = a[ip][ip];
        d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }
    for iter in 1..=100 {
        let mut sm = 0.0;
        for ip in 1..=n - 1 {
            for iq in ip + 1..=n {
                sm += a[ip][iq].abs();
            }
        }
        if sm < XEPS {
            return;
        }
        let tresh = if iter < 4 { 0.2 * sm / ((n * n) as f64) } else { 0.0 };
        for ip in 1..=n - 1 {
            for iq in ip + 1..=n {
                let g = 100.0 * a[ip][iq].abs();
                if iter > 4 && (d[ip].abs() + g) == d[ip].abs() && (d[iq].abs() + g) == d[iq].abs() {
                    a[ip][iq] = 0.0;
                } else if a[ip][iq].abs() > tresh {
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
                    for j in 1..=ip - 1 {
                        rotate(a, j, ip, j, iq, s, tau);
                    }
                    for j in ip + 1..=iq - 1 {
                        rotate(a, ip, j, j, iq, s, tau);
                    }
                    for j in iq + 1..=n {
                        rotate(a, ip, j, iq, j, s, tau);
                    }
                    for j in 1..=n {
                        rotate(v, j, ip, j, iq, s, tau);
                    }
                }
            }
        }
        for ip in 1..=n {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
}

fn rotate(a: &mut [[f64; 4]; 4], i: usize, j: usize, k: usize, l: usize, s: f64, tau: f64) {
    let g = a[i][j];
    let h = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);
}

fn eigsrt(d: &mut [f64; 4], v: &mut [[f64; 4]; 4], n: usize) {
    // Sort eigenvalues into ascending order and rearrange eigenvectors.
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

fn ls_plane(
    bxyz: &[[f64; 4]],
    n: i64,
    pnormal: &mut [f64; 4],
    ppos: &mut [f64; 4],
    odist: &mut f64,
    adist: &mut [f64],
) -> Result<(), String> {
    if n < 3 {
        return Err("ls_plane: too few atoms".to_string());
    }
    let mut cov = [[0.0f64; 4]; 4];
    cov_matrix(bxyz, n, &mut cov);
    let mut v = [[0.0f64; 4]; 4];
    let mut d = [0.0f64; 4];
    jacobi(&mut cov, 3, &mut d, &mut v);
    eigsrt(&mut d, &mut v, 3);

    let mut nml = false;
    for i in 1..=3 {
        for j in 1..=3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            if (v[i][j] - expected).abs() > XEPS {
                nml = true;
                break;
            }
        }
        if nml {
            break;
        }
    }
    if nml {
        // pick the first eigenvector (smallest eigenvalue) as normal.
        for i in 1..=3 {
            pnormal[i] = v[i][1];
        }
    } else {
        // Legacy: if V is identity, force z-axis as normal.
        pnormal[1] = 0.0;
        pnormal[2] = 0.0;
        pnormal[3] = 1.0;
    }
    // average position
    for i in 1..=3 {
        ppos[i] = 0.0;
    }
    for i in 1..=n as usize {
        for j in 1..=3 {
            ppos[j] += bxyz[i][j];
        }
    }
    for j in 1..=3 {
        ppos[j] /= n as f64;
    }

    if pnormal[3] < 0.0 {
        for i in 1..=3 {
            pnormal[i] = -pnormal[i];
        }
    }
    *odist = dot(ppos, pnormal);
    for i in 1..=n as usize {
        adist[i] = dot(&bxyz[i], pnormal) - *odist;
    }
    Ok(())
}
