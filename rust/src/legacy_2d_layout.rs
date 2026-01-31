// Legacy 2D layout (best-pair selection + helix partitioning + XY coordinates).
//
// Ported from legacy C: src/fpair.c, src/fpair_sub.c, src/ps-xy.c, src/ps-xy-sub.c.
//
// This module intentionally uses 1-based indexing (like legacy) to maximize output parity.

use crate::legacy_alg::{compute_base_info, dot, vec_norm, veclen, AtomName4};
use crate::legacy_pairing::check_pairs;
use crate::PairsJson;

const XBIG: f64 = 1.0e18;
const HLXANG: f64 = 0.26;
const MFACTOR: f64 = 10000.0;

#[derive(Debug, Clone)]
pub(crate) struct Layout2d {
    pub(crate) base_xy: Vec<[f64; 3]>, // 1-based: [_, x, y]
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
    seidx: &[[usize; 3]],
    ry: &[i32],
    atom_name: &[AtomName4],
    xyz: &mut [[f64; 4]],
    pairs_all: &PairsJson,
) -> Result<Layout2d, String> {
    if num_residue == 0 {
        return Err("compute_layout_2d: num_residue=0".to_string());
    }

    let base_info = compute_base_info(num_residue, bseq, seidx, ry, atom_name, xyz)?;

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

    let mut matched_idx: Vec<i64> = vec![0; num_residue + 1];
    let mut best_rows: Vec<[i64; 18]> = Vec::new();
    best_rows.push([0; 18]); // 1-based

    for i in 1..=num_residue {
        if let Some(pi) = best_pair(
            i,
            &cand,
            ry,
            &mut matched_idx,
            bseq,
            seidx,
            xyz,
            &base_info.nxyz,
            &base_info.orien,
            &base_info.org,
            atom_name,
            &base_info.bprs,
        )? {
            let j = pi.partner;
            if let Some(pj) = best_pair(
                j,
                &cand,
                ry,
                &mut matched_idx,
                bseq,
                seidx,
                xyz,
                &base_info.nxyz,
                &base_info.orien,
                &base_info.org,
                atom_name,
                &base_info.bprs,
            )? {
                if pj.partner == i {
                    matched_idx[i] = 1;
                    matched_idx[j] = 1;
                    best_rows.push(pi.as_base_pairs_row());
                }
            }
        }
    }

    if best_rows.len() <= 1 {
        return Err("compute_layout_2d: no best pairs".to_string());
    }

    // TODO: port re_ordering + process_2d_fig.
    Err("compute_layout_2d: layout not implemented".to_string())
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
        row[1] = 0; // caller fills
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
    matched_idx: &mut [i64],
    bseq: &[u8],
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
    for &j in cand.get(i).unwrap_or(&Vec::new()) {
        if j == i {
            continue;
        }
        if ry.get(j).copied().unwrap_or(-1) < 0 {
            continue;
        }
        if matched_idx.get(j).copied().unwrap_or(0) != 0 {
            continue;
        }

        let bpid = check_pairs(i, j, bseq, seidx, xyz, nxyz, orien, org, atom_name, bprs, &mut rtn_val, 0);
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

    if let Some(mut out) = best.clone() {
        // Fill the leading residue i into the output row by encoding it at rounded[-]? Caller will set.
        // We keep i separately by convention: the caller knows i.
        // Ensure the row uses legacy rounding; nothing else to do here.
        // The caller will set base_pairs_row[1] to i.
        out.partner = best.unwrap().partner;
        Ok(Some(out))
    } else {
        Ok(None)
    }
}

