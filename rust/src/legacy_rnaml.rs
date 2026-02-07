// Legacy RNAML (XML) writer.
//
// Ported from legacy C: src/rnaxml-new.c.

use crate::{BasePair, Source};
use crate::legacy_alg::{AtomName4, ResName3};
use std::path::Path;

const TAG1: &str = "   ";
const TAG2: &str = "      ";
const TAG3: &str = "         ";
const TAG4: &str = "            ";
const TAG5: &str = "               ";
const TAG6: &str = "                  ";
const TAG7: &str = "                     ";
const TAG8: &str = "                        ";

macro_rules! w {
    ($out:expr, $($args:tt)*) => {{
        std::fmt::Write::write_fmt($out, format_args!($($args)*))
            .map_err(|_| "formatting error".to_string())
    }};
}

macro_rules! wl {
    ($out:expr $(,)?) => {{
        std::fmt::Write::write_str($out, "\n").map_err(|_| "formatting error".to_string())
    }};
    ($out:expr, $($args:tt)*) => {{
        std::fmt::Write::write_fmt($out, format_args!($($args)*))
            .and_then(|_| std::fmt::Write::write_str($out, "\n"))
            .map_err(|_| "formatting error".to_string())
    }};
}

fn resname_eq(r: ResName3, b: [u8; 3]) -> bool {
    r.0 == b
}

fn is_standard_resname(r: ResName3) -> bool {
    resname_eq(r, *b"  A")
        || resname_eq(r, *b"ADE")
        || resname_eq(r, *b"  G")
        || resname_eq(r, *b"GUA")
        || resname_eq(r, *b"  U")
        || resname_eq(r, *b"URA")
        || resname_eq(r, *b"  C")
        || resname_eq(r, *b"CYT")
        || resname_eq(r, *b"  T")
        || resname_eq(r, *b"THY")
}

fn chain_id_for_residue(res_idx: usize, seidx: &[[usize; 3]], chain_id: &[u8]) -> Result<u8, String> {
    let atom_idx = seidx
        .get(res_idx)
        .and_then(|v| v.get(1).copied())
        .ok_or_else(|| format!("missing seidx[{res_idx}][1]"))?;
    chain_id
        .get(atom_idx)
        .copied()
        .ok_or_else(|| format!("missing chain_id[{atom_idx}]"))
}

fn get_chain_idx(num_residue: usize, seidx: &[[usize; 3]], chain_id: &[u8]) -> Result<Vec<[usize; 3]>, String> {
    if num_residue == 0 {
        return Err("get_chain_idx: num_residue=0".to_string());
    }
    let mut chain_idx: Vec<[usize; 3]> = vec![[0; 3]; num_residue + 1];
    chain_idx[1][1] = 1;
    let mut nchain: usize = 1;
    for j in 2..=num_residue {
        let k1 = seidx
            .get(j - 1)
            .and_then(|v| v.get(1).copied())
            .ok_or_else(|| format!("missing seidx[{}][1]", j - 1))?;
        let k2 = seidx
            .get(j)
            .and_then(|v| v.get(1).copied())
            .ok_or_else(|| format!("missing seidx[{j}][1]"))?;
        chain_idx[nchain][2] = j - 1;
        let c1 = chain_id.get(k1).copied().unwrap_or(b' ');
        let c2 = chain_id.get(k2).copied().unwrap_or(b' ');
        if c1 != c2 {
            nchain += 1;
            chain_idx[nchain][1] = j;
        }
    }
    chain_idx[nchain][2] = num_residue;
    chain_idx.truncate(nchain + 1);
    Ok(chain_idx)
}

fn single_continue(single_base: &[i64]) -> Vec<(i64, i64)> {
    let mut out: Vec<(i64, i64)> = Vec::new();
    if single_base.is_empty() {
        return out;
    }
    let mut j: usize = 0;
    while j < single_base.len() {
        let start = single_base[j];
        let mut end = start;
        j += 1;
        while j < single_base.len() && single_base[j] == end + 1 {
            end = single_base[j];
            j += 1;
        }
        out.push((start, end));
    }
    out
}

fn out_edge_orientation(bp: &BasePair) -> (char, char, char) {
    let bang = bp.note.as_deref().is_some_and(|s| s.contains('!')) || bp.text.as_deref().is_some_and(|s| s.contains('!'));
    if bang {
        return ('!', '!', '!');
    }
    let lw = bp.lw.as_deref().unwrap_or("?");
    let edge_5p = lw.chars().next().unwrap_or('?');
    let edge_3p = lw.chars().nth(2).unwrap_or('?');
    let orientation = bp.orientation.as_deref().unwrap_or("?");
    let bond = orientation.chars().next().unwrap_or('?');
    (edge_5p, edge_3p, bond)
}

fn write_base_pair_mol(
    out: &mut String,
    mol_id: i64,
    chain_res: &[i64],
    position_1_n: &[i64],
    base_pairs_out_order: &[BasePair],
) -> Result<(), String> {
    for bp in base_pairs_out_order {
        if bp.kind != "pair" || bp.i <= 0 || bp.j <= 0 {
            continue;
        }
        let nres1 = bp.i as usize;
        let nres2 = bp.j as usize;
        if chain_res.get(nres1).copied().unwrap_or(0) != mol_id || chain_res.get(nres2).copied().unwrap_or(0) != mol_id {
            continue;
        }

        let (edge_5p, edge_3p, bond) = out_edge_orientation(bp);
        let pos1 = position_1_n.get(nres1).copied().unwrap_or(0);
        let pos2 = position_1_n.get(nres2).copied().unwrap_or(0);

        wl!(out, "{TAG5}<base-pair comment=\"?\">")?;
        wl!(out, "{TAG6}<base-id-5p>")?;
        wl!(out, "{TAG7}<base-id><position>{pos1}</position></base-id>")?;
        wl!(out, "{TAG6}</base-id-5p>")?;
        wl!(out, "{TAG6}<base-id-3p>")?;
        wl!(out, "{TAG7}<base-id><position>{pos2}</position></base-id>")?;
        wl!(out, "{TAG6}</base-id-3p>")?;
        wl!(out, "{TAG6}<edge-5p>{edge_5p}</edge-5p>")?;
        wl!(out, "{TAG6}<edge-3p>{edge_3p}</edge-3p>")?;
        wl!(out, "{TAG6}<bond-orientation>{bond}</bond-orientation>")?;
        wl!(out, "{TAG5}</base-pair>")?;
    }
    Ok(())
}

fn write_base_pair_int(
    out: &mut String,
    chain1: i64,
    chain2: i64,
    chain_res: &[i64],
    position_1_n: &[i64],
    base_pairs_out_order: &[BasePair],
) -> Result<(), String> {
    for bp in base_pairs_out_order {
        if bp.kind != "pair" || bp.i <= 0 || bp.j <= 0 {
            continue;
        }
        let nres1 = bp.i as usize;
        let nres2 = bp.j as usize;
        let r1 = chain_res.get(nres1).copied().unwrap_or(0);
        let r2 = chain_res.get(nres2).copied().unwrap_or(0);
        if !((chain1 == r1 && chain2 == r2) || (chain1 == r2 && chain2 == r1)) {
            continue;
        }

        let (edge_5p, edge_3p, bond) = out_edge_orientation(bp);
        let pos1 = position_1_n.get(nres1).copied().unwrap_or(0);
        let pos2 = position_1_n.get(nres2).copied().unwrap_or(0);

        wl!(out, "{TAG5}<base-pair comment=\"?\">")?;
        wl!(out, "{TAG6}<base-id-5p>")?;
        wl!(out, "{TAG7}<base-id>")?;
        wl!(out, "{TAG8}<molecule-id ref=\"{chain1}\"/><position>{pos1}</position>")?;
        wl!(out, "{TAG7}</base-id>")?;
        wl!(out, "{TAG6}</base-id-5p>")?;

        wl!(out, "{TAG6}<base-id-3p>")?;
        wl!(out, "{TAG7}<base-id>")?;
        wl!(out, "{TAG8}<molecule-id ref=\"{chain2}\"/><position>{pos2}</position>")?;
        wl!(out, "{TAG7}</base-id>")?;
        wl!(out, "{TAG6}</base-id-3p>")?;

        wl!(out, "{TAG6}<edge-5p>{edge_5p}</edge-5p>")?;
        wl!(out, "{TAG6}<edge-3p>{edge_3p}</edge-3p>")?;
        wl!(out, "{TAG6}<bond-orientation>{bond}</bond-orientation>")?;
        wl!(out, "{TAG5}</base-pair>")?;
    }
    Ok(())
}

fn write_helix_mol(
    out: &mut String,
    chain1: i64,
    chain2: i64,
    chain_res: &[i64],
    position_1_n: &[i64],
    layout: &crate::legacy_2d_layout::Layout2d,
) -> Result<(), String> {
    let mut nh: i64 = 0;
    for i in 1..=layout.xml_nh {
        let mut k: i64 = 0;
        let mut helix_left: i64 = 0;
        let mut helix_right: i64 = 0;
        let len = layout.xml_helix_len.get(i).copied().unwrap_or(0);
        for j in 0..len {
            let k1 = layout.xml_helix.get(i).and_then(|row| row.get(1).copied()).unwrap_or(0) + j;
            let k2 = layout.xml_helix.get(i).and_then(|row| row.get(2).copied()).unwrap_or(0) - j;
            if k1 <= 0 || k2 <= 0 {
                continue;
            }
            let k1u = k1 as usize;
            let k2u = k2 as usize;
            let r1 = chain_res.get(k1u).copied().unwrap_or(0);
            let r2 = chain_res.get(k2u).copied().unwrap_or(0);

            if chain1 > 0 && chain2 > 0 {
                if (r1 == chain1 && r2 == chain2) || (r2 == chain1 && r1 == chain2) {
                    if k == 0 {
                        helix_left = k1;
                        helix_right = k2;
                    }
                    k += 1;
                }
            } else if chain1 > 0 && chain2 < 0 {
                if (r1 == chain1 && r2 == chain1) || (r2 == chain1 && r1 == chain1) {
                    if k == 0 {
                        helix_left = k1;
                        helix_right = k2;
                    }
                    k += 1;
                }
            }
        }

        if k <= 0 {
            continue;
        }
        nh += 1;
        wl!(out, "{TAG5}<helix id=\"H{nh}\">")?;

        if chain1 > 0 && chain2 > 0 {
            let pos_left = position_1_n.get(helix_left as usize).copied().unwrap_or(0);
            let pos_right = position_1_n.get(helix_right as usize).copied().unwrap_or(0);

            wl!(out, "{TAG6}<base-id-5p>")?;
            wl!(out, "{TAG7}<base-id>")?;
            wl!(out, "{TAG8}<molecule-id ref=\"{chain1}\"/><position>{pos_left}</position>")?;
            wl!(out, "{TAG7}</base-id>")?;
            wl!(out, "{TAG6}</base-id-5p>")?;

            wl!(out, "{TAG6}<base-id-3p>")?;
            wl!(out, "{TAG7}<base-id>")?;
            wl!(out, "{TAG8}<molecule-id ref=\"{chain2}\"/><position>{pos_right}</position>")?;
            wl!(out, "{TAG7}</base-id>")?;
            wl!(out, "{TAG6}</base-id-3p>")?;
        } else {
            let pos_left = position_1_n.get(helix_left as usize).copied().unwrap_or(0);
            let pos_right = position_1_n.get(helix_right as usize).copied().unwrap_or(0);

            wl!(out, "{TAG6}<base-id-5p>")?;
            wl!(out, "{TAG7}<base-id><position>{pos_left}</position></base-id>")?;
            wl!(out, "{TAG6}</base-id-5p>")?;

            wl!(out, "{TAG6}<base-id-3p>")?;
            wl!(out, "{TAG7}<base-id><position>{pos_right}</position></base-id>")?;
            wl!(out, "{TAG6}</base-id-3p>")?;
        }

        wl!(out, "{TAG6}<length>{k}</length>")?;
        wl!(out, "{TAG5}</helix>")?;
    }
    Ok(())
}

fn write_molecule(
    out: &mut String,
    mol_id: i64,
    chain_start: usize,
    chain_end: usize,
    chain_nam: u8,
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    resname: &[ResName3],
    chain_id: &[u8],
    resseq: &[i32],
    bseq: &[u8],
    xyz: &[[f64; 4]],
    sugar_syn: &[i64],
    position_1_n: &[i64],
    chain_res: &[i64],
    layout: &crate::legacy_2d_layout::Layout2d,
    has_any_modify: bool,
    has_any_loop: bool,
    singles: &[(i64, i64)],
    base_pairs_out_order: &[BasePair],
) -> Result<(), String> {
    wl!(out, "{TAG1}<molecule id=\"{mol_id}\">")?;
    wl!(out, "{TAG2}<sequence>")?;

    wl!(out, "{TAG3}<numbering-system id=\"{mol_id}\" used-in-file=\"false\">")?;
    wl!(out, "{TAG4}<numbering-range>")?;
    let pos_start = position_1_n.get(chain_start).copied().unwrap_or(0);
    let pos_end = position_1_n.get(chain_end).copied().unwrap_or(0);
    wl!(out, "{TAG5}<start>{pos_start}</start>")?;
    wl!(out, "{TAG5}<end>{pos_end}</end>")?;
    wl!(out, "{TAG4}</numbering-range>")?;
    wl!(out, "{TAG3}</numbering-system>")?;

    let res_length = (chain_end - chain_start + 1) as i64;
    wl!(
        out,
        "{TAG3}<numbering-table length=\"{res_length}\" comment=\"sequence number in pdb file\">"
    )?;
    out.push_str(TAG3);
    for i in chain_start..=chain_end {
        let atom_idx = seidx.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0);
        let seq = resseq.get(atom_idx).copied().unwrap_or(0) as i64;
        w!(out, "{seq:4} ")?;
        if (i % 10) == 0 {
            out.push('\n');
            out.push_str(TAG3);
        }
    }
    out.push('\n');
    wl!(out, "{TAG3}</numbering-table>")?;

    wl!(out, "{TAG3}<seq-data>")?;
    out.push_str(TAG4);
    for i in chain_start..=chain_end {
        let ch = bseq.get(i).copied().unwrap_or(b'?') as char;
        out.push(ch);
        if (i % 10) == 0 {
            out.push(' ');
        }
        if (i % 60) == 0 {
            out.push('\n');
            out.push_str(TAG4);
        }
    }
    out.push('\n');
    wl!(out, "{TAG3}</seq-data>")?;

    if has_any_modify || has_any_loop {
        wl!(out, "{TAG3}<seq-annotation comment=\"?\">")?;
    }

    if has_any_modify {
        for i in chain_start..=chain_end {
            let atom_idx = seidx.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0);
            let r = resname.get(atom_idx).copied().unwrap_or(ResName3([b' '; 3]));
            if is_standard_resname(r) {
                continue;
            }
            let pos = position_1_n.get(i).copied().unwrap_or(0);
            let r_txt = std::str::from_utf8(&r.0).unwrap_or("???");
            wl!(out, "{TAG4}<modification>")?;
            wl!(out, "{TAG5}<base-id><position>{pos}</position></base-id>")?;
            wl!(out, "{TAG5}<modified-type>{r_txt}</modified-type>")?;
            wl!(out, "{TAG4}</modification>")?;
        }
    }

    if has_any_loop {
        for i in 1..=layout.num_loop {
            let start = layout.loop_.get(i).and_then(|row| row.get(1).copied()).unwrap_or(0);
            let end = layout.loop_.get(i).and_then(|row| row.get(2).copied()).unwrap_or(0);
            if start <= 0 || end <= 0 {
                continue;
            }
            let loop_chain = chain_id_for_residue(start as usize, seidx, chain_id)?;
            if loop_chain != chain_nam {
                continue;
            }
            let pos_start = position_1_n.get(start as usize).copied().unwrap_or(0);
            let pos_end = position_1_n.get(end as usize).copied().unwrap_or(0);

            wl!(out, "{TAG4}<segment>")?;
            wl!(out, "{TAG5}<seg-name>LOOP{i}</seg-name>")?;
            wl!(
                out,
                "{TAG5}<base-id-5p><base-id><position>{pos_start}</position></base-id></base-id-5p>"
            )?;
            wl!(
                out,
                "{TAG5}<base-id-3p><base-id><position>{pos_end}</position></base-id></base-id-3p>"
            )?;
            wl!(out, "{TAG4}</segment>")?;
        }
    }

    if has_any_modify || has_any_loop {
        wl!(out, "{TAG3}</seq-annotation>")?;
    }
    wl!(out, "{TAG2}</sequence>")?;

    wl!(out, "{TAG2}<structure>")?;
    wl!(out, "{TAG3}<model id=\"?\">")?;
    wl!(out, "{TAG4}<model-info>")?;
    wl!(out, "{TAG5}<method>Crystallography ?</method>")?;
    wl!(out, "{TAG5}<resolution>? Angstroms</resolution>")?;
    wl!(out, "{TAG4}</model-info>")?;

    let atom_p = AtomName4::from_bytes(*b" P  ");
    let atom_o3 = AtomName4::from_bytes(*b" O3'");
    for i in chain_start..=chain_end {
        let pos = position_1_n.get(i).copied().unwrap_or(0);
        let base = bseq.get(i).copied().unwrap_or(b'?') as char;
        wl!(out, "{TAG4}<base>")?;
        wl!(out, "{TAG5}<position>{pos}</position>")?;
        wl!(out, "{TAG5}<base-type>{base}</base-type>")?;

        let ib = seidx.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0);
        let ie = seidx.get(i).and_then(|v| v.get(2).copied()).unwrap_or(0);
        for j in ib..=ie {
            let an = atom_name.get(j).copied().unwrap_or(AtomName4([0; 4]));
            if an != atom_p && an != atom_o3 {
                continue;
            }
            let an_txt = std::str::from_utf8(&an.0).unwrap_or("????");
            let x = xyz.get(j).and_then(|v| v.get(1).copied()).unwrap_or(0.0);
            let y = xyz.get(j).and_then(|v| v.get(2).copied()).unwrap_or(0.0);
            let z = xyz.get(j).and_then(|v| v.get(3).copied()).unwrap_or(0.0);

            wl!(out, "{TAG5}<atom serial=\"{j}\">")?;
            wl!(out, "{TAG6}<atom-type>{an_txt}</atom-type>")?;
            wl!(out, "{TAG7}<coordinates>{x:.3} {y:.3} {z:.3}</coordinates>")?;
            wl!(out, "{TAG5}</atom>")?;
        }

        wl!(out, "{TAG4}</base>")?;
    }

    wl!(out, "{TAG4}<str-annotation>")?;
    for i in chain_start..=chain_end {
        if sugar_syn.get(i).copied().unwrap_or(0) <= 0 {
            continue;
        }
        let pos = position_1_n.get(i).copied().unwrap_or(0);
        wl!(out, "{TAG5}<base-conformation>")?;
        wl!(out, "{TAG6}<base-id><position>{pos}</position></base-id>")?;
        wl!(out, "{TAG6}<glycosyl>syn</glycosyl>")?;
        wl!(out, "{TAG5}</base-conformation>")?;
    }

    write_base_pair_mol(out, mol_id, chain_res, position_1_n, base_pairs_out_order)?;
    write_helix_mol(out, mol_id, -1, chain_res, position_1_n, layout)?;

    for (idx0, (seg_start, seg_end)) in singles.iter().enumerate() {
        let seg_id = (idx0 + 1) as i64;
        let pos_start = position_1_n.get(*seg_start as usize).copied().unwrap_or(0);
        let pos_end = position_1_n.get(*seg_end as usize).copied().unwrap_or(0);
        wl!(out, "{TAG5}<single-strand>")?;
        wl!(out, "{TAG6}<segment>")?;
        wl!(out, "{TAG7}<seg-name>SG{seg_id}</seg-name>")?;
        wl!(
            out,
            "{TAG7}<base-id-5p><base-id><position>{pos_start}</position></base-id></base-id-5p>"
        )?;
        wl!(
            out,
            "{TAG7}<base-id-3p><base-id><position>{pos_end}</position></base-id></base-id-3p>"
        )?;
        wl!(out, "{TAG6}</segment>")?;
        wl!(out, "{TAG5}</single-strand>")?;
    }

    wl!(out, "{TAG4}</str-annotation>")?;
    wl!(
        out,
        "{TAG4}<secondary-structure-display comment=\"x,y coodinates\">"
    )?;
    for i in chain_start..=chain_end {
        let pos = position_1_n.get(i).copied().unwrap_or(0);
        let x = layout.base_xy.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0.0);
        let y = layout.base_xy.get(i).and_then(|v| v.get(2).copied()).unwrap_or(0.0);
        wl!(out, "{TAG5}<ss-base-coord>")?;
        wl!(out, "{TAG6}<base-id><position>{pos}</position></base-id>")?;
        wl!(out, "{TAG6}<coordinates>{x:.3} {y:.3}</coordinates>")?;
        wl!(out, "{TAG5}</ss-base-coord>")?;
    }
    wl!(out, "{TAG4}</secondary-structure-display>")?;
    wl!(out, "{TAG3}</model>")?;
    wl!(out, "{TAG2}</structure>")?;
    wl!(out, "{TAG1}</molecule>\n")?;
    Ok(())
}

pub(crate) fn write_rnaml_xml(
    source: Option<&Source>,
    base_path: &Path,
    num_residue: usize,
    bseq: &[u8],
    seidx: &[[usize; 3]],
    atom_name: &[AtomName4],
    resname: &[ResName3],
    chain_id: &[u8],
    resseq: &[i32],
    xyz: &[[f64; 4]],
    sugar_syn: &[i64],
    layout: &crate::legacy_2d_layout::Layout2d,
    base_pairs_out_order: &[BasePair],
) -> Result<String, String> {
    let _ = (source, base_path);

    if num_residue == 0 {
        return Err("write_rnaml_xml: num_residue=0".to_string());
    }
    let chain_idx = get_chain_idx(num_residue, seidx, chain_id)?;
    let nchain = chain_idx.len().saturating_sub(1) as i64;

    let mut position_1_n: Vec<i64> = vec![0; num_residue + 1];
    let mut chain_res: Vec<i64> = vec![0; num_residue + 1];

    for i in 1..=(nchain as usize) {
        let start = chain_idx.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0);
        let end = chain_idx.get(i).and_then(|v| v.get(2).copied()).unwrap_or(0);
        let mut j: i64 = 0;
        for k in start..=end {
            j += 1;
            position_1_n[k] = j;
            chain_res[k] = i as i64;
        }
    }

    let has_any_modify = (1..=num_residue).any(|i| {
        let atom_idx = seidx.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0);
        let r = resname.get(atom_idx).copied().unwrap_or(ResName3([b' '; 3]));
        !is_standard_resname(r)
    });
    let has_any_loop = layout.num_loop > 0;

    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\"?>\n");
    out.push_str("<!DOCTYPE rnaml SYSTEM \"rnaml.dtd\">\n");
    out.push_str("\n<rnaml version=\"1.0\">\n\n");

    for i in 1..=(nchain as usize) {
        let chain_start = chain_idx.get(i).and_then(|v| v.get(1).copied()).unwrap_or(0);
        let chain_end = chain_idx.get(i).and_then(|v| v.get(2).copied()).unwrap_or(0);
        let chain_nam = chain_id_for_residue(chain_start, seidx, chain_id)?;

        let mut xml_bases_mol: Vec<i64> = Vec::new();
        for j in 1..=layout.xml_ns {
            let b = layout.xml_bases.get(j).copied().unwrap_or(0);
            if b <= 0 {
                continue;
            }
            let bu = b as usize;
            if chain_res.get(bu).copied().unwrap_or(0) == i as i64 {
                xml_bases_mol.push(b);
            }
        }
        let singles = single_continue(&xml_bases_mol);

        write_molecule(
            &mut out,
            i as i64,
            chain_start,
            chain_end,
            chain_nam,
            seidx,
            atom_name,
            resname,
            chain_id,
            resseq,
            bseq,
            xyz,
            sugar_syn,
            &position_1_n,
            &chain_res,
            layout,
            has_any_modify,
            has_any_loop,
            &singles,
            base_pairs_out_order,
        )?;
    }

    out.push('\n');
    wl!(&mut out, "{TAG1}<interactions>")?;
    wl!(&mut out, "{TAG4}<str-annotation>")?;
    for i in 1..=(nchain as usize) {
        for j in (i + 1)..=(nchain as usize) {
            write_base_pair_int(
                &mut out,
                i as i64,
                j as i64,
                &chain_res,
                &position_1_n,
                base_pairs_out_order,
            )?;
            write_helix_mol(
                &mut out,
                i as i64,
                j as i64,
                &chain_res,
                &position_1_n,
                layout,
            )?;
        }
    }
    wl!(&mut out, "{TAG4}</str-annotation>")?;
    wl!(&mut out, "{TAG1}</interactions>")?;
    out.push_str("</rnaml>\n");
    Ok(out)
}
