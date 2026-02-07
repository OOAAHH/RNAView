// Legacy RNAML -> PostScript renderer.
//
// Ported from legacy C: src/xml2ps.c.

use std::fmt::Write as _;
use std::fs;
use std::path::Path;

const XBIG: f64 = 1.0e+18;
const PSPIONT: i64 = 12;

#[derive(Debug, Clone)]
struct ParsedForPs {
    nchain: usize,
    chain_idx: Vec<[i64; 2]>, // [start,end] inclusive; 0-based residue indices
    resname: Vec<char>,       // per-residue base letter
    author_seq: Vec<i64>,     // per-residue ResSeq
    sugar_syn: Vec<i64>,      // -1 or 1 (syn)
    xy: Vec<[f64; 3]>,        // per-residue x,y at [0],[1]
    o3_prime_xyz: Vec<[f64; 3]>,
    p_xyz: Vec<[f64; 3]>,
    pair_type: Vec<[char; 3]>, // edge5, edge3, bond
    npair_idx: Vec<[i64; 2]>,  // residue indices (0-based)
    helix_idx: Vec<[i64; 2]>,  // residue indices (0-based)
    helix_length: Vec<i64>,
    sing_st: Vec<i64>,
    sing_end: Vec<i64>,
}

#[derive(Debug, Default, Clone)]
struct PairBuild {
    idx_5p: Option<i64>,
    idx_3p: Option<i64>,
    edge_5p: Option<char>,
    edge_3p: Option<char>,
    bond: Option<char>,
}

#[derive(Debug, Default, Clone)]
struct HelixBuild {
    idx_5p: Option<i64>,
    idx_3p: Option<i64>,
    length: Option<i64>,
}

#[derive(Debug, Default, Clone)]
struct SingleBuild {
    start: Option<i64>,
    end: Option<i64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    FiveP,
    ThreeP,
}

fn text_between<'a>(s: &'a str, open: &str, close: &str) -> Option<&'a str> {
    let start = s.find(open)? + open.len();
    let rest = &s[start..];
    let end = rest.find(close)? + start;
    Some(&s[start..end])
}

fn parse_molecule_id_line(trimmed: &str) -> Option<usize> {
    let rest = trimmed.strip_prefix("<molecule id=\"")?;
    let (num, _) = rest.split_once("\">")?;
    num.parse::<usize>().ok()
}

fn parse_molecule_ref(trimmed: &str) -> Option<usize> {
    let rest = trimmed.split("molecule-id").nth(1)?;
    let rest = rest.split("ref=\"").nth(1)?;
    let (num, _) = rest.split_once('"')?;
    num.parse::<usize>().ok()
}

fn parse_position(trimmed: &str) -> Option<i64> {
    text_between(trimmed, "<position>", "</position>")?.trim().parse::<i64>().ok()
}

fn parse_edge(trimmed: &str, tag: &str) -> Option<char> {
    let open = format!("<{tag}>");
    let close = format!("</{tag}>");
    let v = text_between(trimmed, &open, &close)?;
    v.chars().next()
}

fn parse_length(trimmed: &str) -> Option<i64> {
    text_between(trimmed, "<length>", "</length>")?.trim().parse::<i64>().ok()
}

fn parse_xy_coords(trimmed: &str) -> Option<(f64, f64)> {
    let inner = text_between(trimmed, "<coordinates>", "</coordinates>")?;
    let mut it = inner.split_whitespace();
    let x = it.next()?.parse::<f64>().ok()?;
    let y = it.next()?.parse::<f64>().ok()?;
    Some((x, y))
}

fn parse_xyz_coords(trimmed: &str) -> Option<(f64, f64, f64)> {
    let inner = text_between(trimmed, "<coordinates>", "</coordinates>")?;
    let mut it = inner.split_whitespace();
    let x = it.next()?.parse::<f64>().ok()?;
    let y = it.next()?.parse::<f64>().ok()?;
    let z = it.next()?.parse::<f64>().ok()?;
    Some((x, y, z))
}

fn global_idx_for_position(mol_offsets: &[usize], mol_id: usize, position_1based: i64) -> Option<i64> {
    if mol_id == 0 {
        return None;
    }
    let off = mol_offsets.get(mol_id).copied()?;
    let pos0 = position_1based - 1;
    if pos0 < 0 {
        return None;
    }
    Some(off as i64 + pos0)
}

fn parse_rnaml(xml: &str) -> Result<ParsedForPs, String> {
    let mut mol_offsets: Vec<usize> = vec![0];
    let mut mol_len: Vec<usize> = vec![0];
    let mut chain_idx: Vec<[i64; 2]> = Vec::new();

    let mut current_mol_id: Option<usize> = None;
    let mut current_mol_offset: usize = 0;

    let mut in_numbering_table = false;
    let mut in_seq_data = false;
    let mut seq_buf: Vec<char> = Vec::new();
    let mut auth_buf: Vec<i64> = Vec::new();
    let mut seq_committed = false;

    let mut in_base_block = false;
    let mut in_atom_block = false;
    let mut current_atom_type: Option<String> = None;
    let mut current_p = [9999.0f64; 3];
    let mut current_o3 = [9999.0f64; 3];

    let mut in_ss_display = false;

    let mut in_base_conformation = false;
    let mut base_conf_position: Option<i64> = None;
    let mut base_conf_is_syn = false;

    let mut current_pair: Option<PairBuild> = None;
    let mut pair_side: Option<Side> = None;

    let mut current_helix: Option<HelixBuild> = None;
    let mut helix_side: Option<Side> = None;

    let mut current_single: Option<SingleBuild> = None;

    let mut resname: Vec<char> = Vec::new();
    let mut author_seq: Vec<i64> = Vec::new();
    let mut sugar_syn: Vec<i64> = Vec::new();
    let mut xy: Vec<[f64; 3]> = Vec::new();
    let mut o3_prime_xyz: Vec<[f64; 3]> = Vec::new();
    let mut p_xyz: Vec<[f64; 3]> = Vec::new();

    let mut pair_type: Vec<[char; 3]> = Vec::new();
    let mut npair_idx: Vec<[i64; 2]> = Vec::new();
    let mut helix_idx: Vec<[i64; 2]> = Vec::new();
    let mut helix_length: Vec<i64> = Vec::new();
    let mut sing_st: Vec<i64> = Vec::new();
    let mut sing_end: Vec<i64> = Vec::new();

    for line in xml.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if let Some(mol_id) = parse_molecule_id_line(trimmed) {
            current_mol_id = Some(mol_id);
            current_mol_offset = resname.len();
            seq_buf.clear();
            auth_buf.clear();
            seq_committed = false;

            if mol_offsets.len() <= mol_id {
                mol_offsets.resize(mol_id + 1, 0);
                mol_len.resize(mol_id + 1, 0);
            }
            mol_offsets[mol_id] = current_mol_offset;
            continue;
        }

        if trimmed == "</molecule>" {
            current_mol_id = None;
            continue;
        }

        if trimmed.starts_with("<numbering-table") {
            in_numbering_table = true;
            auth_buf.clear();
            continue;
        }
        if in_numbering_table {
            if trimmed == "</numbering-table>" {
                in_numbering_table = false;
                continue;
            }
            for tok in trimmed.split_whitespace() {
                if let Ok(n) = tok.parse::<i64>() {
                    auth_buf.push(n);
                }
            }
            continue;
        }

        if trimmed == "<seq-data>" {
            in_seq_data = true;
            seq_buf.clear();
            continue;
        }
        if in_seq_data {
            if trimmed == "</seq-data>" {
                in_seq_data = false;
                let mol_id = current_mol_id.ok_or_else(|| "seq-data outside molecule".to_string())?;
                let len = seq_buf.len();
                let off = current_mol_offset;
                if mol_len.len() <= mol_id {
                    mol_len.resize(mol_id + 1, 0);
                }
                mol_len[mol_id] = len;
                if !chain_idx.iter().any(|c| c[0] == off as i64) {
                    let end = off + len;
                    chain_idx.push([off as i64, (end as i64) - 1]);
                }
                for ch in &seq_buf {
                    resname.push(*ch);
                }
                for n in &auth_buf {
                    author_seq.push(*n);
                }
                sugar_syn.extend(std::iter::repeat(-1).take(len));
                seq_committed = true;
                continue;
            }
            for ch in trimmed.chars() {
                if ch.is_whitespace() {
                    continue;
                }
                seq_buf.push(ch);
            }
            continue;
        }

        if trimmed == "<base>" {
            in_base_block = true;
            current_p = [9999.0; 3];
            current_o3 = [9999.0; 3];
            continue;
        }
        if trimmed == "</base>" {
            in_base_block = false;
            o3_prime_xyz.push(current_o3);
            p_xyz.push(current_p);
            continue;
        }

        if in_base_block && trimmed.starts_with("<atom ") {
            in_atom_block = true;
            current_atom_type = None;
            continue;
        }
        if in_atom_block && trimmed == "</atom>" {
            in_atom_block = false;
            current_atom_type = None;
            continue;
        }
        if in_atom_block && trimmed.starts_with("<atom-type>") {
            if let Some(v) = text_between(trimmed, "<atom-type>", "</atom-type>") {
                current_atom_type = Some(v.to_string());
            }
            continue;
        }
        if in_atom_block && trimmed.starts_with("<coordinates>") {
            let Some((x, y, z)) = parse_xyz_coords(trimmed) else {
                continue;
            };
            let atom = current_atom_type.as_deref().unwrap_or("").trim();
            if atom == "P" {
                current_p = [x, y, z];
            } else if atom == "O3'" {
                current_o3 = [x, y, z];
            }
            continue;
        }

        if trimmed.starts_with("<secondary-structure-display") {
            in_ss_display = true;
            continue;
        }
        if trimmed == "</secondary-structure-display>" {
            in_ss_display = false;
            continue;
        }
        if in_ss_display && trimmed.starts_with("<coordinates>") {
            if let Some((x, y)) = parse_xy_coords(trimmed) {
                xy.push([x, y, 0.0]);
            }
            continue;
        }

        if trimmed == "<base-conformation>" {
            in_base_conformation = true;
            base_conf_position = None;
            base_conf_is_syn = false;
            continue;
        }
        if in_base_conformation {
            if trimmed.contains("<position>") {
                base_conf_position = parse_position(trimmed);
                continue;
            }
            if trimmed.starts_with("<glycosyl>") {
                if let Some(v) = text_between(trimmed, "<glycosyl>", "</glycosyl>") {
                    if v.trim().eq_ignore_ascii_case("syn") {
                        base_conf_is_syn = true;
                    }
                }
                continue;
            }
            if trimmed == "</base-conformation>" {
                let mol_id = current_mol_id.ok_or_else(|| "base-conformation outside molecule".to_string())?;
                if !seq_committed {
                    return Err("base-conformation before seq-data committed".to_string());
                }
                if base_conf_is_syn {
                    if let Some(pos) = base_conf_position {
                        if let Some(idx) = global_idx_for_position(&mol_offsets, mol_id, pos) {
                            let iu = idx as usize;
                            if iu < sugar_syn.len() {
                                sugar_syn[iu] = 1;
                            }
                        }
                    }
                }
                in_base_conformation = false;
                base_conf_position = None;
                base_conf_is_syn = false;
                continue;
            }
        }

        if trimmed.starts_with("<base-pair ") {
            current_pair = Some(PairBuild::default());
            pair_side = None;
            continue;
        }
        if trimmed == "</base-pair>" {
            if let Some(p) = current_pair.take() {
                let idx1 = p.idx_5p.ok_or_else(|| "base-pair missing 5p index".to_string())?;
                let idx2 = p.idx_3p.ok_or_else(|| "base-pair missing 3p index".to_string())?;
                let e5 = p.edge_5p.unwrap_or('!');
                let e3 = p.edge_3p.unwrap_or('!');
                let bo = p.bond.unwrap_or('!');
                pair_type.push([e5, e3, bo]);
                npair_idx.push([idx1, idx2]);
            }
            pair_side = None;
            continue;
        }
        if current_pair.is_some() {
            if trimmed == "<base-id-5p>" {
                pair_side = Some(Side::FiveP);
                continue;
            }
            if trimmed == "</base-id-5p>" {
                pair_side = None;
                continue;
            }
            if trimmed == "<base-id-3p>" {
                pair_side = Some(Side::ThreeP);
                continue;
            }
            if trimmed == "</base-id-3p>" {
                pair_side = None;
                continue;
            }
            if trimmed.contains("<position>") {
                let pos = parse_position(trimmed);
                if let (Some(side), Some(pos)) = (pair_side, pos) {
                    let mol_id = parse_molecule_ref(trimmed).or(current_mol_id);
                    if let Some(mol_id) = mol_id {
                        if let Some(idx) = global_idx_for_position(&mol_offsets, mol_id, pos) {
                            if let Some(pair) = current_pair.as_mut() {
                                match side {
                                    Side::FiveP => pair.idx_5p = Some(idx),
                                    Side::ThreeP => pair.idx_3p = Some(idx),
                                }
                            }
                        }
                    }
                }
                continue;
            }
            if let Some(ch) = parse_edge(trimmed, "edge-5p") {
                if let Some(pair) = current_pair.as_mut() {
                    pair.edge_5p = Some(ch);
                }
                continue;
            }
            if let Some(ch) = parse_edge(trimmed, "edge-3p") {
                if let Some(pair) = current_pair.as_mut() {
                    pair.edge_3p = Some(ch);
                }
                continue;
            }
            if let Some(ch) = parse_edge(trimmed, "bond-orientation") {
                if let Some(pair) = current_pair.as_mut() {
                    pair.bond = Some(ch);
                }
                continue;
            }
        }

        if trimmed.starts_with("<helix ") {
            current_helix = Some(HelixBuild::default());
            helix_side = None;
            continue;
        }
        if trimmed == "</helix>" {
            if let Some(h) = current_helix.take() {
                let idx1 = h.idx_5p.ok_or_else(|| "helix missing 5p index".to_string())?;
                let idx2 = h.idx_3p.ok_or_else(|| "helix missing 3p index".to_string())?;
                let len = h.length.ok_or_else(|| "helix missing length".to_string())?;
                helix_idx.push([idx1, idx2]);
                helix_length.push(len);
            }
            helix_side = None;
            continue;
        }
        if current_helix.is_some() {
            if trimmed == "<base-id-5p>" {
                helix_side = Some(Side::FiveP);
                continue;
            }
            if trimmed == "</base-id-5p>" {
                helix_side = None;
                continue;
            }
            if trimmed == "<base-id-3p>" {
                helix_side = Some(Side::ThreeP);
                continue;
            }
            if trimmed == "</base-id-3p>" {
                helix_side = None;
                continue;
            }
            if trimmed.contains("<position>") {
                let pos = parse_position(trimmed);
                if let (Some(side), Some(pos)) = (helix_side, pos) {
                    let mol_id = parse_molecule_ref(trimmed).or(current_mol_id);
                    if let Some(mol_id) = mol_id {
                        if let Some(idx) = global_idx_for_position(&mol_offsets, mol_id, pos) {
                            if let Some(helix) = current_helix.as_mut() {
                                match side {
                                    Side::FiveP => helix.idx_5p = Some(idx),
                                    Side::ThreeP => helix.idx_3p = Some(idx),
                                }
                            }
                        }
                    }
                }
                continue;
            }
            if let Some(len) = parse_length(trimmed) {
                if let Some(helix) = current_helix.as_mut() {
                    helix.length = Some(len);
                }
                continue;
            }
        }

        if trimmed == "<single-strand>" {
            current_single = Some(SingleBuild::default());
            continue;
        }
        if trimmed == "</single-strand>" {
            if let Some(s) = current_single.take() {
                if let (Some(st), Some(en)) = (s.start, s.end) {
                    sing_st.push(st);
                    sing_end.push(en);
                }
            }
            continue;
        }
        if current_single.is_some() {
            if trimmed.starts_with("<base-id-5p>") && trimmed.contains("<position>") {
                let mol_id = current_mol_id.ok_or_else(|| "single-strand base-id outside molecule".to_string())?;
                let pos = parse_position(trimmed).unwrap_or(0);
                if let Some(idx) = global_idx_for_position(&mol_offsets, mol_id, pos) {
                    if let Some(s) = current_single.as_mut() {
                        s.start = Some(idx);
                    }
                }
                continue;
            }
            if trimmed.starts_with("<base-id-3p>") && trimmed.contains("<position>") {
                let mol_id = current_mol_id.ok_or_else(|| "single-strand base-id outside molecule".to_string())?;
                let pos = parse_position(trimmed).unwrap_or(0);
                if let Some(idx) = global_idx_for_position(&mol_offsets, mol_id, pos) {
                    if let Some(s) = current_single.as_mut() {
                        s.end = Some(idx);
                    }
                }
                continue;
            }
        }
    }

    let nres = resname.len();
    if author_seq.len() != nres {
        return Err(format!("xml2ps parse: author_seq len {} != nres {}", author_seq.len(), nres));
    }
    if sugar_syn.len() != nres {
        return Err(format!("xml2ps parse: sugar_syn len {} != nres {}", sugar_syn.len(), nres));
    }
    if xy.len() != nres {
        return Err(format!("xml2ps parse: xy len {} != nres {}", xy.len(), nres));
    }
    if o3_prime_xyz.len() != nres || p_xyz.len() != nres {
        return Err(format!(
            "xml2ps parse: O3/P len {}/{} != nres {}",
            o3_prime_xyz.len(),
            p_xyz.len(),
            nres
        ));
    }

    let nchain = chain_idx.len();

    Ok(ParsedForPs {
        nchain,
        chain_idx,
        resname,
        author_seq,
        sugar_syn,
        xy,
        o3_prime_xyz,
        p_xyz,
        pair_type,
        npair_idx,
        helix_idx,
        helix_length,
        sing_st,
        sing_end,
    })
}

fn dist(a: [f64; 3], b: [f64; 3], n: i64) -> f64 {
    if n == 2 {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    } else if n == 3 {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        (dx * dx + dy * dy).sqrt()
    } else {
        0.0
    }
}

fn get_chain_broken(o3_prime_xyz: &[[f64; 3]], p_xyz: &[[f64; 3]]) -> Vec<i64> {
    let nres = o3_prime_xyz.len().min(p_xyz.len());
    let mut broken: Vec<i64> = vec![0; nres];
    let mut j: i64 = 0;
    for i in 1..nres {
        let dis = dist(o3_prime_xyz[i - 1], p_xyz[i], 3);
        if dis > 2.2 || i == nres - 1 {
            broken[i] = 1;
            j += 1;
        }
        if j > (nres as i64) - 4 {
            // legacy prints a warning; ignore
        }
    }
    broken
}

fn xml_round(d: f64) -> i64 {
    if d > 0.0 {
        (d + 0.5) as i64
    } else {
        (d - 0.5) as i64
    }
}

fn xml_xy4ps(
    xy: &mut [[f64; 3]],
    ps_size: f64,
    n: i64,
    ps: &mut String,
    basepars_ps_image: &Path,
) -> Result<(), String> {
    let num = xy.len();
    if num == 0 {
        return Ok(());
    }

    let mut max_xy = [0.0f64; 3];
    let mut min_xy = [0.0f64; 3];
    max_xy[0] = -XBIG;
    max_xy[1] = -XBIG;
    min_xy[0] = XBIG;
    min_xy[1] = XBIG;
    for p in xy.iter() {
        if p[0] > max_xy[0] {
            max_xy[0] = p[0];
        }
        if p[1] > max_xy[1] {
            max_xy[1] = p[1];
        }
        if p[0] < min_xy[0] {
            min_xy[0] = p[0];
        }
        if p[1] < min_xy[1] {
            min_xy[1] = p[1];
        }
    }

    for p in xy.iter_mut() {
        p[0] -= min_xy[0];
        p[1] -= min_xy[1];
    }

    let temp = (max_xy[0] - min_xy[0]).max(max_xy[1] - min_xy[1]);
    let scale_factor = if temp.abs() < 0.000001 { 1.0 } else { ps_size / temp };
    for p in xy.iter_mut() {
        p[0] *= scale_factor;
        p[1] *= scale_factor;
    }

    if n <= 1 {
        return Ok(());
    }

    // recompute max after scaling
    max_xy[0] = -XBIG;
    max_xy[1] = -XBIG;
    for p in xy.iter() {
        if p[0] > max_xy[0] {
            max_xy[0] = p[0];
        }
        if p[1] > max_xy[1] {
            max_xy[1] = p[1];
        }
    }

    let mut urxy = [0i64; 3];
    urxy[0] = xml_round(max_xy[0]);
    urxy[1] = xml_round(max_xy[1]);

    let paper_size = [8.5f64, 11.0f64];
    let mut llxy = [0i64; 3];
    llxy[0] = xml_round(0.5 * (paper_size[0] * 72.0 - (urxy[0] as f64)));
    llxy[1] = xml_round(0.5 * (paper_size[1] * 72.0 - (urxy[1] as f64)));

    let boundary_offset: i64 = 20;
    let mut bbox = [0i64; 4];
    bbox[0] = llxy[0] - boundary_offset;
    bbox[1] = llxy[1] - boundary_offset;
    bbox[2] = urxy[0] + llxy[0] + boundary_offset;
    bbox[3] = urxy[1] + llxy[1] + boundary_offset;

    ps_head(ps, &bbox, basepars_ps_image)?;
    write!(ps, "{:6}{:6} translate\n\n", llxy[0], llxy[1]).map_err(|e| e.to_string())?;
    Ok(())
}

fn ps_head(ps: &mut String, bbox: &[i64; 4], basepars_ps_image: &Path) -> Result<(), String> {
    ps.push_str("%!PS-Adobe-2.0\n");
    ps.push_str("%%Title: A postscript generated from a RNAML file.\n");
    ps.push_str("%%Orientation: Portrait\n");
    ps.push_str("%%BoundingBox: ");
    for v in bbox {
        write!(ps, "{v:6}").map_err(|e| e.to_string())?;
    }
    ps.push_str("\n\n");
    write!(ps, "/Times-Bold findfont {PSPIONT}  scalefont setfont\n\n").map_err(|e| e.to_string())?;

    let par = fs::read_to_string(basepars_ps_image).map_err(|e| e.to_string())?;
    ps.push_str(&par);
    Ok(())
}

fn label_ps_resname(ps: &mut String, resname: &[char], xy: &[[f64; 3]], sugar_syn: &[i64]) -> Result<(), String> {
    ps.push_str("/center{dup stringwidth pop 2 div neg 0 rmoveto } bind def\n");
    for j in 0..resname.len() {
        let cc = resname[j];
        let x = xy.get(j).map(|v| v[0]).unwrap_or(0.0);
        let y = xy.get(j).map(|v| v[1]).unwrap_or(0.0) - 4.0;
        if sugar_syn.get(j).copied().unwrap_or(-1) <= 0 {
            write!(ps, "{x:.2} {y:.2}  moveto ({cc}) center show \n").map_err(|e| e.to_string())?;
        } else {
            ps.push_str("/center{dup stringwidth pop 2 div neg 0 rmoveto } bind def\n");
            write!(ps, "gsave Al {x:.2} {y:.2}  moveto ({cc}) center show grestore\n").map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

fn slope(k1: i64, k2: i64, xy: &[[f64; 3]]) -> f64 {
    let u1 = k1.max(0) as usize;
    let u2 = k2.max(0) as usize;
    let x1 = xy.get(u1).map(|v| v[0]).unwrap_or(0.0);
    let y1 = xy.get(u1).map(|v| v[1]).unwrap_or(0.0);
    let x2 = xy.get(u2).map(|v| v[0]).unwrap_or(0.0);
    let y2 = xy.get(u2).map(|v| v[1]).unwrap_or(0.0);

    let mut denom = x1 - x2;
    if denom.abs() < 0.00001 {
        denom = 0.00001;
    }
    let mut a = (y1 - y2) / denom;
    if a.abs() < 0.00001 {
        a = 0.00001;
    }
    a
}

fn twopoints(xy0: [f64; 3], a: f64, d: f64) -> ([f64; 3], [f64; 3]) {
    let denom = (1.0 + a * a).sqrt();
    let mut xy1 = [0.0f64; 3];
    let mut xy2 = [0.0f64; 3];
    xy1[1] = xy0[1] + d / denom;
    xy1[2] = xy0[2] + a * d / denom;
    xy2[1] = xy0[1] - d / denom;
    xy2[2] = xy0[2] - a * d / denom;
    (xy1, xy2)
}

fn write_5p_3p(ps: &mut String, k1: i64, k2: i64, a: f64, xy: &[[f64; 3]], label: &str) -> Result<(), String> {
    ps.push_str("/Times findfont 11 scalefont setfont\n");
    ps.push_str("/center{dup stringwidth pop 2 div neg 0 rmoveto } bind def\n");

    let u1 = k1.max(0) as usize;
    let u2 = k2.max(0) as usize;
    let x1 = xy.get(u1).map(|v| v[0]).unwrap_or(0.0);
    let y1 = xy.get(u1).map(|v| v[1]).unwrap_or(0.0);
    let x2 = xy.get(u2).map(|v| v[0]).unwrap_or(0.0);
    let y2 = xy.get(u2).map(|v| v[1]).unwrap_or(0.0);

    let xy0 = [0.0, x1, y1];
    let (p1, p2) = twopoints(xy0, a, 10.0);
    let vec1 = [0.0, x2 - x1, y2 - y1];
    let vec2 = [0.0, p1[1] - x1, p1[2] - y1];
    let sign = vec1[1] * vec2[1] + vec1[2] * vec2[2];

    if sign <= 0.0 {
        write!(ps, "gsave Al {:.2} {:.2} moveto ({label}) center show grestore\n", p1[1], p1[2] - 4.0)
            .map_err(|e| e.to_string())?;
    } else {
        write!(ps, "gsave Al {:.2} {:.2} moveto ({label}) center show grestore\n", p2[1], p2[2] - 4.0)
            .map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn label_5p_3p(ps: &mut String, chain_i: usize, chain_idx: &[[i64; 2]], xy: &[[f64; 3]]) -> Result<(), String> {
    let three_p = "3'";
    let five_p = "5'";

    let start = chain_idx.get(chain_i).map(|v| v[0]).unwrap_or(0);
    let end = chain_idx.get(chain_i).map(|v| v[1]).unwrap_or(0);

    let a5 = slope(start + 1, start, xy);
    write_5p_3p(ps, start, start + 1, a5, xy, five_p)?;

    let a3 = slope(end, end - 1, xy);
    write_5p_3p(ps, end, end - 1, a3, xy, three_p)?;
    Ok(())
}

fn color_line(ps: &mut String, color: &str, xy1: [f64; 3], xy2: [f64; 3]) -> Result<(), String> {
    write!(
        ps,
        "{color} {:.2} {:.2} {:.2} {:.2} LINE\n ",
        xy1[0], xy1[1], xy2[0], xy2[1]
    )
    .map_err(|e| e.to_string())?;
    Ok(())
}

fn link_chain(ps: &mut String, nchain: usize, chain_idx: &[[i64; 2]], xy: &[[f64; 3]], chain_broken: &[i64]) -> Result<(), String> {
    const CHAIN_COLOR: [&str; 20] = [
        "Il", "Tl", "Ul", "Gl", "Am", "Cm", "Gm", "Im", "Xl", "UM", "Um", "AM", "CM", "GM", "IM", "TM", "Al", "Cl",
        "Xl", "Xl",
    ];

    ps.push_str("0.005 setlinewidth 1 setlinejoin 1 setlinecap\n");
    let last = chain_idx
        .get(nchain.saturating_sub(1))
        .map(|v| v[1] - 1)
        .unwrap_or(-1);

    for i in 0..nchain {
        let start = chain_idx.get(i).map(|v| v[0]).unwrap_or(0);
        let end = chain_idx.get(i).map(|v| v[1]).unwrap_or(0);
        if (end - start) < 1 {
            continue;
        }
        label_5p_3p(ps, i, chain_idx, xy)?;
        let color = if i < 19 { CHAIN_COLOR[i] } else { CHAIN_COLOR[19] };

        for j in start..end {
            if j < last {
                let bj = (j + 1) as usize;
                if chain_broken.get(bj).copied().unwrap_or(0) == 1 {
                    continue;
                }
            }
            let u1 = j.max(0) as usize;
            let u2 = (j + 1).max(0) as usize;
            let p1 = xy.get(u1).copied().unwrap_or([0.0; 3]);
            let p2 = xy.get(u2).copied().unwrap_or([0.0; 3]);
            color_line(ps, color, p1, p2)?;
        }
    }
    Ok(())
}

fn label_seq(ps: &mut String, k1: i64, k2: i64, a: f64, xy: &[[f64; 3]], author_seq: &[i64]) -> Result<(), String> {
    let u1 = k1.max(0) as usize;
    let u2 = k2.max(0) as usize;
    let x1 = xy.get(u1).map(|v| v[0]).unwrap_or(0.0);
    let y1 = xy.get(u1).map(|v| v[1]).unwrap_or(0.0);
    let x2 = xy.get(u2).map(|v| v[0]).unwrap_or(0.0);
    let y2 = xy.get(u2).map(|v| v[1]).unwrap_or(0.0);

    let xy0 = [0.0, x1, y1];
    let (mut xy1, mut xy2p) = twopoints(xy0, a, 12.0);

    let vec1 = [0.0, x1 - x2, y1 - y2];
    let vec2 = [0.0, xy1[1] - x1, xy1[2] - y1];
    let sign = vec1[1] * vec2[1] + vec1[2] * vec2[2];

    // Match legacy C logic exactly (including the k1==99 and k1==999 gaps where no label is emitted).
    let pick_d = if k1 < 10 {
        9.0
    } else if k1 >= 10 && k1 < 99 {
        11.0
    } else if k1 >= 100 && k1 < 999 {
        13.0
    } else if k1 > 999 {
        15.0
    } else {
        return Ok(());
    };
    let (p1, p2) = twopoints(xy0, a, pick_d);
    xy1 = p1;
    xy2p = p2;

    let val = author_seq.get(u1).copied().unwrap_or(0);
    if sign >= 0.0 {
        write!(
            ps,
            "gsave Tl {:.2} {:.2} moveto ({val}) center show grestore\n",
            xy1[1],
            xy1[2] - 4.0
        )
        .map_err(|e| e.to_string())?;
    } else {
        write!(
            ps,
            "gsave Tl {:.2} {:.2} moveto ({val}) center show grestore\n",
            xy2p[1],
            xy2p[2] - 4.0
        )
        .map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn label_seq_number(
    ps: &mut String,
    nres: i64,
    helix_idx: &[[i64; 2]],
    helix_length: &[i64],
    sing_st: &[i64],
    sing_end: &[i64],
    xy: &[[f64; 3]],
    author_seq: &[i64],
) -> Result<(), String> {
    ps.push_str("/Times findfont 8 scalefont setfont\n");
    ps.push_str("/center{dup stringwidth pop 2 div neg 0 rmoveto } bind def\n");

    for i in 0..helix_idx.len() {
        let h0 = helix_idx[i][0];
        let h1 = helix_idx[i][1];
        let len = helix_length.get(i).copied().unwrap_or(0);
        for j in 0..len {
            let k1 = h0 + j;
            let k2 = h1 - j;
            let a = slope(k1, k2, xy);
            if ((k1 + 1) % 5) == 0 {
                label_seq(ps, k1, k2, a, xy, author_seq)?;
            }
            if ((k2 + 1) % 5) == 0 {
                label_seq(ps, k2, k1, a, xy, author_seq)?;
            }
        }
    }

    for i in 0..sing_st.len().min(sing_end.len()) {
        let st = sing_st[i];
        let en = sing_end[i];
        for j in st..=en {
            if j == 0 || j == (nres - 1) {
                continue;
            }
            if ((j + 1) % 5) != 0 {
                continue;
            }

            let a = slope(j - 1, j + 1, xy);
            let at = -1.0 / a;

            let uj = j as usize;
            let xj = xy.get(uj).map(|v| v[0]).unwrap_or(0.0);
            let yj = xy.get(uj).map(|v| v[1]).unwrap_or(0.0);
            let xy0 = [0.0, xj, yj];
            let (xy1, xy2) = twopoints(xy0, at, 12.0);

            let xmid = 0.5 * (xy[(j - 1) as usize][0] + xy[(j + 1) as usize][0]);
            let ymid = 0.5 * (xy[(j - 1) as usize][1] + xy[(j + 1) as usize][1]);

            let vec1 = [0.0, xj - xmid, yj - ymid];
            let vec2 = [0.0, xy1[1] - xmid, xy1[2] - ymid];
            let sign = vec1[1] * vec2[1] + vec1[2] * vec2[2];

            let val = author_seq.get(uj).copied().unwrap_or(0);
            if sign >= 0.0 {
                write!(
                    ps,
                    "gsave Tl {:.2} {:.2} moveto ({val}) center show grestore\n",
                    xy1[1],
                    xy1[2] - 4.0
                )
                .map_err(|e| e.to_string())?;
            } else {
                write!(
                    ps,
                    "gsave Tl {:.2} {:.2} moveto ({val}) center show grestore\n",
                    xy2[1],
                    xy2[2] - 4.0
                )
                .map_err(|e| e.to_string())?;
            }
        }
    }
    Ok(())
}

fn line(ps: &mut String, xy1: [f64; 3], xy2: [f64; 3]) -> Result<(), String> {
    write!(ps, "Xl {:.2} {:.2} {:.2} {:.2} LINE\n ", xy1[1], xy1[2], xy2[1], xy2[2]).map_err(|e| e.to_string())?;
    Ok(())
}

fn dashline(ps: &mut String, xy1: [f64; 3], xy2: [f64; 3]) -> Result<(), String> {
    write!(
        ps,
        "Xl {:.2} {:.2} {:.2} {:.2} gsave  DASHLINE grestore\n ",
        xy1[1],
        xy1[2],
        xy2[1],
        xy2[2]
    )
    .map_err(|e| e.to_string())?;
    Ok(())
}

fn dotline(ps: &mut String, xy1: [f64; 3], xy2: [f64; 3]) -> Result<(), String> {
    write!(ps, "Xl {:.2} {:.2} {:.2} {:.2} DOTLINE\n ", xy1[1], xy1[2], xy2[1], xy2[2]).map_err(|e| e.to_string())?;
    Ok(())
}

fn square(ps: &mut String, fill: i64, xy1: [f64; 3], xy2: [f64; 3], d: f64, a: f64) -> Result<(), String> {
    let mut xy0 = [0.0f64; 3];
    xy0[1] = xy1[1];
    xy0[2] = xy1[2];
    let (sxy1, sxy4) = twopoints(xy0, a, d);
    xy0[1] = xy2[1];
    xy0[2] = xy2[2];
    let (sxy2, sxy3) = twopoints(xy0, a, d);

    write!(
        ps,
        "Xl {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} SQUARE\n ",
        sxy1[1],
        sxy1[2],
        sxy2[1],
        sxy2[2],
        sxy3[1],
        sxy3[2],
        sxy4[1],
        sxy4[2]
    )
    .map_err(|e| e.to_string())?;
    if fill == 0 {
        ps.push_str("gsave  grestore stroke\n ");
    } else {
        ps.push_str("gsave FILLBLACK grestore stroke\n ");
    }
    Ok(())
}

fn circle(ps: &mut String, fill: i64, xy1: [f64; 3], xy2: [f64; 3], r: f64) -> Result<(), String> {
    let x0 = 0.5 * (xy1[1] + xy2[1]);
    let y0 = 0.5 * (xy1[2] + xy2[2]);
    write!(ps, "Xl {:.2} {:.2} {:.2} CIRCLE \n ", x0, y0, r).map_err(|e| e.to_string())?;
    if fill == 0 {
        ps.push_str("gsave  grestore stroke\n ");
    } else {
        ps.push_str("gsave FILLBLACK grestore stroke\n ");
    }
    Ok(())
}

fn triangle(ps: &mut String, fill: i64, xy2: [f64; 3], txy3: [f64; 3], d: f64, a: f64) -> Result<(), String> {
    let xy0 = [0.0, xy2[1], xy2[2]];
    let (txy1, txy2p) = twopoints(xy0, a, d);
    write!(
        ps,
        "Xl {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} TRIANGLE\n ",
        txy1[1],
        txy1[2],
        txy2p[1],
        txy2p[2],
        txy3[1],
        txy3[2]
    )
    .map_err(|e| e.to_string())?;
    if fill == 0 {
        ps.push_str("gsave  grestore stroke\n ");
    } else {
        ps.push_str("gsave FILLBLACK grestore stroke\n ");
    }
    Ok(())
}

fn dashline_red(ps: &mut String, xy1: [f64; 3], xy2: [f64; 3]) -> Result<(), String> {
    let dpair = ((xy2[1] - xy1[1]) * (xy2[1] - xy1[1]) + (xy2[2] - xy1[2]) * (xy2[2] - xy1[2])).sqrt();
    let xy0 = [0.0, 0.5 * (xy1[1] + xy2[1]), 0.5 * (xy1[2] + xy2[2])];
    let d = dpair / 2.0 - 6.0;

    let mut denom = xy2[1] - xy1[1];
    if denom.abs() < 0.00001 {
        denom = 0.00001;
    }
    let mut a = (xy2[2] - xy1[2]) / denom;
    if a.abs() < 0.00001 {
        a = 0.00001;
    }
    let (nxy1, nxy2) = twopoints(xy0, a, d);
    write!(
        ps,
        "Al W2 {:.2} {:.2} {:.2} {:.2} gsave  DASHLINE grestore\n ",
        nxy1[1],
        nxy1[2],
        nxy2[1],
        nxy2[2]
    )
    .map_err(|e| e.to_string())?;
    Ok(())
}

fn lw_shapes(
    ps: &mut String,
    bseq: &[char],
    k1: i64,
    k2: i64,
    pair_type: [char; 3],
    x: [f64; 3],
    y: [f64; 3],
    at: f64,
    ratio: f64,
) -> Result<(), String> {
    let dpair = ratio
        * ((x[2] - x[1]) * (x[2] - x[1]) + (y[2] - y[1]) * (y[2] - y[1])).sqrt();
    let xy0 = [0.0, 0.5 * (x[1] + x[2]), 0.5 * (y[1] + y[2])];
    let mut d1 = dpair / 2.0;
    d1 = d1 - 8.0;
    if d1 < 0.0 {
        d1 = 0.0;
    }
    let d2 = 8.0;
    let d3 = 3.0;
    let r = 3.0;

    let a = -1.0 / at;

    let fill = if pair_type[2] == 't' { 0 } else { 1 };
    ps.push_str("NP W1\n ");

    let (mut xy1, mut xy6) = twopoints(xy0, at, d1);
    let dt1 = ((xy1[1] - x[1]) * (xy1[1] - x[1]) + (xy1[2] - y[1]) * (xy1[2] - y[1])).sqrt();
    let dt2 = ((xy6[1] - x[1]) * (xy6[1] - x[1]) + (xy6[2] - y[1]) * (xy6[2] - y[1])).sqrt();

    let (mut xy2, mut xy5) = twopoints(xy0, at, d2);
    let (mut xy3, mut xy4) = twopoints(xy0, at, d3);
    if dt1 >= dt2 {
        std::mem::swap(&mut xy1, &mut xy6);
        std::mem::swap(&mut xy2, &mut xy5);
        std::mem::swap(&mut xy3, &mut xy4);
    }

    if pair_type[0] == '-' && pair_type[1] == '-' {
        line(ps, xy1, xy6)?;
    } else if pair_type[0] == '+' && pair_type[1] == '+' {
        let d = 2.0;
        let (txy1, txy2) = twopoints(xy1, a, d);
        let (txy3, txy4) = twopoints(xy6, a, d);
        line(ps, txy1, txy3)?;
        line(ps, txy2, txy4)?;
    } else if pair_type[0] == 'W' && pair_type[1] == 'W' {
        let c1 = bseq.get(k1.max(0) as usize).copied().unwrap_or('?').to_ascii_uppercase();
        let c2 = bseq.get(k2.max(0) as usize).copied().unwrap_or('?').to_ascii_uppercase();
        if ((c1 == 'G' && c2 == 'U') || (c2 == 'G' && c1 == 'U')) && pair_type[2] == 'c' {
            circle(ps, 0, xy3, xy4, r)?;
        } else {
            line(ps, xy1, xy3)?;
            circle(ps, fill, xy3, xy4, r)?;
            line(ps, xy4, xy6)?;
        }
    } else if pair_type[0] == 'H' && pair_type[1] == 'H' {
        line(ps, xy1, xy3)?;
        square(ps, fill, xy3, xy4, d3, a)?;
        line(ps, xy4, xy6)?;
    } else if pair_type[0] == 'S' && pair_type[1] == 'S' {
        line(ps, xy1, xy3)?;
        triangle(ps, fill, xy3, xy4, d3, a)?;
        line(ps, xy4, xy6)?;
    } else if pair_type[0] == 'W' && pair_type[1] == 'H' {
        line(ps, xy1, xy2)?;
        circle(ps, fill, xy2, xy3, r)?;
        line(ps, xy3, xy4)?;
        square(ps, fill, xy4, xy5, d3, a)?;
        line(ps, xy5, xy6)?;
    } else if pair_type[0] == 'H' && pair_type[1] == 'W' {
        line(ps, xy1, xy2)?;
        square(ps, fill, xy2, xy3, d3, a)?;
        line(ps, xy3, xy4)?;
        circle(ps, fill, xy4, xy5, r)?;
        line(ps, xy5, xy6)?;
    } else if pair_type[0] == 'W' && pair_type[1] == 'S' {
        line(ps, xy1, xy2)?;
        circle(ps, fill, xy2, xy3, r)?;
        line(ps, xy3, xy4)?;
        triangle(ps, fill, xy4, xy5, d3, a)?;
        line(ps, xy5, xy6)?;
    } else if pair_type[0] == 'S' && pair_type[1] == 'W' {
        line(ps, xy1, xy2)?;
        triangle(ps, fill, xy3, xy2, d3, a)?;
        line(ps, xy3, xy4)?;
        circle(ps, fill, xy4, xy5, r)?;
        line(ps, xy5, xy6)?;
    } else if pair_type[0] == 'H' && pair_type[1] == 'S' {
        line(ps, xy1, xy2)?;
        square(ps, fill, xy2, xy3, d3, a)?;
        line(ps, xy3, xy4)?;
        triangle(ps, fill, xy4, xy5, d3, a)?;
        line(ps, xy5, xy6)?;
    } else if pair_type[0] == 'S' && pair_type[1] == 'H' {
        line(ps, xy1, xy2)?;
        triangle(ps, fill, xy3, xy2, d3, a)?;
        line(ps, xy3, xy4)?;
        square(ps, fill, xy4, xy5, d3, a)?;
        line(ps, xy5, xy6)?;
    } else if pair_type[0] == '.' && pair_type[1] == 'W' {
        circle(ps, fill, xy1, xy2, r * 0.1)?;
        dashline(ps, xy2, xy5)?;
        circle(ps, fill, xy5, xy6, r)?;
    } else if pair_type[0] == 'W' && pair_type[1] == '.' {
        circle(ps, fill, xy1, xy2, r * 0.1)?;
        dashline(ps, xy2, xy5)?;
        circle(ps, fill, xy5, xy6, r)?;
    } else if pair_type[0] == '.' && pair_type[1] == 'H' {
        circle(ps, fill, xy1, xy2, r * 0.1)?;
        dashline(ps, xy2, xy5)?;
        square(ps, fill, xy5, xy6, d3, a)?;
    } else if pair_type[0] == 'H' && pair_type[1] == '.' {
        square(ps, fill, xy1, xy2, d3, a)?;
        line(ps, xy2, xy5)?;
        circle(ps, fill, xy5, xy6, r * 0.1)?;
    } else if pair_type[0] == 'H' && pair_type[1] == '.' {
        circle(ps, fill, xy1, xy2, r * 0.1)?;
        line(ps, xy2, xy5)?;
        triangle(ps, fill, xy5, xy6, d3, a)?;
    } else if pair_type[0] == 'H' && pair_type[1] == '.' {
        triangle(ps, fill, xy2, xy1, d3, a)?;
        line(ps, xy2, xy5)?;
        circle(ps, fill, xy5, xy6, r * 0.1)?;
    } else if pair_type[0] == '.' && pair_type[1] == '.' {
        dashline(ps, xy1, xy6)?;
    } else if pair_type[0] == 'X' || pair_type[1] == 'X' {
        line(ps, xy1, xy6)?;
    }

    Ok(())
}

fn draw_lw_diagram(
    ps: &mut String,
    pair_type: &[[char; 3]],
    bseq: &[char],
    npair_idx: &[[i64; 2]],
    xy: &[[f64; 3]],
) -> Result<(), String> {
    for i in 0..pair_type.len().min(npair_idx.len()) {
        if pair_type[i][0] != '!' {
            let k1 = npair_idx[i][0];
            let k2 = npair_idx[i][1];
            let u1 = k1.max(0) as usize;
            let u2 = k2.max(0) as usize;
            let x12 = [0.0, xy.get(u1).map(|v| v[0]).unwrap_or(0.0), xy.get(u2).map(|v| v[0]).unwrap_or(0.0)];
            let y12 = [0.0, xy.get(u1).map(|v| v[1]).unwrap_or(0.0), xy.get(u2).map(|v| v[1]).unwrap_or(0.0)];
            let a = slope(k1, k2, xy);
            lw_shapes(ps, bseq, k1, k2, pair_type[i], x12, y12, a, 1.0)?;
        } else {
            let k1 = npair_idx[i][0].max(0) as usize;
            let k2 = npair_idx[i][1].max(0) as usize;
            let mut xy1 = [0.0f64; 3];
            let mut xy2p = [0.0f64; 3];
            xy1[1] = xy.get(k1).map(|v| v[0]).unwrap_or(0.0);
            xy1[2] = xy.get(k1).map(|v| v[1]).unwrap_or(0.0);
            xy2p[1] = xy.get(k2).map(|v| v[0]).unwrap_or(0.0);
            xy2p[2] = xy.get(k2).map(|v| v[1]).unwrap_or(0.0);
            dashline_red(ps, xy1, xy2p)?;
        }
    }
    Ok(())
}

fn h_width(helix_idx: &[[i64; 2]], xy: &[[f64; 3]]) -> f64 {
    if helix_idx.is_empty() {
        return 0.0;
    }
    let k1 = helix_idx[0][0].max(0) as usize;
    let k2 = helix_idx[0][1].max(0) as usize;
    let a = xy.get(k1).copied().unwrap_or([0.0; 3]);
    let b = xy.get(k2).copied().unwrap_or([0.0; 3]);
    dist(a, b, 2)
}

pub(crate) fn ps_from_rnaml_xml(xml: &str, basepars_ps_image: &Path) -> Result<String, String> {
    let parsed = parse_rnaml(xml)?;
    let nres = parsed.resname.len() as i64;

    let mut ps = String::new();
    let mut xy = parsed.xy.clone();
    let mut xy0 = parsed.xy.clone();

    // rotate 180 on x axis: y -> -y (only for the final xy)
    for p in xy.iter_mut() {
        p[1] = -p[1];
    }

    // Rescale xy0 for helix width estimation (no PS output).
    xml_xy4ps(&mut xy0, 550.0, 1, &mut ps, basepars_ps_image)?;

    let helix_width = h_width(&parsed.helix_idx, &xy0);
    let default_size = 550.0;
    let mut ps_width = helix_width;
    if helix_width > 70.0 {
        ps_width = 70.0;
    } else if helix_width < 30.0 {
        ps_width = 30.0;
    }
    let mut ps_size = if helix_width.abs() < 0.000001 {
        default_size
    } else {
        default_size * ps_width / helix_width
    };
    if ps_size > 650.0 {
        ps_size = 650.0;
    }

    // Rescale xy and emit PS header/translate.
    xml_xy4ps(&mut xy, ps_size, 2, &mut ps, basepars_ps_image)?;

    let chain_broken = get_chain_broken(&parsed.o3_prime_xyz, &parsed.p_xyz);

    label_ps_resname(&mut ps, &parsed.resname, &xy, &parsed.sugar_syn)?;
    link_chain(&mut ps, parsed.nchain, &parsed.chain_idx, &xy, &chain_broken)?;
    label_seq_number(
        &mut ps,
        nres,
        &parsed.helix_idx,
        &parsed.helix_length,
        &parsed.sing_st,
        &parsed.sing_end,
        &xy,
        &parsed.author_seq,
    )?;
    draw_lw_diagram(&mut ps, &parsed.pair_type, &parsed.resname, &parsed.npair_idx, &xy)?;
    ps.push_str("\n showpage \n");
    Ok(ps)
}
