use pdbtbx::{
    ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain,
    Format, ReadOptions, StrictnessLevel,
};
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseResidue {
    pub chain: char,
    pub resseq: i32,
    pub insertion_code: Option<char>,
    pub base: char,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct ResidueKey {
    chain_id: String,
    resseq: i32,
    insertion_code: Option<char>,
    resname: String,
}

#[derive(Debug, Clone)]
struct AtomRec {
    name: String,
    x: f64,
    y: f64,
    z: f64,
}

fn canonical_atom_name(raw: &str) -> String {
    let mut s = raw.trim().to_ascii_uppercase();
    if s.contains('*') {
        s = s.replace('*', "'");
    }
    match s.as_str() {
        "O1'" => "O4'".to_string(),
        "OL" => "O1P".to_string(),
        "OR" => "O2P".to_string(),
        "C5A" => "C5M".to_string(),
        "O5T" => "O5'".to_string(),
        "O3T" => "O3'".to_string(),
        _ => s,
    }
}

fn canonical_residue_name(raw: &str) -> String {
    let s = raw.trim().to_ascii_uppercase();
    if s.len() == 2 && s.starts_with('D') {
        let mut chars = s.chars();
        let _ = chars.next();
        if let Some(second) = chars.next() {
            if matches!(second, 'A' | 'T' | 'G' | 'C') {
                return second.to_string();
            }
        }
    }
    s
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn find_atom_pos(atoms: &[AtomRec], target: &str) -> Option<[f64; 3]> {
    atoms
        .iter()
        .find(|a| a.name == target)
        .map(|a| [a.x, a.y, a.z])
}

fn has_atom(atoms: &[AtomRec], target: &str) -> bool {
    atoms.iter().any(|a| a.name == target)
}

fn residue_ident(atoms: &[AtomRec]) -> i32 {
    let n9 = find_atom_pos(atoms, "N9");
    let n1 = find_atom_pos(atoms, "N1");
    let c2 = find_atom_pos(atoms, "C2");
    let c6 = find_atom_pos(atoms, "C6");

    let dcrt = 2.0;
    let dcrt2 = 3.0;

    if let (Some(n1), Some(c2), Some(c6)) = (n1, c2, c6) {
        let d1 = dist(n1, c2);
        let d2 = dist(n1, c6);
        let d3 = dist(c2, c6);
        if d1 <= dcrt && d2 <= dcrt && d3 <= dcrt2 {
            let mut id = 0; // pyrimidine
            if let Some(n9) = n9 {
                let d = dist(n1, n9);
                if (3.5..=4.5).contains(&d) {
                    id = 1; // purine
                }
            }
            return id;
        }
        return -2;
    }

    let ca = find_atom_pos(atoms, "CA");
    let c = find_atom_pos(atoms, "C").or_else(|| find_atom_pos(atoms, "N"));
    if let (Some(ca), Some(c)) = (ca, c) {
        if dist(ca, c) <= dcrt {
            return -1;
        }
    }

    -2
}

fn identify_uncommon(ry: i32, atoms: &[AtomRec]) -> char {
    if ry == 1 {
        if has_atom(atoms, "N2") {
            return 'g';
        }
        return 'a';
    }
    if ry == 0 {
        let c5m = has_atom(atoms, "C5M");
        let n4 = has_atom(atoms, "N4");
        let o4 = has_atom(atoms, "O4");
        let o2p = has_atom(atoms, "O2'");
        if !o2p && (c5m || (c5m && o4)) {
            return 't';
        }
        if n4 {
            return 'c';
        }
        return 'u';
    }
    '?'
}

fn base_from_resname_or_infer(resname: &str, ry: i32, atoms: &[AtomRec]) -> char {
    match resname {
        "A" | "ADE" => 'A',
        "G" | "GUA" => 'G',
        "U" | "URA" => 'U',
        "C" | "CYT" => 'C',
        "T" | "THY" => 'T',
        "I" | "INO" => 'I',
        "P" | "PSU" => 'P',
        _ => identify_uncommon(ry, atoms),
    }
}

pub(crate) fn parse_structure_nucleic_residues(
    path: &Path,
) -> std::io::Result<Vec<BaseResidue>> {
    let mut opts = ReadOptions::new();
    opts.set_level(StrictnessLevel::Loose)
        .set_discard_hydrogens(true)
        .set_only_first_model(true)
        .set_only_atomic_coords(true)
        .set_capitalise_chains(false);
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        match ext.to_ascii_lowercase().as_str() {
            "pdb" | "pdb1" | "ent" => {
                opts.set_format(Format::Pdb);
            }
            "cif" | "mmcif" => {
                opts.set_format(Format::Mmcif);
            }
            _ => {}
        }
    }

    let path_str = path.to_string_lossy();
    let (pdb, _warnings) = opts.read(path_str.as_ref()).map_err(|errs| {
        let msg = errs
            .iter()
            .map(|e| e.to_string())
            .collect::<Vec<_>>()
            .join("\n");
        std::io::Error::new(std::io::ErrorKind::InvalidData, msg)
    })?;

    let mut groups: HashMap<ResidueKey, BTreeMap<usize, AtomRec>> = HashMap::new();
    for h in pdb.atoms_with_hierarchy() {
        let chain_id = h.chain().id().to_string();
        let resseq = i32::try_from(h.residue().serial_number()).unwrap_or(9999);
        let insertion_code = h
            .residue()
            .insertion_code()
            .and_then(|s| s.chars().next());
        let resname = canonical_residue_name(h.conformer().name());
        if resname == "HOH" || resname == "WAT" {
            continue;
        }

        let key = ResidueKey {
            chain_id,
            resseq,
            insertion_code,
            resname,
        };
        let serial = h.atom().serial_number();
        let atoms = groups.entry(key).or_insert_with(BTreeMap::new);
        atoms.entry(serial).or_insert_with(|| AtomRec {
            name: canonical_atom_name(h.atom().name()),
            x: h.atom().x(),
            y: h.atom().y(),
            z: h.atom().z(),
        });
    }

    let mut by_res_id: HashMap<(String, i32, Option<char>), Vec<ResidueKey>> = HashMap::new();
    for key in groups.keys() {
        by_res_id
            .entry((key.chain_id.clone(), key.resseq, key.insertion_code))
            .or_default()
            .push(key.clone());
    }

    let mut order_key: HashMap<ResidueKey, usize> = HashMap::new();
    for (_res_id, keys) in &by_res_id {
        if keys.len() <= 1 {
            continue;
        }
        let mut serial_count: HashMap<usize, usize> = HashMap::new();
        for key in keys {
            if let Some(atoms) = groups.get(key) {
                for serial in atoms.keys() {
                    *serial_count.entry(*serial).or_insert(0) += 1;
                }
            }
        }

        for key in keys {
            let Some(atoms) = groups.get(key) else {
                continue;
            };
            let Some(overall_min) = atoms.keys().next().copied() else {
                continue;
            };
            let exclusive_min = atoms
                .keys()
                .filter(|s| serial_count.get(s).copied().unwrap_or(0) == 1)
                .min()
                .copied();
            order_key.insert(key.clone(), exclusive_min.unwrap_or(overall_min));
        }

        let mut best_key: Option<ResidueKey> = None;
        let mut best_score: Option<(usize, usize, String)> = None;
        for key in keys {
            let Some(atoms) = groups.get(key) else {
                continue;
            };
            let Some(overall_min) = atoms.keys().next().copied() else {
                continue;
            };
            let exclusive_min = atoms
                .keys()
                .filter(|s| serial_count.get(s).copied().unwrap_or(0) == 1)
                .min()
                .copied();
            let priority = if exclusive_min.is_some() { 0usize } else { 1usize };
            let primary = exclusive_min.unwrap_or(overall_min);
            let score = (priority, primary, key.resname.clone());
            let is_better = match best_score.as_ref() {
                None => true,
                Some(b) => score < *b,
            };
            if is_better {
                best_score = Some(score);
                best_key = Some(key.clone());
            }
        }

        let Some(owner) = best_key else {
            continue;
        };

        for key in keys {
            if key == &owner {
                continue;
            }
            if let Some(atoms) = groups.get_mut(key) {
                atoms.retain(|serial, _| serial_count.get(serial).copied().unwrap_or(0) == 1);
            }
        }
    }

    groups.retain(|_, atoms| !atoms.is_empty());

    let mut entries: Vec<(usize, ResidueKey, Vec<AtomRec>)> = Vec::new();
    for (key, atoms_by_serial) in groups {
        let min_serial = order_key
            .get(&key)
            .copied()
            .unwrap_or_else(|| *atoms_by_serial.keys().next().expect("non-empty"));
        let atoms: Vec<AtomRec> = atoms_by_serial.into_iter().map(|(_, a)| a).collect();
        entries.push((min_serial, key, atoms));
    }
    entries.sort_by(|a, b| {
        (
            a.0,
            &a.1.chain_id,
            a.1.resseq,
            &a.1.resname,
        )
            .cmp(&(
                b.0,
                &b.1.chain_id,
                b.1.resseq,
                &b.1.resname,
            ))
    });

    let mut out: Vec<BaseResidue> = Vec::new();
    for (_min_serial, key, atoms) in entries {
        let ry = residue_ident(&atoms);
        if ry < 0 {
            continue;
        }
        let base = base_from_resname_or_infer(&key.resname, ry, &atoms);
        let chain = key.chain_id.chars().next().unwrap_or(' ');
        out.push(BaseResidue {
            chain,
            resseq: key.resseq,
            insertion_code: key.insertion_code,
            base,
        });
    }

    Ok(out)
}

pub fn parse_structure_bases(path: &Path) -> std::io::Result<Vec<BaseResidue>> {
    let residues = parse_structure_nucleic_residues(path)?;
    let mut out: Vec<BaseResidue> = Vec::new();

    let mut idx = 0usize;
    while idx < residues.len() {
        let chain = residues[idx].chain;
        let start = idx;
        idx += 1;
        while idx < residues.len() && residues[idx].chain == chain {
            idx += 1;
        }
        if idx - start > 1 {
            out.extend_from_slice(&residues[start..idx]);
        }
    }

    Ok(out)
}
