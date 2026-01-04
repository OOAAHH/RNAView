use pdbtbx::{
    ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain,
    Format, ReadOptions, StrictnessLevel,
};
use std::collections::{BTreeMap, HashMap, HashSet};
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
    chain: char,
    resseq: i32,
    insertion_code: Option<char>,
    resname: String,
}

#[derive(Debug, Clone)]
struct AtomRec {
    serial: usize,
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

fn is_explicit_nucleotide_resname(resname: &str) -> bool {
    matches!(
        resname,
        "A" | "ADE" | "G" | "GUA" | "U" | "URA" | "C" | "CYT" | "T" | "THY" | "I" | "INO"
            | "P"
            | "PSU"
    )
}

fn parse_structure_to_nucleic_residues(
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
        let chain = h.chain().id().chars().next().unwrap_or(' ');
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
            chain,
            resseq,
            insertion_code,
            resname,
        };
        let serial = h.atom().serial_number();
        let atoms = groups.entry(key).or_insert_with(BTreeMap::new);
        atoms.entry(serial).or_insert_with(|| AtomRec {
            serial,
            name: canonical_atom_name(h.atom().name()),
            x: h.atom().x(),
            y: h.atom().y(),
            z: h.atom().z(),
        });
    }

    let mut explicit_nuc_serials: HashMap<(char, i32, Option<char>), HashSet<usize>> =
        HashMap::new();
    for (key, atoms) in &groups {
        if !is_explicit_nucleotide_resname(&key.resname) {
            continue;
        }
        explicit_nuc_serials
            .entry((key.chain, key.resseq, key.insertion_code))
            .or_default()
            .extend(atoms.keys().copied());
    }

    let mut entries: Vec<(usize, ResidueKey, Vec<AtomRec>)> = Vec::new();
    for (key, mut atoms_by_serial) in groups {
        if !is_explicit_nucleotide_resname(&key.resname) {
            if let Some(nuc) = explicit_nuc_serials.get(&(key.chain, key.resseq, key.insertion_code))
            {
                if !nuc.is_empty() {
                    atoms_by_serial.retain(|serial, _| !nuc.contains(serial));
                }
            }
        }
        if atoms_by_serial.is_empty() {
            continue;
        }
        let min_serial = *atoms_by_serial.keys().next().expect("non-empty");
        let atoms: Vec<AtomRec> = atoms_by_serial.into_iter().map(|(_, a)| a).collect();
        entries.push((min_serial, key, atoms));
    }
    entries.sort_by_key(|(min_serial, _, _)| *min_serial);

    let mut out: Vec<BaseResidue> = Vec::new();
    for (_min_serial, key, atoms) in entries {
        let ry = residue_ident(&atoms);
        if ry < 0 {
            continue;
        }
        let base = base_from_resname_or_infer(&key.resname, ry, &atoms);
        out.push(BaseResidue {
            chain: key.chain,
            resseq: key.resseq,
            insertion_code: key.insertion_code,
            base,
        });
    }

    Ok(out)
}

pub fn parse_structure_bases(path: &Path) -> std::io::Result<Vec<BaseResidue>> {
    let residues = parse_structure_to_nucleic_residues(path)?;
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
