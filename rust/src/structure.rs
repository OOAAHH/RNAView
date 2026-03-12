use pdbtbx::{
    ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain, Element,
    Format, ReadOptions, StrictnessLevel,
};
use crate::semantics::{ChainIdPolicy, HydrogenPolicy, MissingInsertionCodePolicy, StructurePolicies};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseResidue {
    pub chain: char,
    pub resseq: i32,
    pub insertion_code: Option<char>,
    pub base: char,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum PolymerClass {
    Rna,
    Dna,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum SugarClass {
    Ribose,
    Deoxyribose,
    Unknown,
}

impl SugarClass {
    pub const fn legacy_code(self) -> u8 {
        match self {
            Self::Unknown => 0,
            Self::Ribose => 1,
            Self::Deoxyribose => 2,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct ResidueKey {
    chain_id: String,
    resseq: i32,
    insertion_code: Option<char>,
    resname: String,
}

#[derive(Debug, Clone)]
pub(crate) struct AtomRec {
    pub(crate) name: String,
    pub(crate) x: f64,
    pub(crate) y: f64,
    pub(crate) z: f64,
}

#[derive(Debug, Clone)]
pub(crate) struct ResidueAtoms {
    pub chain: char,
    pub chain_full: String,
    pub resseq: i32,
    pub insertion_code: Option<char>,
    pub resname: String,
    pub raw_resname: String,
    pub ry: i32,
    pub base: char,
    pub polymer_class: PolymerClass,
    pub sugar_class: SugarClass,
    pub is_modified: bool,
    pub atoms: Vec<AtomRec>,
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

fn legacy_mmcif_should_discard_h_atom(raw: &str) -> bool {
    let s = raw.trim();
    if s.is_empty() {
        return false;
    }
    let len = s.chars().count();
    let idx = if len == 4 { 1 } else { 0 };
    let Some(ch) = s.chars().nth(idx) else {
        return false;
    };
    ch.to_ascii_uppercase() == 'H'
}

fn is_hydrogen_atom(raw_name: &str, element: Option<&str>) -> bool {
    if let Some(sym) = element.map(str::trim) {
        if sym.eq_ignore_ascii_case("H") || sym.eq_ignore_ascii_case("D") {
            return true;
        }
    }

    let s = raw_name.trim();
    if s.is_empty() {
        return false;
    }
    let mut it = s.chars();
    while let Some(ch) = it.next() {
        if ch.is_ascii_digit() {
            continue;
        }
        let u = ch.to_ascii_uppercase();
        return u == 'H' || u == 'D';
    }
    false
}

fn normalize_insertion_code(
    is_mmcif: bool,
    insertion_code: Option<char>,
    missing_policy: MissingInsertionCodePolicy,
) -> Option<char> {
    if !is_mmcif {
        return insertion_code;
    }
    match missing_policy {
        MissingInsertionCodePolicy::LegacyQuestionMark => {
            if insertion_code.is_none() || insertion_code == Some('.') {
                Some('?')
            } else {
                insertion_code
            }
        }
        MissingInsertionCodePolicy::None => {
            if insertion_code.is_none() || insertion_code == Some('.') || insertion_code == Some('?') {
                None
            } else {
                insertion_code
            }
        }
    }
}

fn build_chain_id_map<'a>(
    chain_ids: impl IntoIterator<Item = &'a str>,
    policy: ChainIdPolicy,
) -> std::io::Result<HashMap<String, char>> {
    let mut unique: std::collections::BTreeSet<&'a str> = std::collections::BTreeSet::new();
    for id in chain_ids {
        unique.insert(id);
    }

    match policy {
        ChainIdPolicy::Legacy1Char => {
            let mut out = HashMap::new();
            for id in unique {
                out.insert(id.to_string(), id.chars().next().unwrap_or(' '));
            }
            Ok(out)
        }
        ChainIdPolicy::Unique1Char => {
            // `.out` has a 1-char chain field; for science-v1 we still need to avoid collisions.
            // Expand the pool beyond alnum so large mmCIF structures (many distinct chain IDs)
            // can still run deterministically.
            //
            // Intentionally excludes `,` and `:` since those are structural separators in `.out`.
            const POOL: &[u8] =
                b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!#$%&()*+-./;<=>?@[]^_{|}~";

            let mut out: HashMap<String, char> = HashMap::new();
            let mut used: std::collections::HashSet<char> = std::collections::HashSet::new();
            let mut pending: Vec<&'a str> = Vec::new();

            for id in unique {
                if id.chars().count() == 1 {
                    let c = id.chars().next().unwrap_or(' ');
                    if used.insert(c) {
                        out.insert(id.to_string(), c);
                    } else {
                        pending.push(id);
                    }
                } else {
                    pending.push(id);
                }
            }

            let available: Vec<char> = POOL
                .iter()
                .map(|b| *b as char)
                .filter(|c| !used.contains(c))
                .collect();
            let mut pool = available.into_iter();

            for id in pending {
                let Some(c) = pool.next() else {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "too many distinct chains for --chain-id-policy unique-1char",
                    ));
                };
                used.insert(c);
                out.insert(id.to_string(), c);
            }

            Ok(out)
        }
    }
}

fn canonical_residue_name(raw: &str) -> String {
    raw.trim().to_ascii_uppercase()
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

fn is_standard_rna_resname(resname: &str) -> bool {
    matches!(resname, "A" | "ADE" | "G" | "GUA" | "U" | "URA" | "C" | "CYT" | "I" | "INO")
}

fn is_standard_dna_resname(resname: &str) -> bool {
    matches!(resname, "DA" | "DG" | "DC" | "DT" | "DI" | "T" | "THY")
}

fn legacy_resname_from_raw_resname(resname: &str) -> String {
    match resname {
        "DA" | "A" | "ADE" => "A".to_string(),
        "DG" | "G" | "GUA" => "G".to_string(),
        "DC" | "C" | "CYT" => "C".to_string(),
        "DT" | "T" | "THY" => "T".to_string(),
        "DU" | "U" | "URA" => "U".to_string(),
        "DI" | "I" | "INO" => "I".to_string(),
        "PSU" | "P" => "P".to_string(),
        other => other.to_string(),
    }
}

fn base_from_resname_or_infer(resname: &str, ry: i32, atoms: &[AtomRec]) -> char {
    match resname {
        "DA" | "A" | "ADE" => 'A',
        "DG" | "G" | "GUA" => 'G',
        "DU" | "U" | "URA" => 'U',
        "DC" | "C" | "CYT" => 'C',
        "DT" | "T" | "THY" => 'T',
        "DI" | "I" | "INO" => 'I',
        "P" | "PSU" => 'P',
        _ => identify_uncommon(ry, atoms),
    }
}

fn classify_polymer_and_sugar(resname: &str, atoms: &[AtomRec]) -> (PolymerClass, SugarClass) {
    match resname {
        "DA" | "DG" | "DC" | "DT" | "DI" | "T" | "THY" => {
            return (PolymerClass::Dna, SugarClass::Deoxyribose);
        }
        "A" | "ADE" | "G" | "GUA" | "C" | "CYT" | "U" | "URA" | "I" | "INO" | "P" | "PSU" => {
            return (PolymerClass::Rna, SugarClass::Ribose);
        }
        _ => {}
    }

    if has_atom(atoms, "O2'") {
        return (PolymerClass::Rna, SugarClass::Ribose);
    }

    let backbone_markers = ["C1'", "C2'", "C3'", "C4'", "O3'", "O4'"];
    let has_backbone = backbone_markers.iter().any(|name| has_atom(atoms, name));
    if has_backbone {
        return (PolymerClass::Dna, SugarClass::Deoxyribose);
    }

    (PolymerClass::Unknown, SugarClass::Unknown)
}

fn is_modified_residue_name(resname: &str) -> bool {
    !matches!(
        resname,
        "A" | "ADE" | "G" | "GUA" | "C" | "CYT" | "U" | "URA" | "I" | "INO" | "T" | "THY" | "DA" | "DG" | "DC" | "DT" | "DI"
    )
}

pub(crate) fn parse_structure_nucleic_residues(
    path: &Path,
) -> std::io::Result<Vec<BaseResidue>> {
    let structure_policies = StructurePolicies::legacy_v1_defaults();
    parse_structure_nucleic_residues_with_policies(path, &structure_policies)
}

pub(crate) fn parse_structure_nucleic_residues_with_policies(
    path: &Path,
    structure_policies: &StructurePolicies,
) -> std::io::Result<Vec<BaseResidue>> {
    let mut opts = ReadOptions::new();
    opts.set_level(StrictnessLevel::Loose)
        .set_only_first_model(true)
        .set_only_atomic_coords(true)
        .set_capitalise_chains(false);
    let mut is_mmcif = false;
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        match ext.to_ascii_lowercase().as_str() {
            "pdb" | "pdb1" | "ent" => {
                opts.set_format(Format::Pdb);
            }
            "cif" | "mmcif" => {
                opts.set_format(Format::Mmcif);
                is_mmcif = true;
            }
            _ => {}
        }
    }

    let apply_legacy_mmcif_h_bug_filter =
        is_mmcif && structure_policies.hydrogen_policy == HydrogenPolicy::LegacyMmcifBug;
    let discard_hydrogens = match structure_policies.hydrogen_policy {
        HydrogenPolicy::DiscardAll => true,
        HydrogenPolicy::KeepAll => false,
        HydrogenPolicy::LegacyMmcifBug => !is_mmcif,
    };
    opts.set_discard_hydrogens(discard_hydrogens);

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
        let atom = h.atom();
        if structure_policies.hydrogen_policy == HydrogenPolicy::DiscardAll
            && is_hydrogen_atom(atom.name(), atom.element().map(Element::symbol))
        {
            continue;
        }
        let chain_id = h.chain().id().to_string();
        let resseq = i32::try_from(h.residue().serial_number()).unwrap_or(9999);
        let insertion_code = normalize_insertion_code(
            is_mmcif,
            h.residue()
                .insertion_code()
                .and_then(|s| s.chars().next()),
            structure_policies.missing_insertion_code_policy,
        );
        let raw_resname = canonical_residue_name(h.conformer().name());
        if raw_resname == "HOH" || raw_resname == "WAT" {
            continue;
        }

        let key = ResidueKey {
            chain_id,
            resseq,
            insertion_code,
            resname: raw_resname,
        };
        let serial = atom.serial_number();
        if apply_legacy_mmcif_h_bug_filter && legacy_mmcif_should_discard_h_atom(atom.name()) {
            continue;
        }
        let name = canonical_atom_name(atom.name());
        let atoms = groups.entry(key).or_insert_with(BTreeMap::new);
        if atoms.values().any(|a| a.name == name) {
            continue;
        }
        atoms.entry(serial).or_insert_with(|| AtomRec {
            name,
            x: atom.x(),
            y: atom.y(),
            z: atom.z(),
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

    // Chain mapping only needs to consider residues that actually participate in base indexing
    // (i.e., nucleic residues). Including protein chains can exceed the 1-char pool under
    // `unique-1char` for large complexes (e.g. ribosomes).
    let mut nucleic: Vec<(ResidueKey, char)> = Vec::new();
    for (_min_serial, key, atoms) in entries {
        let ry = residue_ident(&atoms);
        if ry < 0 {
            continue;
        }
        let base = base_from_resname_or_infer(&key.resname, ry, &atoms);
        nucleic.push((key, base));
    }

    let chain_map = build_chain_id_map(
        nucleic.iter().map(|(key, _base)| key.chain_id.as_str()),
        structure_policies.chain_id_policy,
    )?;

    let mut out: Vec<BaseResidue> = Vec::new();
    for (key, base) in nucleic {
        let chain = chain_map.get(key.chain_id.as_str()).copied().unwrap_or(' ');
        out.push(BaseResidue {
            chain,
            resseq: key.resseq,
            insertion_code: key.insertion_code,
            base,
        });
    }

    Ok(out)
}

pub(crate) fn parse_structure_nucleic_residues_with_atoms(path: &Path) -> std::io::Result<Vec<ResidueAtoms>> {
    let structure_policies = StructurePolicies::legacy_v1_defaults();
    parse_structure_nucleic_residues_with_atoms_with_policies(path, &structure_policies)
}

pub(crate) fn parse_structure_nucleic_residues_with_atoms_with_policies(
    path: &Path,
    structure_policies: &StructurePolicies,
) -> std::io::Result<Vec<ResidueAtoms>> {
    let mut opts = ReadOptions::new();
    opts.set_level(StrictnessLevel::Loose)
        .set_only_first_model(true)
        .set_only_atomic_coords(true)
        .set_capitalise_chains(false);
    let mut is_mmcif = false;
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        match ext.to_ascii_lowercase().as_str() {
            "pdb" | "pdb1" | "ent" => {
                opts.set_format(Format::Pdb);
            }
            "cif" | "mmcif" => {
                opts.set_format(Format::Mmcif);
                is_mmcif = true;
            }
            _ => {}
        }
    }

    let apply_legacy_mmcif_h_bug_filter =
        is_mmcif && structure_policies.hydrogen_policy == HydrogenPolicy::LegacyMmcifBug;
    let discard_hydrogens = match structure_policies.hydrogen_policy {
        HydrogenPolicy::DiscardAll => true,
        HydrogenPolicy::KeepAll => false,
        HydrogenPolicy::LegacyMmcifBug => !is_mmcif,
    };
    opts.set_discard_hydrogens(discard_hydrogens);
    if apply_legacy_mmcif_h_bug_filter {
        // Legacy mmCIF parsing attempted to remove hydrogens but effectively kept some due to an
        // indexing bug. To reproduce that behavior, we must read H atoms and apply the same
        // buggy filter ourselves.
        opts.set_discard_hydrogens(false);
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
        let atom = h.atom();
        if structure_policies.hydrogen_policy == HydrogenPolicy::DiscardAll
            && is_hydrogen_atom(atom.name(), atom.element().map(Element::symbol))
        {
            continue;
        }
        let chain_id = h.chain().id().to_string();
        let resseq = i32::try_from(h.residue().serial_number()).unwrap_or(9999);
        let insertion_code = normalize_insertion_code(
            is_mmcif,
            h.residue()
                .insertion_code()
                .and_then(|s| s.chars().next()),
            structure_policies.missing_insertion_code_policy,
        );
        let raw_resname = canonical_residue_name(h.conformer().name());
        if raw_resname == "HOH" || raw_resname == "WAT" {
            continue;
        }

        let key = ResidueKey {
            chain_id,
            resseq,
            insertion_code,
            resname: raw_resname,
        };
        let serial = atom.serial_number();
        if apply_legacy_mmcif_h_bug_filter && legacy_mmcif_should_discard_h_atom(atom.name()) {
            continue;
        }
        let name = canonical_atom_name(atom.name());
        let atoms = groups.entry(key).or_insert_with(BTreeMap::new);
        if atoms.values().any(|a| a.name == name) {
            continue;
        }
        atoms.entry(serial).or_insert_with(|| AtomRec {
            name,
            x: atom.x(),
            y: atom.y(),
            z: atom.z(),
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

    // Chain mapping only needs to consider residues that actually participate in base indexing
    // (i.e., nucleic residues). Including protein chains can exceed the 1-char pool under
    // `unique-1char` for large complexes (e.g. ribosomes).
    let mut nucleic: Vec<(ResidueKey, Vec<AtomRec>, i32, char)> = Vec::new();
    for (_min_serial, key, atoms) in entries {
        let ry = residue_ident(&atoms);
        if ry < 0 {
            continue;
        }
        let base = base_from_resname_or_infer(&key.resname, ry, &atoms);
        nucleic.push((key, atoms, ry, base));
    }

    let chain_map = build_chain_id_map(
        nucleic.iter().map(|(key, _atoms, _ry, _base)| key.chain_id.as_str()),
        structure_policies.chain_id_policy,
    )?;

    let mut out: Vec<ResidueAtoms> = Vec::new();
    for (key, atoms, ry, base) in nucleic {
        let chain = chain_map.get(key.chain_id.as_str()).copied().unwrap_or(' ');
        let raw_resname = key.resname;
        let (polymer_class, sugar_class) = classify_polymer_and_sugar(&raw_resname, &atoms);
        let legacy_resname = legacy_resname_from_raw_resname(&raw_resname);
        let is_modified = is_modified_residue_name(&raw_resname);
        out.push(ResidueAtoms {
            chain,
            chain_full: key.chain_id.clone(),
            resseq: key.resseq,
            insertion_code: key.insertion_code,
            resname: legacy_resname,
            raw_resname,
            ry,
            base,
            polymer_class,
            sugar_class,
            is_modified,
            atoms,
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

pub(crate) fn parse_structure_bases_with_atoms(path: &Path) -> std::io::Result<Vec<ResidueAtoms>> {
    let structure_policies = StructurePolicies::legacy_v1_defaults();
    parse_structure_bases_with_atoms_with_policies(path, &structure_policies)
}

pub(crate) fn parse_structure_bases_with_atoms_with_policies(
    path: &Path,
    structure_policies: &StructurePolicies,
) -> std::io::Result<Vec<ResidueAtoms>> {
    let residues = parse_structure_nucleic_residues_with_atoms_with_policies(path, structure_policies)?;
    let mut out: Vec<ResidueAtoms> = Vec::new();

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_mmcif_h_bug_filter_is_bug_compatible() {
        assert!(legacy_mmcif_should_discard_h_atom("H1"));
        assert!(legacy_mmcif_should_discard_h_atom(" H2 "));
        assert!(!legacy_mmcif_should_discard_h_atom("H5''"));
        assert!(!legacy_mmcif_should_discard_h_atom("1H2"));
    }

    #[test]
    fn hydrogen_detection_uses_element_or_name() {
        assert!(is_hydrogen_atom("H5''", None));
        assert!(is_hydrogen_atom("1H2", None));
        assert!(is_hydrogen_atom("D1", None));
        assert!(is_hydrogen_atom("CA", Some("H")));
        assert!(is_hydrogen_atom("CA", Some("D")));
        assert!(!is_hydrogen_atom("CA", Some("C")));
        assert!(!is_hydrogen_atom("CA", None));
    }

    #[test]
    fn insertion_code_normalization_matches_policies() {
        assert_eq!(
            normalize_insertion_code(true, None, MissingInsertionCodePolicy::LegacyQuestionMark),
            Some('?')
        );
        assert_eq!(
            normalize_insertion_code(true, Some('.'), MissingInsertionCodePolicy::LegacyQuestionMark),
            Some('?')
        );
        assert_eq!(
            normalize_insertion_code(true, Some('A'), MissingInsertionCodePolicy::LegacyQuestionMark),
            Some('A')
        );

        assert_eq!(
            normalize_insertion_code(true, None, MissingInsertionCodePolicy::None),
            None
        );
        assert_eq!(
            normalize_insertion_code(true, Some('.'), MissingInsertionCodePolicy::None),
            None
        );
        assert_eq!(
            normalize_insertion_code(true, Some('?'), MissingInsertionCodePolicy::None),
            None
        );
        assert_eq!(
            normalize_insertion_code(true, Some('A'), MissingInsertionCodePolicy::None),
            Some('A')
        );

        assert_eq!(
            normalize_insertion_code(false, Some('A'), MissingInsertionCodePolicy::None),
            Some('A')
        );
        assert_eq!(
            normalize_insertion_code(false, None, MissingInsertionCodePolicy::LegacyQuestionMark),
            None
        );
    }

    #[test]
    fn chain_id_mapping_unique_1char_is_deterministic_and_collision_free() {
        let ids = ["AB", "AA", "A", "B"];

        let legacy = build_chain_id_map(ids.iter().copied(), ChainIdPolicy::Legacy1Char).expect("legacy map");
        assert_eq!(legacy.get("AA").copied(), Some('A'));
        assert_eq!(legacy.get("AB").copied(), Some('A'));

        let unique = build_chain_id_map(ids.iter().copied(), ChainIdPolicy::Unique1Char).expect("unique map");
        let aa = unique.get("AA").copied().unwrap_or('?');
        let ab = unique.get("AB").copied().unwrap_or('?');
        assert_ne!(aa, ab, "unique mapping must avoid collisions");

        let unique2 = build_chain_id_map(ids.iter().copied(), ChainIdPolicy::Unique1Char).expect("unique map 2");
        assert_eq!(unique, unique2, "unique mapping must be deterministic");
    }

    #[test]
    fn raw_residue_name_is_preserved_while_legacy_name_is_collapsed() {
        assert_eq!(canonical_residue_name("dt"), "DT");
        assert_eq!(legacy_resname_from_raw_resname("DT"), "T");
        assert_eq!(legacy_resname_from_raw_resname("DA"), "A");
        assert_eq!(legacy_resname_from_raw_resname("PSU"), "P");
    }

    #[test]
    fn base_inference_understands_standard_dna_residue_names() {
        let atoms = vec![
            AtomRec { name: "N1".to_string(), x: 0.0, y: 0.0, z: 0.0 },
            AtomRec { name: "C2".to_string(), x: 1.3, y: 0.0, z: 0.0 },
            AtomRec { name: "C6".to_string(), x: -1.3, y: 0.0, z: 0.0 },
        ];
        assert_eq!(base_from_resname_or_infer("DA", 0, &atoms), 'A');
        assert_eq!(base_from_resname_or_infer("DT", 0, &atoms), 'T');
        assert_eq!(base_from_resname_or_infer("DC", 0, &atoms), 'C');
        assert_eq!(base_from_resname_or_infer("DG", 1, &atoms), 'G');
    }

    #[test]
    fn polymer_classifier_uses_explicit_names_before_atom_fallback() {
        let dna_atoms = vec![
            AtomRec { name: "C1'".to_string(), x: 0.0, y: 0.0, z: 0.0 },
            AtomRec { name: "C2'".to_string(), x: 0.0, y: 0.0, z: 0.0 },
        ];
        let rna_atoms = vec![
            AtomRec { name: "C1'".to_string(), x: 0.0, y: 0.0, z: 0.0 },
            AtomRec { name: "O2'".to_string(), x: 0.0, y: 0.0, z: 0.0 },
        ];

        assert_eq!(
            classify_polymer_and_sugar("DT", &rna_atoms),
            (PolymerClass::Dna, SugarClass::Deoxyribose)
        );
        assert_eq!(
            classify_polymer_and_sugar("PSU", &rna_atoms),
            (PolymerClass::Rna, SugarClass::Ribose)
        );
        assert_eq!(
            classify_polymer_and_sugar("XAA", &dna_atoms),
            (PolymerClass::Dna, SugarClass::Deoxyribose)
        );
    }

    #[test]
    fn modified_flag_distinguishes_standard_and_modified_residues() {
        assert!(!is_modified_residue_name("DT"));
        assert!(!is_modified_residue_name("A"));
        assert!(is_modified_residue_name("PSU"));
        assert!(is_modified_residue_name("1MG"));
    }
}
