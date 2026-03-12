use crate::legacy_alg::{compute_base_info, AtomName4, ResName3};
use crate::legacy_pairing::{all_pairs, bp_network_multiplets, pair_type_statistics, uncommon_lines};
use crate::out_full::{OutEol, OutFull};
use crate::structure::parse_structure_bases_with_atoms_with_policies;
use crate::{ChemistryPolicies, StructurePolicies};
use std::path::Path;

#[derive(Debug, Clone)]
pub(crate) struct LegacyArrays {
    pub(crate) num_residue: usize,
    pub(crate) num_atoms: usize,
    pub(crate) seidx: Vec<[usize; 3]>, // 1-based: [_, start, end]
    pub(crate) atom_name: Vec<AtomName4>, // 1-based
    pub(crate) xyz: Vec<[f64; 4]>,        // 1-based: [_, x, y, z]
    pub(crate) chain_id: Vec<u8>,         // 1-based
    pub(crate) resseq: Vec<i32>,          // 1-based
    pub(crate) icode: Vec<u8>,            // 1-based
    pub(crate) resname: Vec<ResName3>,    // 1-based (legacy field)
    pub(crate) bseq: Vec<u8>,             // 1-based
    pub(crate) ry: Vec<i32>,              // 1-based
    pub(crate) sugar_class: Vec<u8>,      // 1-based: 0 unknown, 1 ribose, 2 deoxyribose
}

fn atom_name_field_legacy(name: &str) -> AtomName4 {
    let mut s = name.trim().to_ascii_uppercase();
    if s.len() > 4 {
        s = s.chars().take(4).collect();
    }
    let mut out = String::with_capacity(4);
    if s.len() < 4 {
        out.push(' ');
    }
    out.push_str(&s);
    while out.len() < 4 {
        out.push(' ');
    }
    let b = out.as_bytes();
    AtomName4([b[0], b[1], b[2], b[3]])
}

fn residue_name_field_legacy(name: &str) -> ResName3 {
    let mut s = name.trim().to_ascii_uppercase();
    if s.len() > 3 {
        s = s.chars().take(3).collect();
    }
    let field = format!("{s:>3}");
    let b = field.as_bytes();
    ResName3([b[0], b[1], b[2]])
}

pub(crate) fn build_legacy_arrays(path: &Path, structure_policies: &StructurePolicies) -> std::io::Result<LegacyArrays> {
    let residues = parse_structure_bases_with_atoms_with_policies(path, structure_policies)?;
    let num_residue = residues.len();

    let mut seidx: Vec<[usize; 3]> = vec![[0; 3]; num_residue + 1];
    let mut atom_name: Vec<AtomName4> = vec![AtomName4([0; 4])];
    let mut xyz: Vec<[f64; 4]> = vec![[0.0; 4]];
    let mut chain_id: Vec<u8> = vec![b' '];
    let mut resseq: Vec<i32> = vec![0];
    let mut icode: Vec<u8> = vec![b' '];
    let mut resname: Vec<ResName3> = vec![ResName3([b' '; 3])];

    let mut bseq: Vec<u8> = vec![b'?'; num_residue + 1];
    let mut ry: Vec<i32> = vec![-1; num_residue + 1];
    let mut sugar_class: Vec<u8> = vec![0; num_residue + 1];

    for (idx0, residue) in residues.iter().enumerate() {
        let i = idx0 + 1;
        bseq[i] = residue.base as u8;
        ry[i] = residue.ry;
        sugar_class[i] = residue.sugar_class.legacy_code();

        let start = atom_name.len();
        let resname3 = residue_name_field_legacy(&residue.resname);
        for atom in &residue.atoms {
            atom_name.push(atom_name_field_legacy(&atom.name));
            xyz.push([0.0, atom.x, atom.y, atom.z]);
            chain_id.push(residue.chain as u8);
            resseq.push(residue.resseq);
            icode.push(residue.insertion_code.unwrap_or(' ') as u8);
            resname.push(resname3);
        }
        let end = atom_name.len() - 1;
        seidx[i][1] = start;
        seidx[i][2] = end;
    }

    Ok(LegacyArrays {
        num_residue,
        num_atoms: atom_name.len() - 1,
        seidx,
        atom_name,
        xyz,
        chain_id,
        resseq,
        icode,
        resname,
        bseq,
        ry,
        sugar_class,
    })
}

pub(crate) fn compute_out_full_from_structure(
    input_path: &Path,
    pdb_data_file_name: String,
    structure_policies: &StructurePolicies,
    chemistry_policies: Option<&ChemistryPolicies>,
) -> Result<OutFull, String> {
    let arrays = build_legacy_arrays(input_path, structure_policies).map_err(|e| e.to_string())?;
    let base_info = compute_base_info(
        arrays.num_residue,
        &arrays.bseq,
        &arrays.seidx,
        &arrays.ry,
        &arrays.atom_name,
        &arrays.xyz,
    )?;
    let uncommon_lines = uncommon_lines(
        arrays.num_residue,
        &arrays.seidx,
        &arrays.resname,
        &arrays.chain_id,
        &arrays.resseq,
        &arrays.icode,
        &arrays.bseq,
        &arrays.ry,
    );

    let all = all_pairs(
        arrays.num_residue,
        &arrays.ry,
        &base_info.nxyz,
        &base_info.orien,
        &base_info.org,
        &base_info.bprs,
        &arrays.seidx,
        &arrays.xyz,
        &arrays.atom_name,
        &arrays.chain_id,
        &arrays.resseq,
        &arrays.bseq,
        &arrays.sugar_class,
        chemistry_policies,
    )?;

    let multiplets = bp_network_multiplets(
        arrays.num_residue,
        &arrays.ry,
        &arrays.seidx,
        &arrays.atom_name,
        &arrays.chain_id,
        &arrays.resseq,
        &arrays.xyz,
        &arrays.bseq,
        &arrays.sugar_class,
        &all.pair_info,
        &base_info.nxyz,
        &base_info.orien,
        &base_info.org,
        &base_info.bprs,
        chemistry_policies,
    )?;
    let has_multiplets = multiplets.is_some();

    let (type_counts_1_to_7, type_counts_8_to_13) = pair_type_statistics(&all.pair_type_codes);

    Ok(OutFull {
        eol: OutEol::Lf,
        trailing_newline: true,
        pdb_data_file_name,
        uncommon_lines,
        bprs: Some([
            base_info.bprs[1],
            base_info.bprs[2],
            base_info.bprs[3],
            base_info.bprs[4],
            base_info.bprs[5],
            base_info.bprs[6],
        ]),
        base_pairs: all.base_pairs,
        multiplets,
        blank_lines_after_end_base_pair: if has_multiplets { 1 } else { 0 },
        blank_lines_after_end_multiplets: if has_multiplets { 1 } else { 0 },
        total_base_pairs: all.num_pair_tot as i64,
        total_bases: arrays.num_residue as i64,
        type_counts_1_to_7,
        type_counts_8_to_13,
    })
}
