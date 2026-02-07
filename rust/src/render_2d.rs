use crate::semantics::StructurePolicies;
use crate::PairsJson;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Render2dOutput {
    pub xml: String,
    pub ps: String,
}

pub fn render_2d_from_pairs_json(
    pairs: &PairsJson,
    source_path: &Path,
    structure_policies: &StructurePolicies,
    basepars_ps_image: &Path,
) -> Result<Render2dOutput, String> {
    let mut arrays =
        crate::noc_engine::build_legacy_arrays(source_path, structure_policies)
            .map_err(|e| e.to_string())?;

    let sugar_syn = crate::legacy_pairing::syn_or_anti(
        arrays.num_residue,
        &arrays.atom_name,
        &arrays.seidx,
        &arrays.xyz,
        &arrays.ry,
    );

    let layout = crate::legacy_2d_layout::compute_layout_2d(
        arrays.num_residue,
        &arrays.bseq,
        &arrays.seidx,
        &arrays.ry,
        &arrays.atom_name,
        &arrays.chain_id,
        &mut arrays.xyz,
        pairs,
    )?;

    let mut base_pairs_out_order = pairs.core.base_pairs.clone();
    base_pairs_out_order.sort_by_key(|bp| bp.out_index.unwrap_or(i32::MAX));

    let xml = crate::legacy_rnaml::write_rnaml_xml(
        pairs.source.as_ref(),
        source_path,
        arrays.num_residue,
        &arrays.bseq,
        &arrays.seidx,
        &arrays.atom_name,
        &arrays.resname,
        &arrays.chain_id,
        &arrays.resseq,
        &arrays.xyz,
        &sugar_syn,
        &layout,
        &base_pairs_out_order,
    )?;
    let ps = crate::legacy_xml2ps::ps_from_rnaml_xml(&xml, basepars_ps_image)?;
    Ok(Render2dOutput { xml, ps })
}
