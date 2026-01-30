use crate::semantics::StructurePolicies;
use crate::PairsJson;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Render2dOutput {
    pub xml: String,
    pub ps: String,
}

pub fn render_2d_from_pairs_json(
    _pairs: &PairsJson,
    _source_path: &Path,
    _structure_policies: &StructurePolicies,
    _basepars_ps_image: &Path,
) -> Result<Render2dOutput, String> {
    Err("render_2d_from_pairs_json: not implemented".to_string())
}

