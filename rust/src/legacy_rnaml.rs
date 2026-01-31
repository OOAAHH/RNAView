// Legacy RNAML (XML) writer.
//
// Ported from legacy C: src/rnaxml-new.c.

use crate::{BasePair, Source};
use crate::legacy_alg::{AtomName4, ResName3};
use std::path::Path;

pub(crate) fn write_rnaml_xml(
    _source: Option<&Source>,
    _base_path: &Path,
    _num_residue: usize,
    _bseq: &[u8],
    _seidx: &[[usize; 3]],
    _atom_name: &[AtomName4],
    _resname: &[ResName3],
    _chain_id: &[u8],
    _resseq: &[i32],
    _xyz: &[[f64; 4]],
    _sugar_syn: &[i64],
    _layout: &crate::legacy_2d_layout::Layout2d,
    _base_pairs_out_order: &[BasePair],
) -> Result<String, String> {
    Err("write_rnaml_xml: not implemented".to_string())
}

