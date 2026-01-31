// Legacy RNAML -> PostScript renderer.
//
// Ported from legacy C: src/xml2ps.c.

use std::path::Path;

pub(crate) fn ps_from_rnaml_xml(_xml: &str, _basepars_ps_image: &Path) -> Result<String, String> {
    Err("ps_from_rnaml_xml: not implemented".to_string())
}

