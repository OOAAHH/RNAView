use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;

mod structure;
#[cfg(feature = "legacy-ffi")]
mod legacy_ffi;
pub use structure::{parse_structure_bases, BaseResidue};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Source {
    pub path: String,
    pub format: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub id_scheme: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub model: Option<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PairsJson {
    pub schema_version: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source: Option<Source>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub options: Option<serde_json::Value>,
    pub core: Core,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Core {
    #[serde(default)]
    pub base_pairs: Vec<BasePair>,
    #[serde(default)]
    pub multiplets: Vec<Multiplet>,
    #[serde(default)]
    pub stats: Stats,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Stats {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_pairs: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_bases: Option<u32>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub pair_type_counts: BTreeMap<String, i32>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BasePair {
    pub i: i32,
    pub j: i32,
    pub chain_i: String,
    pub resseq_i: i32,
    pub base_i: String,
    pub base_j: String,
    pub resseq_j: i32,
    pub chain_j: String,
    pub kind: String, // "pair" | "stacked" | "unknown"
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lw: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub orientation: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub syn: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub text: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Multiplet {
    pub indices: Vec<i32>,
    pub text: String,
}

fn norm_ws(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

fn extract_block<'a>(lines: &'a [&'a str], begin: &str, end: &str) -> Vec<&'a str> {
    let mut in_block = false;
    let mut out = Vec::new();
    for line in lines {
        if !in_block {
            if line.trim() == begin {
                in_block = true;
            }
            continue;
        }
        if line.trim() == end {
            break;
        }
        out.push(*line);
    }
    out
}

fn parse_pair_rest(tokens: &[&str]) -> (String, Option<String>, Option<String>, i32, Option<String>) {
    if tokens.is_empty() {
        return ("pair".to_string(), None, None, 0, None);
    }

    if tokens
        .last()
        .is_some_and(|t| t.eq_ignore_ascii_case("stacked"))
    {
        let syn_count = tokens.iter().filter(|t| t.eq_ignore_ascii_case("syn")).count() as i32;
        return ("stacked".to_string(), None, None, syn_count, None);
    }

    let lw = Some(tokens[0].to_string());
    let mut orientation: Option<String> = None;
    let mut syn_count = 0i32;
    let mut note_tokens: Vec<&str> = Vec::new();

    for t in tokens.iter().skip(1) {
        if t.eq_ignore_ascii_case("cis") {
            orientation = Some("cis".to_string());
            continue;
        }
        if t.to_ascii_lowercase().starts_with("tran") {
            orientation = Some("tran".to_string());
            continue;
        }
        if t.eq_ignore_ascii_case("syn") {
            syn_count += 1;
            continue;
        }
        note_tokens.push(*t);
    }

    let note = if note_tokens.is_empty() {
        None
    } else {
        Some(note_tokens.join(" "))
    };
    ("pair".to_string(), lw, orientation, syn_count, note)
}

fn parse_base_pair_line(line: &str) -> Option<BasePair> {
    let s = line.trim();
    if s.is_empty() {
        return None;
    }

    let (ij_part, rest_part) = s.split_once(',')?;
    let (i_str, j_str) = ij_part.trim().split_once('_')?;
    let i = i_str.trim().parse::<i32>().ok()?;
    let j = j_str.trim().parse::<i32>().ok()?;

    let rest = rest_part.trim_start();
    let mut chars = rest.chars();
    let chain_i = chars.next()?.to_string();
    if chars.next()? != ':' {
        return None;
    }

    let after_chain = chars.as_str();
    let tokens: Vec<&str> = after_chain.split_whitespace().collect();
    if tokens.len() < 4 {
        return None;
    }

    let resseq_i = tokens[0].parse::<i32>().ok()?;
    let (base_i, base_j) = tokens[1].split_once('-')?;
    if base_i.len() != 1 || base_j.len() != 1 {
        return None;
    }
    let resseq_j = tokens[2].parse::<i32>().ok()?;

    let chain_j_token = tokens[3];
    let chain_j = chain_j_token
        .trim_end_matches(':')
        .chars()
        .next()
        .unwrap_or(' ')
        .to_string();

    let (kind, lw, orientation, syn_count, note) = parse_pair_rest(&tokens[4..]);
    let syn = if syn_count > 0 { Some(syn_count) } else { None };

    if kind == "stacked" {
        return Some(BasePair {
            i,
            j,
            chain_i,
            resseq_i,
            base_i: base_i.to_string(),
            base_j: base_j.to_string(),
            resseq_j,
            chain_j,
            kind,
            lw: None,
            orientation: None,
            syn,
            note: None,
            text: None,
        });
    }

    Some(BasePair {
        i,
        j,
        chain_i,
        resseq_i,
        base_i: base_i.to_string(),
        base_j: base_j.to_string(),
        resseq_j,
        chain_j,
        kind,
        lw,
        orientation,
        syn,
        note,
        text: None,
    })
}

fn extract_stats(lines: &[&str]) -> Stats {
    let total_re = Regex::new(r"(?i)^\s*The total base pairs\s*=\s*(\d+)\s*\(from\s*(\d+)\s*bases\)\s*$")
        .expect("regex total");
    let sep_re = Regex::new(r"^\s*-{5,}\s*$").expect("regex sep");

    let mut total_pairs: Option<u32> = None;
    let mut total_bases: Option<u32> = None;
    let mut total_idx: Option<usize> = None;

    for (idx, line) in lines.iter().enumerate() {
        if let Some(caps) = total_re.captures(line) {
            total_pairs = caps.get(1).and_then(|m| m.as_str().parse().ok());
            total_bases = caps.get(2).and_then(|m| m.as_str().parse().ok());
            total_idx = Some(idx);
            break;
        }
    }

    let Some(start_idx) = total_idx else {
        return Stats::default();
    };
    let mut stats = Stats::default();
    stats.total_pairs = total_pairs;
    stats.total_bases = total_bases;

    let mut pair_type_counts: BTreeMap<String, i32> = BTreeMap::new();
    let mut pending_header: Option<Vec<String>> = None;

    for line in lines.iter().skip(start_idx + 1) {
        if sep_re.is_match(line) {
            continue;
        }
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.is_empty() {
            continue;
        }
        if tokens.iter().any(|t| t.contains("--")) {
            pending_header = Some(tokens.iter().map(|t| t.to_string()).collect());
            continue;
        }
        let Some(header) = pending_header.as_ref() else {
            continue;
        };
        let all_int = tokens.iter().all(|t| {
            let t = t.strip_prefix('-').unwrap_or(t);
            !t.is_empty() && t.chars().all(|c| c.is_ascii_digit())
        });
        if !all_int {
            continue;
        }
        if tokens.len() != header.len() {
            pending_header = None;
            continue;
        }
        for (k, v) in header.iter().zip(tokens.iter()) {
            if let Ok(n) = v.parse::<i32>() {
                pair_type_counts.insert(k.clone(), n);
            }
        }
        pending_header = None;
    }

    stats.pair_type_counts = pair_type_counts;
    stats
}

fn base_pair_sort_key(bp: &BasePair) -> (i32, i32, &str, &str, i32, &str, &str, i32, &str, (&str, &str, i32, &str, &str)) {
    (
        bp.i,
        bp.j,
        bp.kind.as_str(),
        bp.chain_i.as_str(),
        bp.resseq_i,
        bp.base_i.as_str(),
        bp.base_j.as_str(),
        bp.resseq_j,
        bp.chain_j.as_str(),
        (
            bp.lw.as_deref().unwrap_or(""),
            bp.orientation.as_deref().unwrap_or(""),
            bp.syn.unwrap_or(0),
            bp.note.as_deref().unwrap_or(""),
            bp.text.as_deref().unwrap_or(""),
        ),
    )
}

pub fn extract_core_from_out_str(text: &str) -> Core {
    let raw_lines: Vec<&str> = text.lines().collect();

    let base_pair_lines = extract_block(&raw_lines, "BEGIN_base-pair", "END_base-pair");
    let mut base_pairs: Vec<BasePair> = Vec::new();
    for raw in base_pair_lines {
        if raw.trim().is_empty() {
            continue;
        }
        if let Some(parsed) = parse_base_pair_line(raw) {
            base_pairs.push(parsed);
        } else {
            base_pairs.push(BasePair {
                i: -1,
                j: -1,
                chain_i: "".to_string(),
                resseq_i: -1,
                base_i: "".to_string(),
                base_j: "".to_string(),
                resseq_j: -1,
                chain_j: "".to_string(),
                kind: "unknown".to_string(),
                lw: None,
                orientation: None,
                syn: None,
                note: None,
                text: Some(norm_ws(raw)),
            });
        }
    }

    let multiplet_lines = extract_block(&raw_lines, "BEGIN_multiplets", "END_multiplets");
    let mut multiplets: Vec<Multiplet> = Vec::new();
    for raw in multiplet_lines {
        let s = raw.trim();
        if s.is_empty() {
            continue;
        }
        if let Some((left, right)) = s.split_once('|') {
            let left = left.trim().trim_end_matches('_').trim();
            let mut idxs: Vec<i32> = Vec::new();
            for part in left.split('_') {
                let part = part.trim();
                if part.is_empty() {
                    continue;
                }
                if let Ok(n) = part.parse::<i32>() {
                    idxs.push(n);
                }
            }
            multiplets.push(Multiplet {
                indices: idxs,
                text: norm_ws(right),
            });
        } else {
            multiplets.push(Multiplet {
                indices: Vec::new(),
                text: norm_ws(s),
            });
        }
    }

    base_pairs.sort_by(|a, b| base_pair_sort_key(a).cmp(&base_pair_sort_key(b)));
    multiplets.sort_by(|a, b| (&a.indices, a.text.as_str()).cmp(&(&b.indices, b.text.as_str())));

    Core {
        base_pairs,
        multiplets,
        stats: extract_stats(&raw_lines),
    }
}

pub fn extract_core_from_out_path(path: &Path) -> std::io::Result<Core> {
    let bytes = std::fs::read(path)?;
    Ok(extract_core_from_out_str(&String::from_utf8_lossy(&bytes)))
}

pub fn write_out_core(core: &Core) -> String {
    let mut base_pairs = core.base_pairs.clone();
    base_pairs.sort_by(|a, b| {
        base_pair_sort_key(a).cmp(&base_pair_sort_key(b))
    });

    let mut multiplets = core.multiplets.clone();
    multiplets.sort_by(|a, b| (&a.indices, a.text.as_str()).cmp(&(&b.indices, b.text.as_str())));

    let mut out = String::new();
    out.push_str("BEGIN_base-pair\n");
    for bp in base_pairs {
        let line = if bp.kind == "unknown" {
            bp.text.unwrap_or_default().trim().to_string()
        } else if bp.kind == "stacked" {
            let mut tokens: Vec<String> = Vec::new();
            if let Some(syn) = bp.syn {
                if syn > 0 {
                    tokens.extend(std::iter::repeat("syn".to_string()).take(syn as usize));
                }
            }
            tokens.push("stacked".to_string());
            let rest = tokens.join(" ");
            format!(
                "{}_{}, {}: {} {}-{} {} {}: {}",
                bp.i, bp.j, bp.chain_i, bp.resseq_i, bp.base_i, bp.base_j, bp.resseq_j, bp.chain_j, rest
            )
        } else {
            let mut tokens: Vec<String> = Vec::new();
            if let Some(lw) = bp.lw {
                tokens.push(lw);
            }
            if let Some(ori) = bp.orientation {
                tokens.push(ori);
            }
            if let Some(syn) = bp.syn {
                if syn > 0 {
                    tokens.extend(std::iter::repeat("syn".to_string()).take(syn as usize));
                }
            }
            if let Some(note) = bp.note {
                if !note.trim().is_empty() {
                    tokens.push(note);
                }
            }
            let rest = tokens.join(" ");
            format!(
                "{}_{}, {}: {} {}-{} {} {}: {}",
                bp.i, bp.j, bp.chain_i, bp.resseq_i, bp.base_i, bp.base_j, bp.resseq_j, bp.chain_j, rest
            )
        };

        if !line.trim().is_empty() {
            out.push_str(line.trim_end());
            out.push('\n');
        }
    }
    out.push_str("END_base-pair\n\n");

    out.push_str("Summary of triplets and higher multiplets\n");
    out.push_str("BEGIN_multiplets\n");
    for m in multiplets {
        if m.indices.is_empty() {
            if !m.text.trim().is_empty() {
                out.push_str(m.text.trim_end());
                out.push('\n');
            }
            continue;
        }
        let idx = m
            .indices
            .iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join("_");
        out.push_str(&format!("{}_| {}\n", idx, m.text.trim_end()));
    }
    out.push_str("END_multiplets\n\n");

    if let (Some(total_pairs), Some(total_bases)) = (core.stats.total_pairs, core.stats.total_bases) {
        out.push_str(&format!(
            "  The total base pairs = {:3} (from {:4} bases)\n",
            total_pairs, total_bases
        ));
        if !core.stats.pair_type_counts.is_empty() {
            out.push_str("------------------------------------------------\n");
            let keys = core
                .stats
                .pair_type_counts
                .keys()
                .map(|k| k.as_str())
                .collect::<Vec<_>>();
            out.push_str(&format!("{}\n", keys.join(" ")));
            let vals = keys
                .iter()
                .map(|k| core.stats.pair_type_counts.get(*k).copied().unwrap_or(0))
                .map(|v| v.to_string())
                .collect::<Vec<_>>();
            out.push_str(&format!("{}\n", vals.join(" ")));
            out.push_str("------------------------------------------------\n");
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn repo_root() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("repo root")
            .to_path_buf()
    }

    fn structure_path_for_out(
        repo: &std::path::Path,
        out_path: &std::path::Path,
    ) -> std::path::PathBuf {
        let text = std::fs::read_to_string(out_path).expect("read .out text");
        let first_line = text.lines().next().unwrap_or_default();
        let prefix = "PDB data file name:";
        let input = first_line
            .strip_prefix(prefix)
            .unwrap_or_else(|| panic!("missing header '{prefix}' in {out_path:?}"))
            .trim();

        let bases = [
            repo.join(input),
            out_path
                .parent()
                .unwrap_or_else(|| std::path::Path::new("."))
                .join(input),
        ];

        for candidate in bases {
            if candidate.exists() {
                return candidate;
            }

            if let Some(file_name) = candidate.file_name().and_then(|n| n.to_str()) {
                if let Some(stripped) = file_name.strip_suffix("_new") {
                    let alt = candidate.with_file_name(stripped);
                    if alt.exists() {
                        return alt;
                    }
                }

                if let Some(stripped) = file_name.strip_suffix("_new_sort") {
                    let alt = candidate.with_file_name(stripped);
                    if alt.exists() {
                        return alt;
                    }
                }
            }
        }

        repo.join(input)
    }

    #[test]
    fn extract_core_matches_frozen_golden() {
        let repo = repo_root();
        let manifest_path = repo.join("test/golden_core/manifest.json");
        let manifest: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&manifest_path).expect("read manifest")).expect("json");

        let entries = manifest
            .get("entries")
            .and_then(|v| v.as_array())
            .expect("manifest entries");
        assert!(!entries.is_empty(), "manifest entries empty");

        for entry in entries {
            let out_rel = entry.get("out").and_then(|v| v.as_str()).expect("out");
            let core_rel = entry.get("core_json").and_then(|v| v.as_str()).expect("core_json");

            let out_path = repo.join(out_rel);
            let core_path = repo.join(core_rel);

            let parsed = extract_core_from_out_path(&out_path).expect("parse .out");
            let golden: Core = serde_json::from_str(&std::fs::read_to_string(&core_path).expect("read core"))
                .expect("json core");

            assert_eq!(parsed, golden, "core mismatch for {out_rel}");
        }
    }

    #[test]
    fn writer_roundtrip_preserves_core() {
        let repo = repo_root();
        let core_path = repo.join("test/golden_core/pdb/tr0001/tr0001.pdb.core.json");
        let golden: Core = serde_json::from_str(&std::fs::read_to_string(&core_path).expect("read core")).expect("json core");
        let out_text = write_out_core(&golden);
        let parsed = extract_core_from_out_str(&out_text);
        assert_eq!(parsed, golden);
    }

    #[test]
    fn structure_parsing_matches_legacy_out_indices() {
        let repo = repo_root();
        let manifest_path = repo.join("test/golden_core/manifest.json");
        let manifest: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&manifest_path).expect("read manifest"))
                .expect("json");

        let entries = manifest
            .get("entries")
            .and_then(|v| v.as_array())
            .expect("manifest entries");
        assert!(!entries.is_empty(), "manifest entries empty");

        let mut validated = 0usize;
        let mut skipped = 0usize;

        for entry in entries {
            let out_rel = entry.get("out").and_then(|v| v.as_str()).expect("out");
            let out_path = repo.join(out_rel);
            let structure_path = structure_path_for_out(&repo, &out_path);
            if !structure_path.exists() {
                skipped += 1;
                continue;
            }
            validated += 1;

            let bases = parse_structure_bases(&structure_path).expect("parse structure bases");
            let core = extract_core_from_out_path(&out_path).expect("parse .out core");
            let total_bases = core.stats.total_bases.expect("stats.total_bases");
            if bases.len() != total_bases as usize {
                let unfiltered = super::structure::parse_structure_nucleic_residues(&structure_path)
                    .map(|v| v.len())
                    .unwrap_or(0);
                let tail = bases
                    .iter()
                    .rev()
                    .take(5)
                    .map(|b| format!("{}:{} {}", b.chain, b.resseq, b.base))
                    .collect::<Vec<_>>();
                let has_b107 = bases.iter().any(|b| b.chain == 'B' && b.resseq == 107);
                panic!(
                    "base count mismatch for {out_rel}; got {} want {}; unfiltered={}; has B:107? {}; tail={:?}",
                    bases.len(),
                    total_bases,
                    unfiltered,
                    has_b107,
                    tail
                );
            }

            for bp in &core.base_pairs {
                if bp.i <= 0 || bp.j <= 0 {
                    continue;
                }
                let i = (bp.i - 1) as usize;
                let j = (bp.j - 1) as usize;
                assert!(i < bases.len(), "i index out of range for {out_rel}");
                assert!(j < bases.len(), "j index out of range for {out_rel}");

                let bi = &bases[i];
                let bj = &bases[j];

                assert_eq!(
                    bi.chain,
                    bp.chain_i.chars().next().unwrap_or(' '),
                    "chain_i mismatch for {out_rel} (i={})",
                    bp.i
                );
                assert_eq!(
                    bj.chain,
                    bp.chain_j.chars().next().unwrap_or(' '),
                    "chain_j mismatch for {out_rel} (j={})",
                    bp.j
                );
                assert_eq!(bi.resseq, bp.resseq_i, "resseq_i mismatch for {out_rel} (i={})", bp.i);
                assert_eq!(bj.resseq, bp.resseq_j, "resseq_j mismatch for {out_rel} (j={})", bp.j);
                assert_eq!(
                    bi.base,
                    bp.base_i.chars().next().unwrap_or(' '),
                    "base_i mismatch for {out_rel} (i={})",
                    bp.i
                );
                assert_eq!(
                    bj.base,
                    bp.base_j.chars().next().unwrap_or(' '),
                    "base_j mismatch for {out_rel} (j={})",
                    bp.j
                );
            }
        }

        assert!(validated > 0, "no golden entries validated");
        assert!(validated + skipped == entries.len(), "lost entries during iteration");
    }
}
