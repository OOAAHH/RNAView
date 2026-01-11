use crate::BasePair;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutEol {
    Lf,
    Crlf,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutPairKind {
    Pair,
    Stacked,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutBasePairLine {
    pub i: i32,
    pub j: i32,
    pub chain_i: char,
    pub resseq_i: i32,
    pub base_i: char,
    pub base_j: char,
    pub resseq_j: i32,
    pub chain_j: char,
    pub kind: OutPairKind,
    pub type_field: Option<String>, // e.g. "+/+ cis" (width=7)
    pub syn_i: bool,
    pub syn_j: bool,
    pub note: Option<String>, // Saenger roman numeral / n/a / legacy notes starting with '!'
}

fn parse_work_num(field: &str) -> Result<(i32, i32), String> {
    let s = field.trim();
    let (a, b) = s
        .split_once('_')
        .ok_or_else(|| format!("invalid work_num: {s:?}"))?;
    let i = a.parse::<i32>().map_err(|_| format!("invalid i in work_num: {s:?}"))?;
    let j = b.parse::<i32>().map_err(|_| format!("invalid j in work_num: {s:?}"))?;
    Ok((i, j))
}

fn parse_fixed_fields(line: &str) -> Result<(String, char, i32, char, char, i32, char), String> {
    if line.len() < 33 {
        return Err(format!("line too short: {line:?}"));
    }
    let work_num = line
        .get(0..9)
        .ok_or_else(|| format!("missing work_num: {line:?}"))?
        .to_string();
    let chain_i = line.chars().nth(11).ok_or_else(|| format!("missing chain_i: {line:?}"))?;
    let resseq_i = line
        .get(14..19)
        .ok_or_else(|| format!("missing resseq_i: {line:?}"))?
        .trim()
        .parse::<i32>()
        .map_err(|_| format!("invalid resseq_i: {line:?}"))?;
    let base_i = line.chars().nth(20).ok_or_else(|| format!("missing base_i: {line:?}"))?;
    let base_j = line.chars().nth(22).ok_or_else(|| format!("missing base_j: {line:?}"))?;
    let resseq_j = line
        .get(24..29)
        .ok_or_else(|| format!("missing resseq_j: {line:?}"))?
        .trim()
        .parse::<i32>()
        .map_err(|_| format!("invalid resseq_j: {line:?}"))?;
    let chain_j = line.chars().nth(30).ok_or_else(|| format!("missing chain_j: {line:?}"))?;
    Ok((work_num, chain_i, resseq_i, base_i, base_j, resseq_j, chain_j))
}

fn parse_syn_flags(mut rest: &str) -> Result<(bool, bool, &str), String> {
    let syn_i = rest.starts_with("syn ");
    if syn_i {
        rest = rest.get(4..).ok_or_else(|| "truncated syn_i".to_string())?;
    } else if rest.starts_with("  ") {
        rest = rest.get(2..).ok_or_else(|| "truncated syn_i spaces".to_string())?;
    } else {
        return Err(format!("unexpected syn_i field: {rest:?}"));
    }

    let syn_j = rest.starts_with("syn ");
    if syn_j {
        rest = rest.get(4..).ok_or_else(|| "truncated syn_j".to_string())?;
    } else if rest.starts_with("  ") {
        rest = rest.get(2..).ok_or_else(|| "truncated syn_j spaces".to_string())?;
    } else {
        return Err(format!("unexpected syn_j field: {rest:?}"));
    }

    Ok((syn_i, syn_j, rest))
}

fn parse_syn_flags_suffix(s: &str) -> Result<(bool, bool, &str), String> {
    if let Some(prefix) = s.strip_suffix("syn syn ") {
        return Ok((true, true, prefix));
    }
    if let Some(prefix) = s.strip_suffix("syn   ") {
        return Ok((true, false, prefix));
    }
    if let Some(prefix) = s.strip_suffix("  syn ") {
        return Ok((false, true, prefix));
    }
    if let Some(prefix) = s.strip_suffix("    ") {
        return Ok((false, false, prefix));
    }
    Err(format!("unexpected syn fields at end: {s:?}"))
}

pub fn parse_out_base_pair_line(line: &str) -> Result<OutBasePairLine, String> {
    let (work_num_field, chain_i, resseq_i, base_i, base_j, resseq_j, chain_j) =
        parse_fixed_fields(line)?;
    let (i, j) = parse_work_num(&work_num_field)?;

    let is_stacked = line.trim_end().ends_with("stacked");
    if is_stacked {
        let rest = line.get(33..).ok_or_else(|| format!("missing rest: {line:?}"))?;
        let rest2 = rest
            .strip_suffix("stacked")
            .ok_or_else(|| format!("expected stacked suffix: {line:?}"))?;
        let syn_part = rest2
            .strip_suffix(' ')
            .ok_or_else(|| format!("expected space before stacked: {line:?}"))?;
        let (syn_i, syn_j, tail) = parse_syn_flags(syn_part)?;
        if !tail.is_empty() {
            return Err(format!("unexpected tail before stacked: {tail:?}"));
        }
        return Ok(OutBasePairLine {
            i,
            j,
            chain_i,
            resseq_i,
            base_i,
            base_j,
            resseq_j,
            chain_j,
            kind: OutPairKind::Stacked,
            type_field: None,
            syn_i,
            syn_j,
            note: None,
        });
    }

    let rest = line.get(33..).ok_or_else(|| format!("missing rest: {line:?}"))?;
    let (before_note, note) = rest
        .rsplit_once(' ')
        .ok_or_else(|| format!("missing space before note: {line:?}"))?;
    let (syn_i, syn_j, before_syn) = parse_syn_flags_suffix(before_note)?;
    let type_field = before_syn
        .strip_suffix("   ")
        .ok_or_else(|| format!("expected 3-space delimiter before syn fields: {line:?}"))?
        .to_string();
    let note = if note.is_empty() { None } else { Some(note.to_string()) };

    Ok(OutBasePairLine {
        i,
        j,
        chain_i,
        resseq_i,
        base_i,
        base_j,
        resseq_j,
        chain_j,
        kind: OutPairKind::Pair,
        type_field: Some(type_field),
        syn_i,
        syn_j,
        note,
    })
}

pub fn format_out_base_pair_line(bp: &OutBasePairLine) -> Result<String, String> {
    let syn_i = if bp.syn_i { "syn " } else { "  " };
    let syn_j = if bp.syn_j { "syn " } else { "  " };
    let work_num = format!("{}_{}", bp.i, bp.j);

    match bp.kind {
        OutPairKind::Stacked => Ok(format!(
            "{work_num:>9}, {chain_i}: {resseq_i:>5} {base_i}-{base_j} {resseq_j:>5} {chain_j}: {syn_i}{syn_j} stacked",
            chain_i = bp.chain_i,
            resseq_i = bp.resseq_i,
            base_i = bp.base_i,
            base_j = bp.base_j,
            resseq_j = bp.resseq_j,
            chain_j = bp.chain_j,
            syn_i = syn_i,
            syn_j = syn_j,
        )),
        OutPairKind::Pair => {
            let type_field = bp
                .type_field
                .as_deref()
                .ok_or_else(|| "missing type_field for Pair".to_string())?;
            let note = bp.note.as_deref().unwrap_or("");
            Ok(format!(
                "{work_num:>9}, {chain_i}: {resseq_i:>5} {base_i}-{base_j} {resseq_j:>5} {chain_j}: {type_field:>7}   {syn_i}{syn_j} {note}",
                chain_i = bp.chain_i,
                resseq_i = bp.resseq_i,
                base_i = bp.base_i,
                base_j = bp.base_j,
                resseq_j = bp.resseq_j,
                chain_j = bp.chain_j,
                type_field = type_field,
                syn_i = syn_i,
                syn_j = syn_j,
                note = note,
            ))
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct OutFull {
    pub eol: OutEol,
    pub trailing_newline: bool,
    pub pdb_data_file_name: String,
    pub uncommon_lines: Vec<String>,
    pub bprs: Option<[f64; 6]>,
    pub base_pairs: Vec<OutBasePairLine>,
    pub multiplets: Option<Vec<String>>,
    pub blank_lines_after_end_base_pair: usize,
    pub blank_lines_after_end_multiplets: usize,
    pub total_base_pairs: i64,
    pub total_bases: i64,
    pub type_counts_1_to_7: [i64; 7],
    pub type_counts_8_to_13: [i64; 6],
}

fn parse_bprs(lines: &[String], start_idx: usize) -> Result<([f64; 6], usize), String> {
    let mut bprs = [0.0f64; 6];
    let mut idx = start_idx;
    for i in 0..6 {
        let line = lines.get(idx).ok_or_else(|| "unexpected EOF reading BPRS".to_string())?;
        let num = line
            .split_whitespace()
            .next()
            .ok_or_else(|| format!("missing bprs value: {line:?}"))?
            .parse::<f64>()
            .map_err(|_| format!("invalid bprs value: {line:?}"))?;
        bprs[i] = num;
        idx += 1;
    }
    Ok((bprs, idx))
}

fn parse_total_line(line: &str) -> Result<(i64, i64), String> {
    // Matches: "  The total base pairs =%4ld (from %4ld bases)"
    let prefix = "  The total base pairs =";
    if !line.starts_with(prefix) {
        return Err(format!("unexpected total line: {line:?}"));
    }
    let rest = &line[prefix.len()..];
    let (pairs_part, rest2) = rest
        .split_once(" (from")
        .ok_or_else(|| format!("unexpected total line: {line:?}"))?;
    let pairs = pairs_part.trim().parse::<i64>().map_err(|_| format!("invalid pair count: {line:?}"))?;
    let rest2 = rest2
        .strip_suffix(" bases)")
        .ok_or_else(|| format!("unexpected total line: {line:?}"))?;
    let bases = rest2.trim().parse::<i64>().map_err(|_| format!("invalid base count: {line:?}"))?;
    Ok((pairs, bases))
}

fn parse_counts(line: &str, want: usize) -> Result<Vec<i64>, String> {
    let nums: Vec<i64> = line
        .split_whitespace()
        .filter(|t| !t.is_empty())
        .map(|t| t.parse::<i64>().map_err(|_| format!("invalid count token: {t:?}")))
        .collect::<Result<Vec<_>, _>>()?;
    if nums.len() != want {
        return Err(format!("expected {want} counts, got {}: {line:?}", nums.len()));
    }
    Ok(nums)
}

pub fn parse_out_full(text: &str) -> Result<OutFull, String> {
    let eol = if text.contains('\r') {
        OutEol::Crlf
    } else {
        OutEol::Lf
    };
    let trailing_newline = match eol {
        OutEol::Lf => text.ends_with('\n'),
        OutEol::Crlf => text.ends_with("\r\n"),
    };

    let lines: Vec<String> = text.lines().map(|s| s.to_string()).collect();
    if lines.is_empty() {
        return Err("empty .out".to_string());
    }

    let pdb_data_file_name = lines[0].to_string();
    let mut idx = 1usize;

    let mut uncommon_lines: Vec<String> = Vec::new();
    while idx < lines.len() && lines[idx].starts_with("uncommon ") {
        uncommon_lines.push(lines[idx].to_string());
        idx += 1;
    }

    let mut bprs: Option<[f64; 6]> = None;
    if idx < lines.len() && lines[idx] != "BEGIN_base-pair" {
        if idx >= lines.len() || lines[idx] != "-----------------------------------------------------------" {
            return Err(format!("expected preamble separator at line {}", idx + 1));
        }
        idx += 1;

        // Criteria header
        if idx >= lines.len() || lines[idx] != "CRITERIA USED TO GENERATE BASE-PAIR: " {
            return Err(format!("expected criteria header at line {}", idx + 1));
        }
        idx += 1;

        let (bprs_parsed, idx2) = parse_bprs(&lines, idx)?;
        bprs = Some(bprs_parsed);
        idx = idx2;

        if idx >= lines.len() || lines[idx] != "-----------------------------------------------------------" {
            return Err(format!("expected separator after BPRS at line {}", idx + 1));
        }
        idx += 1;

        // Skip instructions until we hit BEGIN_base-pair.
        while idx < lines.len() && lines[idx] != "BEGIN_base-pair" {
            idx += 1;
        }
    }
    if idx >= lines.len() || lines[idx] != "BEGIN_base-pair" {
        return Err("missing BEGIN_base-pair".to_string());
    }
    idx += 1;

    let mut base_pairs: Vec<OutBasePairLine> = Vec::new();
    while idx < lines.len() && lines[idx] != "END_base-pair" {
        let line = lines[idx].to_string();
        if !line.trim().is_empty() {
            base_pairs.push(parse_out_base_pair_line(&line)?);
        }
        idx += 1;
    }
    if idx >= lines.len() || lines[idx] != "END_base-pair" {
        return Err("missing END_base-pair".to_string());
    }
    idx += 1;

    let mut multiplets: Option<Vec<String>> = None;
    let mut blank_lines_after_end_base_pair = 0usize;
    while idx < lines.len() && lines[idx].is_empty() {
        blank_lines_after_end_base_pair += 1;
        idx += 1;
    }
    let mut blank_lines_after_end_multiplets = 0usize;
    if idx < lines.len() && lines[idx] == "Summary of triplets and higher multiplets" {
        idx += 1;
        if idx >= lines.len() || lines[idx] != "BEGIN_multiplets" {
            return Err("missing BEGIN_multiplets".to_string());
        }
        idx += 1;
        let mut out: Vec<String> = Vec::new();
        while idx < lines.len() && lines[idx] != "END_multiplets" {
            out.push(lines[idx].to_string());
            idx += 1;
        }
        if idx >= lines.len() || lines[idx] != "END_multiplets" {
            return Err("missing END_multiplets".to_string());
        }
        idx += 1;
        while idx < lines.len() && lines[idx].is_empty() {
            blank_lines_after_end_multiplets += 1;
            idx += 1;
        }
        multiplets = Some(out);
    }

    let total_line = lines.get(idx).ok_or_else(|| "missing total line".to_string())?;
    let (total_base_pairs, total_bases) = parse_total_line(total_line)?;
    idx += 1;

    if idx >= lines.len() || lines[idx] != "------------------------------------------------" {
        return Err("missing stats separator".to_string());
    }
    idx += 1;
    if idx >= lines.len() || lines[idx] != " Standard  WW--cis  WW-tran  HH--cis  HH-tran  SS--cis  SS-tran" {
        return Err("missing stats header 1".to_string());
    }
    idx += 1;
    let counts1 = lines.get(idx).ok_or_else(|| "missing stats counts 1".to_string())?;
    let c1 = parse_counts(counts1, 7)?;
    idx += 1;
    if idx >= lines.len() || lines[idx] != "  WH--cis  WH-tran  WS--cis  WS-tran  HS--cis  HS-tran" {
        return Err("missing stats header 2".to_string());
    }
    idx += 1;
    let counts2 = lines.get(idx).ok_or_else(|| "missing stats counts 2".to_string())?;
    let c2 = parse_counts(counts2, 6)?;
    idx += 1;

    if idx >= lines.len() || lines[idx] != "------------------------------------------------" {
        return Err("missing final separator".to_string());
    }

    let mut type_counts_1_to_7 = [0i64; 7];
    for (i, v) in c1.into_iter().enumerate() {
        type_counts_1_to_7[i] = v;
    }
    let mut type_counts_8_to_13 = [0i64; 6];
    for (i, v) in c2.into_iter().enumerate() {
        type_counts_8_to_13[i] = v;
    }

    Ok(OutFull {
        eol,
        trailing_newline,
        pdb_data_file_name,
        uncommon_lines,
        bprs,
        base_pairs,
        multiplets,
        blank_lines_after_end_base_pair,
        blank_lines_after_end_multiplets,
        total_base_pairs,
        total_bases,
        type_counts_1_to_7,
        type_counts_8_to_13,
    })
}

pub fn write_out_full(out: &OutFull) -> Result<String, String> {
    let mut s = String::new();

    s.push_str(&out.pdb_data_file_name);
    s.push('\n');

    for line in &out.uncommon_lines {
        s.push_str(line);
        s.push('\n');
    }

    if let Some(bprs) = &out.bprs {
        s.push_str("-----------------------------------------------------------\n");
        s.push_str("CRITERIA USED TO GENERATE BASE-PAIR: \n");
        s.push_str(&format!(
            "{:6.2} --> upper H-bond length limits (ON..ON).\n",
            bprs[0]
        ));
        s.push_str(&format!(
            "{:6.2} --> max. distance between paired base origins.\n",
            bprs[1]
        ));
        s.push_str(&format!(
            "{:6.2} --> max. vertical distance between paired base origins.\n",
            bprs[2]
        ));
        s.push_str(&format!(
            "{:6.2} --> max. angle between paired bases [0-90].\n",
            bprs[3]
        ));
        s.push_str(&format!(
            "{:6.2} --> min. distance between RN9/YN1 atoms.\n",
            bprs[4]
        ));
        s.push_str(&format!(
            "{:6.2} --> max. distance criterion for helix break[0-12]\n",
            bprs[5]
        ));
        s.push_str("-----------------------------------------------------------\n");
        s.push_str("BASE-PAIR INSTRUCTIONS: \n");
        s.push_str("Column 1 is rnaview assigned base numbers n1_n2, start from 1.\n");
        s.push_str("Column 2 & 3 are chain ID & residue number in input PDB file.\n");
        s.push_str("Column 4 is for base pair. The left & right are the bases as \n");
        s.push_str("         identified by column 2 & 3 and 5 & 6.\n");
        s.push_str("Column 5 & 6 are residue number & chain ID in input PDB file.\n");
        s.push_str("Column 7 is for base pair annotation. The standard Watson-Crick\n");
        s.push_str("         (W.C.) pairs are annotated as -/- (AU,AT) or +/+ (GC).\n");
        s.push_str("         Other pairs are annotated as Leontis_Westhof Classification.\n");
        s.push_str("         The three edges: W(Watson-Crick); H(Hoogsteen); S(suger).\n");
        s.push_str("         e.g. W/H means the pair is edge of Watson-Crick & Hoogsteen.\n");
        s.push_str("Column 8 is glycosidic bond orientation (either cis or trans).\n");
        s.push_str("         e.g. 'W/H cis' means the pair is interaction on Watson-Crick\n");
        s.push_str("         and Hoogsteen side, glycosidic bond orientation is 'cis'.\n");
        s.push_str("Column 9 corresponds to Saenger Classification.\n\n");
        s.push_str("Other columns: \n");
        s.push_str("        Syn sugar-base conformations are annotated as (syn).\n");
        s.push_str("        Stacked base pairs are annotated as (stack).\n");
        s.push_str("        Non-identified edges are annotated as (.) or (?)\n");
        s.push_str("        Tertiary interactions are marked by (!) in the line.\n");
        s.push_str("Reference:\n");
        s.push_str("Yang et al (2003) Nucleic Acids Research, Vol31,No13,p3450-3461.\n");
        s.push_str("-----------------------------------------------------------\n");
    }

    s.push_str("BEGIN_base-pair\n");
    for bp in &out.base_pairs {
        let line = format_out_base_pair_line(bp)?;
        s.push_str(&line);
        s.push('\n');
    }
    s.push_str("END_base-pair\n");

    for _ in 0..out.blank_lines_after_end_base_pair {
        s.push('\n');
    }

    if let Some(lines) = &out.multiplets {
        s.push_str("Summary of triplets and higher multiplets\n");
        s.push_str("BEGIN_multiplets\n");
        for line in lines {
            s.push_str(line);
            s.push('\n');
        }
        s.push_str("END_multiplets\n");
        for _ in 0..out.blank_lines_after_end_multiplets {
            s.push('\n');
        }
    }

    s.push_str(&format!(
        "  The total base pairs ={total_pairs:4} (from {total_bases:4} bases)\n",
        total_pairs = out.total_base_pairs,
        total_bases = out.total_bases
    ));
    s.push_str("------------------------------------------------\n");
    s.push_str(" Standard  WW--cis  WW-tran  HH--cis  HH-tran  SS--cis  SS-tran\n");
    for v in out.type_counts_1_to_7 {
        s.push_str(&format!("{v:9}"));
    }
    s.push('\n');
    s.push_str("  WH--cis  WH-tran  WS--cis  WS-tran  HS--cis  HS-tran\n");
    for v in out.type_counts_8_to_13 {
        s.push_str(&format!("{v:9}"));
    }
    s.push('\n');
    s.push_str("------------------------------------------------\n");

    if !out.trailing_newline && s.ends_with('\n') {
        s.pop();
    }
    if out.eol == OutEol::Crlf {
        s = s.replace('\n', "\r\n");
    }
    Ok(s)
}

impl From<OutBasePairLine> for BasePair {
    fn from(v: OutBasePairLine) -> Self {
        BasePair {
            i: v.i,
            j: v.j,
            chain_i: v.chain_i.to_string(),
            resseq_i: v.resseq_i,
            base_i: v.base_i.to_string(),
            base_j: v.base_j.to_string(),
            resseq_j: v.resseq_j,
            chain_j: v.chain_j.to_string(),
            kind: match v.kind {
                OutPairKind::Pair => "pair".to_string(),
                OutPairKind::Stacked => "stacked".to_string(),
            },
            lw: v.type_field.as_ref().map(|s| s.split_whitespace().next().unwrap_or("").to_string()),
            orientation: None,
            syn: Some((v.syn_i as i32) + (v.syn_j as i32)),
            note: v.note,
            text: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn repo_root() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).parent().unwrap().to_path_buf()
    }

    #[test]
    fn base_pair_line_roundtrip_matches_golden_text() {
        let sample = "     1_72, A:     1 G-C    72 A: +/+ cis         XIX";
        let parsed = parse_out_base_pair_line(sample).expect("parse");
        let rendered = format_out_base_pair_line(&parsed).expect("format");
        assert_eq!(rendered, sample);
    }

    #[test]
    fn out_full_roundtrip_is_byte_exact_for_manifest_entries() {
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
            let out_path = repo.join(out_rel);
            let bytes = std::fs::read(&out_path).expect("read .out");
            let text = String::from_utf8_lossy(&bytes);
            let parsed = parse_out_full(&text).expect("parse full");
            let rendered = write_out_full(&parsed).expect("write");
            if rendered.as_bytes() != bytes {
                panic!("byte mismatch for {out_rel}");
            }
        }
    }
}
