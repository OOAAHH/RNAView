use rnaview_hotcore::{
    compute_out_full_from_structure_with_policies, extract_core_from_out_path_with_out_index,
    extract_core_from_out_str_with_out_index, parse_structure_bases, render_vrml_from_pairs_json,
    write_out_core, write_out_full, BaseResidue, ChainIdPolicy, Core, HydrogenPolicy,
    MissingInsertionCodePolicy, PairsJson, Semantics, SemanticsConfig, Source,
};
use pdbtbx::{
    ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain, Element,
    Format, ReadOptions, StrictnessLevel,
};
use std::io::BufRead;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Stdio;

fn usage() -> ! {
    eprintln!(
        "Usage:\n  rnaview-hotcore from-out <file.out> [-o pairs.json]\n  rnaview-hotcore from-structure <file.pdb|file.cif>\n      [--format pdb|cif]\n      [--oracle legacy|out|compute]\n      [--mmcif-parser legacy|pdbtbx]\n      [--semantics legacy-v1|science-v1]\n      [--hydrogen-policy legacy-mmcif-bug|discard-all|keep-all]\n      [--missing-insertion-code legacy-question-mark|none]\n      [--chain-id-policy legacy-1char|unique-1char]\n      [-o pairs.json]\n      [--emit-out file.out]\n  rnaview-hotcore write-out <pairs.json> [-o candidate.out]\n  rnaview-hotcore render-wrl <pairs.json> [--source <file.pdb|file.cif>] [-o file.wrl]\n  rnaview-hotcore render-2d <pairs.json> [--source <file.pdb|file.cif>] [--out-xml file.xml] [--out-ps file.ps]"
    );
    std::process::exit(2);
}

fn infer_format(path: &Path) -> Option<String> {
    let ext = path.extension()?.to_string_lossy().to_ascii_lowercase();
    match ext.as_str() {
        "pdb" | "ent" => Some("pdb".to_string()),
        "cif" => Some("cif".to_string()),
        _ => None,
    }
}

fn default_rnaview_root(legacy_bin: &Path) -> Option<PathBuf> {
    let bin_dir = legacy_bin.parent()?;
    if bin_dir.file_name().is_some_and(|n| n == "bin") {
        return Some(bin_dir.parent()?.to_path_buf());
    }
    None
}

fn default_legacy_bin(rnaview_root: &Path) -> PathBuf {
    rnaview_root.join("bin").join("rnaview")
}

fn mmcif_resname_to_legacy(resname: &str) -> String {
    let s = resname.trim().to_ascii_uppercase();
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

fn atom_name_field_legacy(name: &str) -> String {
    let mut s = name.trim().to_ascii_uppercase();
    if s.len() > 4 {
        s = s.chars().take(4).collect();
    }
    if s.len() == 4 {
        return s;
    }
    let mut out = String::with_capacity(4);
    out.push(' ');
    out.push_str(&s);
    while out.len() < 4 {
        out.push(' ');
    }
    out
}

fn residue_name_field_legacy(name: &str) -> String {
    let mut s = mmcif_resname_to_legacy(name);
    if s.len() > 3 {
        s = s.chars().take(3).collect();
    }
    format!("{s:>3}")
}

fn pdb_charge_field(charge: &str) -> String {
    let s = charge.trim();
    if s.len() >= 2 {
        return s.chars().rev().take(2).collect::<Vec<_>>().into_iter().rev().collect();
    }
    format!("{s:>2}")
}

fn mmcif_simple_value_after_tag(line: &str, tag: &str) -> Option<String> {
    let s = line.trim_start();
    if !s.starts_with(tag) {
        return None;
    }
    let rest = s[tag.len()..].trim_start();
    if rest.is_empty() {
        return Some(String::new());
    }

    let first = rest.chars().next().unwrap_or(' ');
    if first == '\'' || first == '"' {
        let quote = first;
        let mut chars = rest.chars();
        let _ = chars.next();
        let mut out = String::new();
        for ch in chars {
            if ch == quote {
                break;
            }
            out.push(ch);
        }
        return Some(out);
    }

    let token = rest.split_whitespace().next().unwrap_or("").to_string();
    Some(token)
}

fn mmcif_best_model_if_nmr(path: &Path) -> std::io::Result<Option<usize>> {
    let f = std::fs::File::open(path)?;
    let r = std::io::BufReader::new(f);

    let mut is_nmr: Option<bool> = None;
    let mut conformer_id: Option<String> = None;
    let mut representative: Option<String> = None;

    for (idx, line) in r.lines().enumerate() {
        if idx > 20000 {
            break;
        }
        let line = line?;

        if is_nmr.is_none() {
            if let Some(v) = mmcif_simple_value_after_tag(&line, "_exptl.method") {
                is_nmr = Some(v.to_ascii_uppercase().contains("NMR"));
            }
        }
        if conformer_id.is_none() {
            if let Some(v) = mmcif_simple_value_after_tag(&line, "_pdbx_nmr_representative.conformer_id") {
                conformer_id = Some(v);
            }
        }
        if representative.is_none() {
            if let Some(v) =
                mmcif_simple_value_after_tag(&line, "_pdbx_nmr_ensemble.representative_conformer")
            {
                representative = Some(v);
            }
        }

        if line.trim_start().starts_with("_atom_site.") {
            break;
        }
    }

    if is_nmr != Some(true) {
        return Ok(None);
    }

    let pick = |v: Option<String>| -> Option<String> {
        let s = v?.trim().to_string();
        if s.is_empty() || s == "?" || s == "." {
            None
        } else {
            Some(s)
        }
    };

    let best = pick(conformer_id)
        .or_else(|| pick(representative))
        .unwrap_or_else(|| "1".to_string());
    Ok(Some(best.parse::<usize>().unwrap_or(1)))
}

fn write_legacy_pdb_from_pdbtbx(
    pdb: &pdbtbx::PDB,
    out_path: &Path,
) -> std::io::Result<()> {
    struct Row {
        serial: usize,
        hetero: bool,
        atom_name: String,
        alt_loc: char,
        resname: String,
        chain: char,
        resseq: isize,
        icode: char,
        x: f64,
        y: f64,
        z: f64,
        occupancy: f64,
        b_factor: f64,
        element: String,
        charge: String,
    }

    let mut rows: Vec<Row> = Vec::new();
    for h in pdb.atoms_with_hierarchy() {
        let atom = h.atom();
        let chain_id = h.chain().id();
        let chain = chain_id.chars().next().unwrap_or(' ');
        let residue = h.residue();
        let resseq = residue.serial_number();
        let icode = residue
            .insertion_code()
            .and_then(|s| s.chars().next())
            .unwrap_or(' ');
        let conformer = h.conformer();
        let resname = residue_name_field_legacy(conformer.name());
        let alt_loc = conformer
            .alternative_location()
            .and_then(|s| s.chars().next())
            .unwrap_or(' ');
        let atom_name = atom_name_field_legacy(atom.name());
        let element = atom
            .element()
            .map(Element::symbol)
            .unwrap_or("")
            .to_string();
        let charge = atom.pdb_charge();

        rows.push(Row {
            serial: atom.serial_number(),
            hetero: atom.hetero(),
            atom_name,
            alt_loc,
            resname,
            chain,
            resseq,
            icode,
            x: atom.x(),
            y: atom.y(),
            z: atom.z(),
            occupancy: atom.occupancy(),
            b_factor: atom.b_factor(),
            element,
            charge,
        });
    }
    rows.sort_by_key(|r| r.serial);

    let file = std::fs::File::create(out_path)?;
    let mut w = std::io::BufWriter::new(file);

    for row in rows {
        let record = if row.hetero { "HETATM" } else { "ATOM  " };
        let serial = row.serial % 100000;
        let element_field = format!("{:>2}", row.element);
        let charge_field = pdb_charge_field(&row.charge);

        writeln!(
            w,
            "{record}{serial:>5} {atom_name}{alt_loc}{resname} {chain}{resseq:>4}{icode}   {x:>8.3}{y:>8.3}{z:>8.3}{occupancy:>6.2}{b_factor:>6.2}          {element_field}{charge_field}",
            atom_name = row.atom_name,
            alt_loc = row.alt_loc,
            resname = row.resname,
            chain = row.chain,
            resseq = row.resseq,
            icode = row.icode,
            x = row.x,
            y = row.y,
            z = row.z,
            occupancy = row.occupancy,
            b_factor = row.b_factor,
        )?;
    }
    writeln!(w, "END")?;
    w.flush()
}

struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn create(prefix: &str) -> std::io::Result<Self> {
        let base = std::env::temp_dir();
        let pid = std::process::id();
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        for i in 0..100u32 {
            let cand = base.join(format!("{prefix}-{pid}-{now}-{i}"));
            match std::fs::create_dir(&cand) {
                Ok(()) => return Ok(Self { path: cand }),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e),
            }
        }
        Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            "failed to create temp dir",
        ))
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        usage();
    }

    let cmd = args.remove(0);
    match cmd.as_str() {
        "from-out" => {
            if args.is_empty() {
                usage();
            }
            let input = PathBuf::from(args.remove(0));
            let mut output: Option<PathBuf> = None;
            while !args.is_empty() {
                let flag = args.remove(0);
                if flag == "-o" || flag == "--output" {
                    if args.is_empty() {
                        usage();
                    }
                    output = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                usage();
            }

            let core: Core = extract_core_from_out_path_with_out_index(&input)?;
            let pairs = PairsJson {
                schema_version: 1,
                source: Some(Source {
                    path: input.to_string_lossy().to_string(),
                    format: "out".to_string(),
                    id_scheme: None,
                    model: None,
                }),
                options: Some(serde_json::json!({"engine":"rust","source":"out"})),
                core,
            };

            let json_text = serde_json::to_string(&serde_json::to_value(&pairs)?)? + "\n";
            if let Some(out_path) = output {
                std::fs::write(out_path, json_text)?;
            } else {
                print!("{json_text}");
            }
            Ok(())
        }
        "from-structure" => {
            if args.is_empty() {
                usage();
            }
            let input = PathBuf::from(args.remove(0));

            let mut format: Option<String> = None;
            let mut output: Option<PathBuf> = None;
            let mut legacy_bin: Option<PathBuf> = None;
            let mut rnaview_root: Option<PathBuf> = None;
            let mut emit_out: Option<PathBuf> = None;
            let mut mmcif_parser: String = "legacy".to_string();
            let mut oracle: String = "legacy".to_string();
            let mut semantics: Option<String> = None;
            let mut hydrogen_policy: Option<String> = None;
            let mut missing_insertion_code: Option<String> = None;
            let mut chain_id_policy: Option<String> = None;

            while !args.is_empty() {
                let flag = args.remove(0);
                if flag == "-o" || flag == "--output" {
                    if args.is_empty() {
                        usage();
                    }
                    output = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                if flag == "--format" {
                    if args.is_empty() {
                        usage();
                    }
                    format = Some(args.remove(0));
                    continue;
                }
                if flag == "--mmcif-parser" {
                    if args.is_empty() {
                        usage();
                    }
                    mmcif_parser = args.remove(0);
                    continue;
                }
                if flag == "--oracle" {
                    if args.is_empty() {
                        usage();
                    }
                    oracle = args.remove(0);
                    continue;
                }
                if flag == "--semantics" {
                    if args.is_empty() {
                        usage();
                    }
                    semantics = Some(args.remove(0));
                    continue;
                }
                if flag == "--hydrogen-policy" {
                    if args.is_empty() {
                        usage();
                    }
                    hydrogen_policy = Some(args.remove(0));
                    continue;
                }
                if flag == "--missing-insertion-code" {
                    if args.is_empty() {
                        usage();
                    }
                    missing_insertion_code = Some(args.remove(0));
                    continue;
                }
                if flag == "--chain-id-policy" {
                    if args.is_empty() {
                        usage();
                    }
                    chain_id_policy = Some(args.remove(0));
                    continue;
                }
                if flag == "--legacy-bin" {
                    if args.is_empty() {
                        usage();
                    }
                    legacy_bin = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                if flag == "--rnaview-root" {
                    if args.is_empty() {
                        usage();
                    }
                    rnaview_root = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                if flag == "--emit-out" || flag == "--emit-legacy-out" {
                    if args.is_empty() {
                        usage();
                    }
                    emit_out = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                usage();
            }

            let fmt = match format {
                Some(s) => s,
                None => infer_format(&input).ok_or_else(|| {
                    std::io::Error::new(std::io::ErrorKind::InvalidInput, "unknown input format; pass --format")
                })?,
            };
            if fmt != "pdb" && fmt != "cif" {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "format must be pdb|cif",
                )));
            }

            if mmcif_parser != "legacy" && mmcif_parser != "pdbtbx" {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "mmcif-parser must be legacy|pdbtbx",
                )));
            }

            if oracle != "legacy" && oracle != "out" && oracle != "compute" {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "oracle must be legacy|out|compute",
                )));
            }
            if oracle != "compute"
                && (semantics.is_some()
                    || hydrogen_policy.is_some()
                    || missing_insertion_code.is_some()
                    || chain_id_policy.is_some())
            {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "--semantics/--hydrogen-policy/--missing-insertion-code/--chain-id-policy require --oracle compute",
                )));
            }

            if oracle == "out" {
                let oracle_out = PathBuf::from(format!("{}.out", input.to_string_lossy()));
                if !oracle_out.exists() {
                    return Err(Box::new(std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("missing oracle .out next to input: {oracle_out:?}"),
                    )));
                }

                let bases: Vec<BaseResidue> = parse_structure_bases(&input)?;
                let core: Core = extract_core_from_out_path_with_out_index(&oracle_out)?;

                let want = core
                    .stats
                    .total_bases
                    .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "missing stats.total_bases"))?
                    as usize;
                if bases.len() != want {
                    return Err(Box::new(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("base count mismatch: got {} want {}", bases.len(), want),
                    )));
                }

                for bp in &core.base_pairs {
                    if bp.i <= 0 || bp.j <= 0 {
                        continue;
                    }
                    let i = (bp.i - 1) as usize;
                    let j = (bp.j - 1) as usize;
                    if i >= bases.len() || j >= bases.len() {
                        return Err(Box::new(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!("base pair index out of range: {}_{}", bp.i, bp.j),
                        )));
                    }
                    let bi = &bases[i];
                    let bj = &bases[j];
                    if bi.chain.to_string() != bp.chain_i
                        || bj.chain.to_string() != bp.chain_j
                        || bi.resseq != bp.resseq_i
                        || bj.resseq != bp.resseq_j
                        || bi.base.to_string() != bp.base_i
                        || bj.base.to_string() != bp.base_j
                    {
                        return Err(Box::new(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!("base-pair identity mismatch at {}_{}", bp.i, bp.j),
                        )));
                    }
                }

                if let Some(dst) = emit_out {
                    std::fs::copy(&oracle_out, dst)?;
                }

                let pairs = PairsJson {
                    schema_version: 1,
                    source: Some(Source {
                        path: input.to_string_lossy().to_string(),
                        format: fmt.clone(),
                        id_scheme: if fmt == "cif" { Some("auth".to_string()) } else { None },
                        model: None,
                    }),
                    options: Some(serde_json::json!({"engine":"rust","oracle":"out"})),
                    core,
                };

                let json_text = serde_json::to_string(&serde_json::to_value(&pairs)?)? + "\n";
                if let Some(out_path) = output {
                    std::fs::write(out_path, json_text)?;
                } else {
                    print!("{json_text}");
                }
                return Ok(());
            }

            if oracle == "compute" {
                let semantics = semantics
                    .as_deref()
                    .map(Semantics::parse_cli)
                    .transpose()
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?
                    .unwrap_or(Semantics::LegacyV1);
                let mut structure_policies = SemanticsConfig::defaults(semantics).policies.structure;
                if let Some(s) = hydrogen_policy.as_deref() {
                    structure_policies.hydrogen_policy = HydrogenPolicy::parse_cli(s)
                        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;
                }
                if let Some(s) = missing_insertion_code.as_deref() {
                    structure_policies.missing_insertion_code_policy =
                        MissingInsertionCodePolicy::parse_cli(s)
                            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;
                }
                if let Some(s) = chain_id_policy.as_deref() {
                    structure_policies.chain_id_policy = ChainIdPolicy::parse_cli(s)
                        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;
                }

                let root_opt = rnaview_root
                    .clone()
                    .or_else(|| std::env::var_os("RNAVIEW").map(PathBuf::from));
                let input_abs = input.canonicalize().unwrap_or_else(|_| input.clone());
                let input_arg = if let Some(root) = &root_opt {
                    let root_abs = root.canonicalize().unwrap_or_else(|_| root.clone());
                    input_abs
                        .strip_prefix(&root_abs)
                        .ok()
                        .map(|p| p.to_string_lossy().to_string())
                        .unwrap_or_else(|| input_abs.to_string_lossy().to_string())
                } else {
                    input_abs.to_string_lossy().to_string()
                };

                let header_path = if fmt == "pdb" {
                    format!("{input_arg}_new")
                } else {
                    input_arg
                };
                let pdb_data_file_name = format!("PDB data file name: {header_path}");

                let out_full =
                    compute_out_full_from_structure_with_policies(&input, pdb_data_file_name, &structure_policies)?;
                let out_text = write_out_full(&out_full)
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

                if let Some(dst) = emit_out {
                    std::fs::write(dst, out_text.as_bytes())?;
                }

                let core: Core = extract_core_from_out_str_with_out_index(&out_text);
                let pairs = PairsJson {
                    schema_version: 1,
                    source: Some(Source {
                        path: input.to_string_lossy().to_string(),
                        format: fmt.clone(),
                        id_scheme: if fmt == "cif" { Some("auth".to_string()) } else { None },
                        model: None,
                    }),
                    options: Some(serde_json::json!({
                        "engine":"rust",
                        "oracle":"compute",
                        "semantics": semantics.as_str(),
                        "policies": { "structure": structure_policies },
                    })),
                    core,
                };

                let json_text = serde_json::to_string(&serde_json::to_value(&pairs)?)? + "\n";
                if let Some(out_path) = output {
                    std::fs::write(out_path, json_text)?;
                } else {
                    print!("{json_text}");
                }
                return Ok(());
            }

            let legacy_bin = match legacy_bin {
                Some(p) => p,
                None => {
                    let root = rnaview_root
                        .clone()
                        .or_else(|| std::env::var_os("RNAVIEW").map(PathBuf::from))
                        .ok_or_else(|| {
                            std::io::Error::new(
                                std::io::ErrorKind::InvalidInput,
                                "missing RNAVIEW env; pass --rnaview-root or --legacy-bin",
                            )
                        })?;
                    default_legacy_bin(&root)
                }
            };
            let rnaview_root = match rnaview_root
                .clone()
                .or_else(|| std::env::var_os("RNAVIEW").map(PathBuf::from))
                .or_else(|| default_rnaview_root(&legacy_bin))
            {
                Some(p) => p,
                None => {
                    return Err(Box::new(std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "missing RNAVIEW root; pass --rnaview-root or set RNAVIEW",
                    )))
                }
            };

            let work = TempDir::create("rnaview-hotcore")?;
            let input_abs = input.canonicalize().unwrap_or_else(|_| input.clone());
            let root_abs = rnaview_root.canonicalize().unwrap_or_else(|_| rnaview_root.clone());

            let input_rel = input_abs.strip_prefix(&root_abs).ok().map(|p| p.to_path_buf());
            let (input_arg, local_input) = if let Some(rel) = &input_rel {
                (rel.to_string_lossy().to_string(), work.path.join(rel))
            } else {
                let input_name = input_abs.file_name().ok_or_else(|| {
                    std::io::Error::new(std::io::ErrorKind::InvalidInput, "input has no file name")
                })?;
                (input_name.to_string_lossy().to_string(), work.path.join(input_name))
            };
            if let Some(parent) = local_input.parent() {
                std::fs::create_dir_all(parent)?;
            }
            std::fs::copy(&input_abs, &local_input)?;

            let mut legacy_flag: &str = if fmt == "pdb" { "--pdb" } else { "--cif" };
            let mut source_model: Option<u32> = None;
            if fmt == "cif" && mmcif_parser == "pdbtbx" {
                let best_model = mmcif_best_model_if_nmr(&local_input)?;
                source_model = best_model.map(|m| m as u32);
                if best_model.is_none() {
                    let mut opts = ReadOptions::new();
                    opts.set_level(StrictnessLevel::Loose)
                        .set_format(Format::Mmcif)
                        .set_discard_hydrogens(true)
                        .set_only_first_model(true)
                        .set_only_atomic_coords(true)
                        .set_capitalise_chains(false);

                    let path_str = local_input.to_string_lossy();
                    let (pdb, _warnings) = opts.read(path_str.as_ref()).map_err(|errs| {
                        let msg = errs
                            .iter()
                            .map(|e| e.to_string())
                            .collect::<Vec<_>>()
                            .join("\n");
                        std::io::Error::new(std::io::ErrorKind::InvalidData, msg)
                    })?;

                    write_legacy_pdb_from_pdbtbx(&pdb, &local_input)?;
                    legacy_flag = "--pdb";
                }
            }

            let log_path = work.path.join("legacy.log");
            let log = std::fs::File::create(&log_path)?;
            let log_err = log.try_clone()?;

            let mut cmd = std::process::Command::new(&legacy_bin);
            cmd.arg(legacy_flag);
            cmd.arg(&input_arg);
            cmd.current_dir(&work.path);
            cmd.env("RNAVIEW", &rnaview_root);
            cmd.stdout(Stdio::from(log));
            cmd.stderr(Stdio::from(log_err));

            let status = cmd.status()?;
            if !status.success() {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("legacy rnaview failed (status={status}); see {log_path:?}"),
                )));
            }

            let produced_out = work.path.join(format!("{input_arg}.out"));
            if !produced_out.exists() {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("legacy rnaview produced no .out; see {log_path:?}"),
                )));
            }

            if let Some(dst) = emit_out {
                std::fs::copy(&produced_out, dst)?;
            }

            let core: Core = extract_core_from_out_path_with_out_index(&produced_out)?;
            let pairs = PairsJson {
                schema_version: 1,
                source: Some(Source {
                    path: input.to_string_lossy().to_string(),
                    format: fmt.clone(),
                    id_scheme: if fmt == "cif" { Some("auth".to_string()) } else { None },
                    model: source_model,
                }),
                options: Some(serde_json::json!({"engine":"rust","oracle":"legacy","mmcif_parser":mmcif_parser})),
                core,
            };

            let json_text = serde_json::to_string(&serde_json::to_value(&pairs)?)? + "\n";
            if let Some(out_path) = output {
                std::fs::write(out_path, json_text)?;
            } else {
                print!("{json_text}");
            }
            Ok(())
        }
        "write-out" => {
            if args.is_empty() {
                usage();
            }
            let input = PathBuf::from(args.remove(0));
            let mut output: Option<PathBuf> = None;
            while !args.is_empty() {
                let flag = args.remove(0);
                if flag == "-o" || flag == "--output" {
                    if args.is_empty() {
                        usage();
                    }
                    output = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                usage();
            }

            let pairs: PairsJson = serde_json::from_str(&std::fs::read_to_string(&input)?)?;
            let out_text = write_out_core(&pairs.core);
            if let Some(out_path) = output {
                std::fs::write(out_path, out_text)?;
            } else {
                print!("{out_text}");
            }
            Ok(())
        }
        "render-wrl" => {
            if args.is_empty() {
                usage();
            }
            let input = PathBuf::from(args.remove(0));
            let mut output: Option<PathBuf> = None;
            let mut source_override: Option<PathBuf> = None;

            while !args.is_empty() {
                let flag = args.remove(0);
                if flag == "-o" || flag == "--output" {
                    if args.is_empty() {
                        usage();
                    }
                    output = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                if flag == "--source" {
                    if args.is_empty() {
                        usage();
                    }
                    source_override = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                usage();
            }

            let pairs: PairsJson = serde_json::from_str(&std::fs::read_to_string(&input)?)?;
            let source_path = if let Some(p) = source_override {
                p
            } else if let Some(src) = pairs.source.as_ref() {
                PathBuf::from(&src.path)
            } else {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "pairs.json is missing .source.path; pass --source",
                )));
            };

            let rnaview_root: PathBuf = match std::env::var_os("RNAVIEW") {
                Some(p) => PathBuf::from(p),
                None => PathBuf::from("."),
            };
            let basepars_vrml = rnaview_root.join("BASEPARS").join("vrml_image.par");

            let policies = SemanticsConfig::defaults(Semantics::LegacyV1).policies.structure;
            let wrl = render_vrml_from_pairs_json(&pairs, &source_path, &policies, &basepars_vrml)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

            if let Some(out_path) = output {
                std::fs::write(out_path, wrl)?;
            } else {
                print!("{wrl}");
            }
            Ok(())
        }
        "render-2d" => {
            if args.is_empty() {
                usage();
            }
            let input = PathBuf::from(args.remove(0));
            let mut out_xml: Option<PathBuf> = None;
            let mut out_ps: Option<PathBuf> = None;
            let mut source_override: Option<PathBuf> = None;

            while !args.is_empty() {
                let flag = args.remove(0);
                if flag == "--source" {
                    if args.is_empty() {
                        usage();
                    }
                    source_override = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                if flag == "--out-xml" {
                    if args.is_empty() {
                        usage();
                    }
                    out_xml = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                if flag == "--out-ps" {
                    if args.is_empty() {
                        usage();
                    }
                    out_ps = Some(PathBuf::from(args.remove(0)));
                    continue;
                }
                usage();
            }

            if out_xml.is_none() && out_ps.is_none() {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "render-2d requires --out-xml and/or --out-ps",
                )));
            }

            let pairs: PairsJson = serde_json::from_str(&std::fs::read_to_string(&input)?)?;
            let source_path = if let Some(p) = source_override {
                p
            } else if let Some(src) = pairs.source.as_ref() {
                PathBuf::from(&src.path)
            } else {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "pairs.json is missing .source.path; pass --source",
                )));
            };

            let rnaview_root: PathBuf = match std::env::var_os("RNAVIEW") {
                Some(p) => PathBuf::from(p),
                None => PathBuf::from("."),
            };
            let basepars_ps_image = rnaview_root.join("BASEPARS").join("ps_image.par");

            let structure_policies = SemanticsConfig::defaults(Semantics::LegacyV1).policies.structure;
            let out = rnaview_hotcore::render_2d_from_pairs_json(
                &pairs,
                &source_path,
                &structure_policies,
                &basepars_ps_image,
            )
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

            if let Some(p) = out_xml {
                std::fs::write(p, out.xml)?;
            }
            if let Some(p) = out_ps {
                std::fs::write(p, out.ps)?;
            }
            Ok(())
        }
        _ => usage(),
    }
}
