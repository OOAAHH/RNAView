use rnaview_hotcore::{extract_core_from_out_path, write_out_core, Core, PairsJson, Source};
use std::path::{Path, PathBuf};
use std::process::Stdio;

fn usage() -> ! {
    eprintln!(
        "Usage:\n  rnaview-hotcore from-out <file.out> [-o pairs.json]\n  rnaview-hotcore from-structure <file.pdb|file.cif> [--format pdb|cif] [-o pairs.json]\n  rnaview-hotcore write-out <pairs.json> [-o candidate.out]"
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

            let core: Core = extract_core_from_out_path(&input)?;
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
            let mut emit_legacy_out: Option<PathBuf> = None;

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
                if flag == "--emit-legacy-out" {
                    if args.is_empty() {
                        usage();
                    }
                    emit_legacy_out = Some(PathBuf::from(args.remove(0)));
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
            let input_name = input
                .file_name()
                .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidInput, "input has no file name"))?;
            let local_input = work.path.join(input_name);
            std::fs::copy(&input, &local_input)?;

            let log_path = work.path.join("legacy.log");
            let log = std::fs::File::create(&log_path)?;
            let log_err = log.try_clone()?;

            let mut cmd = std::process::Command::new(&legacy_bin);
            if fmt == "pdb" {
                cmd.arg("--pdb");
            } else if fmt == "cif" {
                cmd.arg("--cif");
            }
            cmd.arg(input_name);
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

            let produced_out = work.path.join(format!("{}.out", input_name.to_string_lossy()));
            if !produced_out.exists() {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("legacy rnaview produced no .out; see {log_path:?}"),
                )));
            }

            if let Some(dst) = emit_legacy_out {
                std::fs::copy(&produced_out, dst)?;
            }

            let core: Core = extract_core_from_out_path(&produced_out)?;
            let pairs = PairsJson {
                schema_version: 1,
                source: Some(Source {
                    path: input.to_string_lossy().to_string(),
                    format: fmt.clone(),
                    id_scheme: None,
                    model: None,
                }),
                options: Some(serde_json::json!({"engine":"rust","oracle":"legacy"})),
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
        _ => usage(),
    }
}
