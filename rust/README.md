# Rust hot core (scaffold)

This crate is the starting point for the Rust “hot core engine”.

For Phase 0 / Phase 1 it focuses on **I/O contracts**:

- Emit `pairs.json` (schema v1) from a legacy RNAVIEW `.out` file (core sections only).
- Render a minimal RNAVIEW `.out` containing only the core sections from `pairs.json` (`pairs.json -> .out(core)` as a pure function).

## Phase 2 "No-C" target

The Phase 2 end-goal is **byte-exact legacy `.out` compatibility without any C/legacy runtime dependency**
(pure Rust core + Python orchestration).

That means:

- `from-structure` must not shell out to the legacy `bin/rnaview` (legacy is allowed only as a test oracle).
- Writing `FILEOUT.out` must not rely on the legacy C writer (a No-C `.out(full)` writer is required).

## Build

```bash
cargo build --manifest-path rust/Cargo.toml
```

If your environment is missing a system linker/dev libs (common in minimal containers), use the vendored sysroot toolchain:

```bash
bash tools/cargo_sysroot.sh build --manifest-path rust/Cargo.toml
bash tools/cargo_sysroot.sh test --manifest-path rust/Cargo.toml
```

## Usage

### `PDB/mmCIF -> pairs.json` (Phase 2 bootstrap)

By default this shells out to the legacy `bin/rnaview` as an oracle, then parses the generated `.out` into `pairs.json`.

```bash
RNAVIEW="$(pwd)" bash tools/cargo_sysroot.sh run --manifest-path rust/Cargo.toml -- \
  from-structure test/pdb/tr0001/tr0001.pdb \
  -o /tmp/tr0001.pairs.json \
  --emit-out /tmp/tr0001.out
```

There is also a No-C dev mode which reads a sidecar oracle at `<input>.out` (e.g. `test/pdb/tr0001/tr0001.pdb.out`):

```bash
bash tools/cargo_sysroot.sh run --manifest-path rust/Cargo.toml -- \
  from-structure test/pdb/tr0001/tr0001.pdb \
  --oracle out \
  -o /tmp/tr0001.pairs.json \
  --emit-out /tmp/tr0001.out
```

### Semantics & policies (Phase 3+)

For `--oracle compute` you can choose an explicit semantics preset, and override individual policies.
These are written to `pairs.json.options` for traceability.

```bash
bash tools/cargo_sysroot.sh run --manifest-path rust/Cargo.toml -- \
  from-structure test/mmcif/nmr_structure/8if5/8if5.cif \
  --oracle compute \
  --semantics science-v1 \
  --hydrogen-policy discard-all \
  -o /tmp/8if5.pairs.json \
  --emit-out /tmp/8if5.out
```

### `.out -> pairs.json`

```bash
cargo run --manifest-path rust/Cargo.toml -- from-out test/pdb/tr0001/tr0001.pdb.out -o /tmp/tr0001.pairs.json
```

### `pairs.json -> .out(core)`

```bash
cargo run --manifest-path rust/Cargo.toml -- write-out /tmp/tr0001.pairs.json -o /tmp/tr0001.core.out
```

### Render (Phase 4)

Render commands take `pairs.json` as the contract and may re-read the source structure file (`--source`) for coordinates/labels.
For determinism and compatibility, pass the same `--semantics` (and policy overrides, if any) that were used to generate `pairs.json`.

```bash
bash tools/cargo_sysroot.sh run --manifest-path rust/Cargo.toml -- \
  render-2d /tmp/tr0001.pairs.json \
  --source test/pdb/tr0001/tr0001.pdb \
  --semantics legacy-v1 \
  --out-xml /tmp/tr0001.xml --out-ps /tmp/tr0001.ps

bash tools/cargo_sysroot.sh run --manifest-path rust/Cargo.toml -- \
  render-wrl /tmp/tr0001.pairs.json \
  --source test/pdb/tr0001/tr0001.pdb \
  --semantics legacy-v1 \
  -o /tmp/tr0001.wrl
```

### Reuse the existing core regression

```bash
python3 tools/rnaview_out_core.py compare test/pdb/tr0001/tr0001.pdb.out /tmp/tr0001.core.out
```

Or validate the writer against the frozen golden set:

```bash
python3 tools/rnaview_pairs_json.py validate-golden
```

## Phase 2 regressions

- Gate A (legacy + rustcore + rust; byte-exact `.out` regression): `bash test_phase2.sh`
- Gate B (No-C; Rust+Python only): `bash test_phase2_noc.sh`

## Phase 3/4 regressions

- Gate C (science diff + allowlist): `bash test_phase3_gate_c.sh`
- science-v1 frozen core regression: `bash test_phase3_science.sh`
- Gate D (render canonical regression): `bash test_phase4_gate_d.sh`
- Gate NA (DNA/RNA hybrid, science-v1 self-oracle): `bash test_phase4_gate_na.sh`
