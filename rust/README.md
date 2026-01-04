# Rust hot core (scaffold)

This crate is the starting point for the Rust “hot core engine”.

For Phase 0 / Phase 1 it focuses on **I/O contracts**:

- Emit `pairs.json` (schema v1) from a legacy RNAVIEW `.out` file (core sections only).
- Render a minimal RNAVIEW `.out` containing only the core sections from `pairs.json` (`pairs.json -> .out(core)` as a pure function).

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

### `PDB/mmCIF -> pairs.json` (Phase2 bootstrap)

This currently shells out to the legacy `bin/rnaview` as an oracle, then parses the generated `.out` into `pairs.json`.

```bash
RNAVIEW="$(pwd)" bash tools/cargo_sysroot.sh run --manifest-path rust/Cargo.toml -- \
  from-structure test/pdb/tr0001/tr0001.pdb -o /tmp/tr0001.pairs.json
```

### `.out -> pairs.json`

```bash
cargo run --manifest-path rust/Cargo.toml -- from-out test/pdb/tr0001/tr0001.pdb.out -o /tmp/tr0001.pairs.json
```

### `pairs.json -> .out(core)`

```bash
cargo run --manifest-path rust/Cargo.toml -- write-out /tmp/tr0001.pairs.json -o /tmp/tr0001.core.out
```

### Reuse the existing core regression

```bash
python3 tools/rnaview_out_core.py compare test/pdb/tr0001/tr0001.pdb.out /tmp/tr0001.core.out
```

Or validate the writer against the frozen golden set:

```bash
python3 tools/rnaview_pairs_json.py validate-golden
```
