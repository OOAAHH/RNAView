#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export RNAVIEW="$ROOT_DIR"

# Phase 3 (science-v1) acceptance:
# - No-C: pure Rust core + Python orchestration
# - core regression against frozen science-v1 golden set

bash "$ROOT_DIR/tools/cargo_sysroot.sh" build --manifest-path "$ROOT_DIR/rust/Cargo.toml"
bash "$ROOT_DIR/tools/cargo_sysroot.sh" test --manifest-path "$ROOT_DIR/rust/Cargo.toml"
python3 -m unittest discover -s "$ROOT_DIR/tools" -p "test_*.py"

echo "== science-v1 engine (No-C: structure -> Rust compute -> core regression) ==" >&2
OUT_DIR="$(mktemp -d)"
if python3 "$ROOT_DIR/tools/rnaview_batch.py" run \
  test/pdb \
  test/mmcif \
  --out-dir "$OUT_DIR" \
  --engine rust \
  --rust-oracle compute \
  --semantics science-v1 \
  --regress \
  --manifest "$ROOT_DIR/test/golden_science_core/manifest.json" \
  --regress-mode core \
  --keep-going; then
  rm -rf "$OUT_DIR"
  exit 0
fi

echo "FAILED: outputs kept at $OUT_DIR" >&2
exit 1

