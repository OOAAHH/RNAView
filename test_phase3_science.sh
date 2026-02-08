#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export RNAVIEW="$ROOT_DIR"

# Phase 3 (science-v1) acceptance:
# - No-C: pure Rust core + Python orchestration
# - core regression against frozen science-v1 golden set

if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
  bash "$ROOT_DIR/tools/cargo_sysroot.sh" build --manifest-path "$ROOT_DIR/rust/Cargo.toml"
  bash "$ROOT_DIR/tools/cargo_sysroot.sh" test --manifest-path "$ROOT_DIR/rust/Cargo.toml"
  python3 -m unittest discover -s "$ROOT_DIR/tools" -p "test_*.py"
fi

echo "== science-v1 engine (No-C: structure -> Rust compute -> core regression) ==" >&2
OUT_DIR="${OUT_DIR:-}"
_CLEAN_OUT_DIR_ON_SUCCESS=0
if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="$(mktemp -d)"
  _CLEAN_OUT_DIR_ON_SUCCESS=1
else
  mkdir -p "$OUT_DIR"
fi
echo "science out dir: $OUT_DIR" >&2
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
  if [[ "${_CLEAN_OUT_DIR_ON_SUCCESS}" -eq 1 ]]; then
    rm -rf "$OUT_DIR"
  fi
  exit 0
fi

echo "FAILED: outputs kept at $OUT_DIR" >&2
exit 1
