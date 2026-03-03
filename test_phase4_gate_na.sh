#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

export RNAVIEW="$ROOT"

OUT_DIR="${OUT_DIR:-out_phase4_gate_na}"
SEMANTICS="${SEMANTICS:-science-v1}"
GOLDEN_DIR="${GOLDEN_DIR:-test/golden_na}"

# Ensure hotcore exists (Gate NA is self-oracle; no legacy/C required).
if [[ ! -x "$ROOT/rust/target/release/rnaview-hotcore" && ! -x "$ROOT/rust/target/debug/rnaview-hotcore" ]]; then
  bash "$ROOT/tools/cargo_sysroot.sh" build --manifest-path "$ROOT/rust/Cargo.toml"
fi

python3 tools/rnaview_gate_na.py compare \
  --out-dir "$OUT_DIR" \
  --golden-dir "$GOLDEN_DIR" \
  --semantics "$SEMANTICS" \
  test/hybrid

echo "Gate NA report: $OUT_DIR/report.md"
echo "Gate NA summary: $OUT_DIR/summary.json"

