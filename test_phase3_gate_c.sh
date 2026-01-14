#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

OUT_DIR="${OUT_DIR:-out_phase3_gate_c}"
WORKERS="${WORKERS:-2}"

rm -rf "$OUT_DIR"
python3 tools/rnaview_gate_c.py run test/pdb test/mmcif --out-dir "$OUT_DIR" --workers "$WORKERS"
echo "Gate C report: $OUT_DIR/report.md"
echo "Gate C summary: $OUT_DIR/summary.json"
