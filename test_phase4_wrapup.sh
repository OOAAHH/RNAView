#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export RNAVIEW="$ROOT_DIR"

echo "== phase4: build + unit tests ==" >&2
bash "$ROOT_DIR/tools/cargo_sysroot.sh" test --manifest-path "$ROOT_DIR/rust/Cargo.toml"
python3 -m unittest discover -s "$ROOT_DIR/tools" -p "test_*.py"

echo "== phase4: Gate D (render canonical regression) ==" >&2
bash "$ROOT_DIR/test_phase4_gate_d.sh"

echo "OK" >&2

