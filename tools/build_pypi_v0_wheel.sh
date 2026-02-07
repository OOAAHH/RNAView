#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if ! command -v cargo >/dev/null 2>&1; then
  echo "missing cargo; cannot build rnaview-hotcore" >&2
  exit 2
fi

if command -v cc >/dev/null 2>&1; then
  cargo build --release --manifest-path rust/Cargo.toml
else
  bash tools/cargo_sysroot.sh build --release --manifest-path rust/Cargo.toml
fi

HOTCORE_BIN="rust/target/release/rnaview-hotcore"
if [[ ! -x "$HOTCORE_BIN" ]]; then
  echo "missing built binary: $HOTCORE_BIN" >&2
  exit 2
fi

install -d rnaview/_bin dist
install -m 0755 "$HOTCORE_BIN" rnaview/_bin/rnaview-hotcore

python3 -m pip wheel . -w dist --no-deps
echo "Wheels written to: dist/"

