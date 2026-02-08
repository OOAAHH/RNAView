set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export RNAVIEW="$ROOT_DIR"

# Gate B (No-C) acceptance: do not build/link any legacy C binaries.
# This runs the Rust toolchain + Python orchestration only, and performs
# byte-exact .out regression against the frozen golden set.

if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
  bash "$ROOT_DIR/tools/cargo_sysroot.sh" build --manifest-path "$ROOT_DIR/rust/Cargo.toml"
  bash "$ROOT_DIR/tools/cargo_sysroot.sh" test --manifest-path "$ROOT_DIR/rust/Cargo.toml"
  python3 -m unittest discover -s "$ROOT_DIR/tools" -p "test_*.py"
fi

echo "== rust engine (No-C: structure -> Rust compute -> .out(full)) ==" >&2
OUT_DIR="${OUT_DIR:-}"
_CLEAN_OUT_DIR_ON_SUCCESS=0
if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="$(mktemp -d)"
  _CLEAN_OUT_DIR_ON_SUCCESS=1
else
  mkdir -p "$OUT_DIR"
fi
echo "Gate B out dir: $OUT_DIR" >&2
if python3 "$ROOT_DIR/tools/rnaview_batch.py" run \
  test/pdb/pdb1nvy/pdb1nvy.pdb \
  test/pdb/test1/test1.pdb \
  test/pdb/tr0001/tr0001.pdb \
  test/pdb/url064/url064.pdb \
  test/pdb/urx053/urx053.pdb \
  test/mmcif/insertion_code/1EFW/1EFW.cif \
  test/mmcif/insertion_code/1VVJ/1VVJ.cif \
  test/mmcif/insertion_code/4ARC/4ARC.cif \
  test/mmcif/nmr_structure/8if5/8if5.cif \
  test/mmcif/other/6pom/6pom.cif \
  test/mmcif/x-ray/3P4J/assembly-1/3p4j-assembly1.cif \
  test/mmcif/x-ray/434D/assembly-1/434d-assembly1.cif \
  test/mmcif/x-ray/434D/assembly-2/434d-assembly2.cif \
  test/mmcif/x-ray/4NMG/assembly-1/4nmg-assembly1.cif \
  --out-dir "$OUT_DIR" \
  --engine rust \
  --rust-oracle compute \
  --regress \
  --regress-mode out \
  --keep-going; then
  if [[ "${_CLEAN_OUT_DIR_ON_SUCCESS}" -eq 1 ]]; then
    rm -rf "$OUT_DIR"
  fi
  exit 0
fi

echo "FAILED: outputs kept at $OUT_DIR" >&2
exit 1
