set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export RNAVIEW="$ROOT_DIR"
export PATH="$ROOT_DIR/bin:$PATH"

if [[ ! -x "$ROOT_DIR/bin/rnaview" ]]; then
  bash "$ROOT_DIR/tools/build_legacy_rnaview.sh"
fi

OUT_DIR="$(mktemp -d)"
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
  --regress \
  --regress-mode out \
  --keep-going; then
  rm -rf "$OUT_DIR"
  exit 0
fi

echo "FAILED: outputs kept at $OUT_DIR" >&2
exit 1
