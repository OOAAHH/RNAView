#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export RNAVIEW="$ROOT_DIR"

usage() {
  cat >&2 <<'EOF'
Usage:
  bash test_phase3_wrapup.sh [options] [extra_inputs...]

Default (no args):
  - Runs Phase 3 closeout on the repository targets:
    - Gate B (legacy-v1, No-C, .out byte-exact)
    - Gate C (science diff + allowlist)
    - science-v1 regression (frozen core goldens)

Extra inputs (files/dirs/globs):
  - Additionally runs legacy-v1 + science-v1 on the given inputs and produces:
    - per-semantics outputs (pairs.json + .out)
    - a Gate C diff report (legacy-v1 vs science-v1)

Options:
  --only-extra         Skip repository closeout gates; run only extra inputs
  --strict-extra       Fail if extra-input Gate C reports unapproved diffs
  --bench              Run legacy vs rustcore benchmark (non-gating)
  --workers N          Parallel workers for Gate C (default: env WORKERS or 2)
  --out-extra-dir DIR  Where to write extra-input outputs (default: mktemp under repo)
  -h, --help           Show this help

Examples:
  bash test_phase3_wrapup.sh
  bash test_phase3_wrapup.sh /path/to/new.pdb
  bash test_phase3_wrapup.sh test/mmcif/x-ray/434D/assembly-1/434d-assembly1.cif
  bash test_phase3_wrapup.sh --only-extra --out-extra-dir out/mycase ./new.cif
EOF
}

ONLY_EXTRA=0
STRICT_EXTRA=0
RUN_BENCH=0
WORKERS="${WORKERS:-2}"
OUT_EXTRA_DIR="${OUT_EXTRA_DIR:-}"

extra_inputs=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only-extra)
      ONLY_EXTRA=1
      shift
      ;;
    --strict-extra)
      STRICT_EXTRA=1
      shift
      ;;
    --bench)
      RUN_BENCH=1
      shift
      ;;
    --workers)
      WORKERS="${2:-}"
      if [[ -z "$WORKERS" ]]; then
        echo "missing value for --workers" >&2
        exit 2
      fi
      shift 2
      ;;
    --out-extra-dir)
      OUT_EXTRA_DIR="${2:-}"
      if [[ -z "$OUT_EXTRA_DIR" ]]; then
        echo "missing value for --out-extra-dir" >&2
        exit 2
      fi
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      extra_inputs+=("$@")
      break
      ;;
    *)
      extra_inputs+=("$1")
      shift
      ;;
  esac
done

if [[ "${ONLY_EXTRA}" -eq 0 ]]; then
  echo "== phase3: build + unit tests ==" >&2
  bash "$ROOT_DIR/tools/cargo_sysroot.sh" build --manifest-path "$ROOT_DIR/rust/Cargo.toml"
  bash "$ROOT_DIR/tools/cargo_sysroot.sh" test --manifest-path "$ROOT_DIR/rust/Cargo.toml"
  python3 -m unittest discover -s "$ROOT_DIR/tools" -p "test_*.py"

  echo "== phase3: Gate B (legacy-v1, No-C, .out byte-exact) ==" >&2
  SKIP_BUILD=1 bash "$ROOT_DIR/test_phase2_noc.sh"

  echo "== phase3: Gate C (science diff + allowlist) ==" >&2
  WORKERS="$WORKERS" bash "$ROOT_DIR/test_phase3_gate_c.sh"

  echo "== phase3: science-v1 regression (frozen core goldens) ==" >&2
  SKIP_BUILD=1 bash "$ROOT_DIR/test_phase3_science.sh"

  if [[ "${RUN_BENCH}" -eq 1 ]]; then
    echo "== phase3: bench (non-gating) ==" >&2
    bash "$ROOT_DIR/tools/build_legacy_rnaview.sh"
    bash "$ROOT_DIR/tools/build_rnaview_rustcore_release.sh"
    python3 "$ROOT_DIR/tools/rnaview_bench.py" compare \
      --suite phase2 \
      --runs 3 \
      --legacy-bin "$ROOT_DIR/bin/rnaview" \
      --rustcore-bin "$ROOT_DIR/bin/rnaview_rustcore_release" \
      -o "$ROOT_DIR/out/bench_phase3.json" || true
    echo "bench json: $ROOT_DIR/out/bench_phase3.json" >&2
  fi
fi

if [[ "${#extra_inputs[@]}" -gt 0 ]]; then
  if [[ -z "${OUT_EXTRA_DIR}" ]]; then
    OUT_EXTRA_DIR="$(mktemp -d "$ROOT_DIR/out_phase3_extra.XXXXXX")"
  else
    mkdir -p "$OUT_EXTRA_DIR"
    OUT_EXTRA_DIR="$(cd "$OUT_EXTRA_DIR" && pwd)"
  fi

  echo "== phase3: extra inputs ==" >&2
  echo "extra out dir: $OUT_EXTRA_DIR" >&2

  echo "== extra: legacy-v1 (No-C compute + .out(full)) ==" >&2
  python3 "$ROOT_DIR/tools/rnaview_batch.py" run \
    "${extra_inputs[@]}" \
    --out-dir "$OUT_EXTRA_DIR/legacy-v1" \
    --engine rust \
    --rust-oracle compute \
    --semantics legacy-v1 \
    --keep-going

  echo "== extra: science-v1 (No-C compute + .out(full)) ==" >&2
  python3 "$ROOT_DIR/tools/rnaview_batch.py" run \
    "${extra_inputs[@]}" \
    --out-dir "$OUT_EXTRA_DIR/science-v1" \
    --engine rust \
    --rust-oracle compute \
    --semantics science-v1 \
    --keep-going

  echo "== extra: Gate C diff (legacy-v1 vs science-v1) ==" >&2
  set +e
  python3 "$ROOT_DIR/tools/rnaview_gate_c.py" run \
    "${extra_inputs[@]}" \
    --out-dir "$OUT_EXTRA_DIR/gate_c" \
    --workers "$WORKERS"
  gate_c_status=$?
  set -e

  echo "extra Gate C report: $OUT_EXTRA_DIR/gate_c/report.md" >&2
  echo "extra Gate C summary: $OUT_EXTRA_DIR/gate_c/summary.json" >&2

  if [[ "$gate_c_status" -ne 0 ]]; then
    echo "WARNING: extra Gate C had unapproved diffs or failures; inspect the report above." >&2
    if [[ "${STRICT_EXTRA}" -eq 1 ]]; then
      exit "$gate_c_status"
    fi
  fi
fi

echo "OK" >&2
