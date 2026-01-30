#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

export RNAVIEW="$ROOT"

OUT_DIR="${OUT_DIR:-out_phase4_gate_d}"
SEMANTICS="${SEMANTICS:-legacy-v1}"
RENDERER_VERSION="${RENDERER_VERSION:-dev}"
EMIT_GOLDEN="${EMIT_GOLDEN:-none}"
GATE_D_FORMATS="${GATE_D_FORMATS:-}"
CANDIDATE_ENGINE="${CANDIDATE_ENGINE:-command}"
CANDIDATE_CMD="${CANDIDATE_CMD:-}"
CANDIDATE_BACKEND="${CANDIDATE_BACKEND:-}"

if [[ "$CANDIDATE_ENGINE" == "legacy" ]]; then
  if [[ ! -x "$ROOT/bin/rnaview" ]]; then
    bash "$ROOT/tools/build_legacy_rnaview.sh"
  fi
else
  # candidate_engine=command
  if [[ -z "$CANDIDATE_CMD" && -n "$CANDIDATE_BACKEND" ]]; then
    CANDIDATE_CMD="python3 tools/rnaview_render.py render --backend ${CANDIDATE_BACKEND}"
  fi
  if [[ -z "$CANDIDATE_BACKEND" && -n "$CANDIDATE_CMD" ]]; then
    if [[ "$CANDIDATE_CMD" == *"--backend rustcore-release"* ]]; then
      CANDIDATE_BACKEND="rustcore-release"
    elif [[ "$CANDIDATE_CMD" == *"--backend rustcore"* ]]; then
      CANDIDATE_BACKEND="rustcore"
    elif [[ "$CANDIDATE_CMD" == *"--backend pairs-json"* ]]; then
      CANDIDATE_BACKEND="pairs-json"
    elif [[ "$CANDIDATE_CMD" == *"--backend pairs-out-noc3d"* ]]; then
      CANDIDATE_BACKEND="pairs-out-noc3d"
    elif [[ "$CANDIDATE_CMD" == *"--backend pairs-out"* ]]; then
      CANDIDATE_BACKEND="pairs-out"
    elif [[ "$CANDIDATE_CMD" == *"--backend legacy"* ]]; then
      CANDIDATE_BACKEND="legacy"
    fi
  fi
  case "${CANDIDATE_BACKEND:-legacy}" in
    legacy)
      if [[ -z "$CANDIDATE_CMD" || "$CANDIDATE_CMD" == *"tools/rnaview_render.py render"* ]]; then
        if [[ ! -x "$ROOT/bin/rnaview" ]]; then
          bash "$ROOT/tools/build_legacy_rnaview.sh"
        fi
      fi
      ;;
    rustcore)
      if [[ ! -x "$ROOT/bin/rnaview_rustcore" ]]; then
        bash "$ROOT/tools/build_rnaview_rustcore.sh"
      fi
      ;;
    rustcore-release)
      if [[ ! -x "$ROOT/bin/rnaview_rustcore_release" ]]; then
        bash "$ROOT/tools/build_rnaview_rustcore_release.sh"
      fi
      ;;
    pairs-out)
      if [[ ! -x "$ROOT/bin/rnaview_rustcore_release" ]]; then
        bash "$ROOT/tools/build_rnaview_rustcore_release.sh"
      fi
      ;;
    pairs-json)
      if [[ ! -x "$ROOT/bin/rnaview_rustcore_release" ]]; then
        bash "$ROOT/tools/build_rnaview_rustcore_release.sh"
      fi
      ;;
    pairs-out-noc3d)
      if [[ ! -x "$ROOT/bin/rnaview_rustcore_release" ]]; then
        bash "$ROOT/tools/build_rnaview_rustcore_release.sh"
      fi
      ;;
    *)
      echo "unknown CANDIDATE_BACKEND: ${CANDIDATE_BACKEND}" >&2
      exit 2
      ;;
  esac
fi

cmd=(python3 tools/rnaview_gate_d.py compare \
  --out-dir "$OUT_DIR" \
  --golden-dir test/golden_render \
  --allowlist test/gate_d_allowlist.yaml \
  --candidate-engine "$CANDIDATE_ENGINE" \
  --semantics "$SEMANTICS" \
  --renderer-version "$RENDERER_VERSION" \
  --emit-golden "$EMIT_GOLDEN")

if [[ -n "$GATE_D_FORMATS" ]]; then
  cmd+=(--formats "$GATE_D_FORMATS")
fi

if [[ -n "$CANDIDATE_CMD" ]]; then
  read -r -a CANDIDATE_CMD_ARR <<<"$CANDIDATE_CMD"
  cmd+=(--candidate-cmd "${CANDIDATE_CMD_ARR[@]}")
fi

"${cmd[@]}"

echo "Gate D report: $OUT_DIR/report.md"
echo "Gate D summary: $OUT_DIR/summary.json"
