#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
SYSROOT="${HOME}/.cache/rnaview-toolchain/root"

if [[ ! -x "${SYSROOT}/usr/bin/cc" ]]; then
  bash "${REPO_ROOT}/tools/build_legacy_rnaview.sh"
fi

export PATH="${SYSROOT}/usr/bin:${SYSROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${SYSROOT}/usr/lib/x86_64-linux-gnu:${SYSROOT}/lib/x86_64-linux-gnu:${SYSROOT}/usr/lib:${SYSROOT}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

echo "== build rust staticlib (legacy-ffi) ==" >&2
bash "${REPO_ROOT}/tools/cargo_sysroot.sh" build --manifest-path "${REPO_ROOT}/rust/Cargo.toml" --lib --features legacy-ffi

RUST_LIB="${REPO_ROOT}/rust/target/debug/librnaview_hotcore.a"
if [[ ! -f "${RUST_LIB}" ]]; then
  echo "missing rust staticlib: ${RUST_LIB}" >&2
  exit 2
fi

OBJ_DIR="${REPO_ROOT}/obj/rustcore"
BIN_OUT="${REPO_ROOT}/bin/rnaview_rustcore"
mkdir -p "${OBJ_DIR}" "${REPO_ROOT}/bin"

CC="${SYSROOT}/usr/bin/gcc-14"
CFLAGS=(
  "--sysroot=${SYSROOT}"
  "-isystem" "${SYSROOT}/usr/include"
  "-isystem" "${SYSROOT}/usr/include/x86_64-linux-gnu"
  "-I${REPO_ROOT}/include"
  "-g"
  "-Wall"
  "-DRNAVIEW_RUST_CHECK_PAIRS"
  "-DRNAVIEW_RUST_HBOND_PAIR"
  "-DRNAVIEW_RUST_LW_PAIR_TYPE"
)

sources=(
  "rnaview.c"
  "fpair.c"
  "fpair_sub.c"
  "pair_type.c"
  "nrutil.c"
  "ps-xy.c"
  "ps-xy-sub.c"
  "vrml.c"
  "rnaxml-new.c"
  "analyze.c"
  "pattern.c"
  "xml2ps.c"
  "multiple.c"
  "statistics.c"
)

objs=()
for src in "${sources[@]}"; do
  obj="${OBJ_DIR}/${src%.c}.o"
  objs+=("${obj}")
  "${CC}" "${CFLAGS[@]}" -c "${REPO_ROOT}/src/${src}" -o "${obj}"
done

echo "== link ${BIN_OUT} ==" >&2
"${CC}" \
  "--sysroot=${SYSROOT}" \
  -g -Wall \
  -o "${BIN_OUT}" \
  "${objs[@]}" \
  "${RUST_LIB}" \
  -lm -ldl -lpthread

echo "built: ${BIN_OUT}" >&2
