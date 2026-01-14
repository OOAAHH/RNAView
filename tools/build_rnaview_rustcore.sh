#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
SYSROOT="${HOME}/.cache/rnaview-toolchain/root"

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

CC=""
LD_FLAGS=()
SYS_CC=""
if [[ -x "${SYSROOT}/usr/bin/gcc-14" ]]; then
  SYS_CC="${SYSROOT}/usr/bin/gcc-14"
elif [[ -x "${SYSROOT}/usr/bin/cc" ]]; then
  SYS_CC="${SYSROOT}/usr/bin/cc"
fi

if [[ -n "${SYS_CC}" ]]; then
  CC="${SYS_CC}"
  export PATH="${SYSROOT}/usr/bin:${SYSROOT}/bin:${PATH}"
  export LD_LIBRARY_PATH="${SYSROOT}/usr/lib/x86_64-linux-gnu:${SYSROOT}/lib/x86_64-linux-gnu:${SYSROOT}/usr/lib:${SYSROOT}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
  LD_FLAGS=("--sysroot=${SYSROOT}")
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
    "-DRNAVIEW_RUST_CANDIDATE_PAIRS"
  )
else
  if command -v cc >/dev/null 2>&1; then
    CC="$(command -v cc)"
  elif command -v gcc >/dev/null 2>&1; then
    CC="$(command -v gcc)"
  elif command -v clang >/dev/null 2>&1; then
    CC="$(command -v clang)"
  else
    echo "missing C compiler (need cc/gcc/clang or sysroot toolchain)" >&2
    exit 2
  fi
  CFLAGS=(
    "-I${REPO_ROOT}/include"
    "-g"
    "-Wall"
    "-DRNAVIEW_RUST_CHECK_PAIRS"
    "-DRNAVIEW_RUST_HBOND_PAIR"
    "-DRNAVIEW_RUST_LW_PAIR_TYPE"
    "-DRNAVIEW_RUST_CANDIDATE_PAIRS"
  )
fi

sources=(
  "rnaview.c"
  "fpair.c"
  "fpair_sub.c"
  "pair_type.c"
  "rnaview_profile.c"
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
  "${LD_FLAGS[@]}" \
  -g -Wall \
  -o "${BIN_OUT}" \
  "${objs[@]}" \
  "${RUST_LIB}" \
  -lm -ldl -lpthread

echo "built: ${BIN_OUT}" >&2
