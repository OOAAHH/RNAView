#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
SYSROOT="${HOME}/.cache/rnaview-toolchain/root"

if [[ ! -x "${SYSROOT}/usr/bin/cc" ]]; then
  bash "${REPO_ROOT}/tools/build_legacy_rnaview.sh"
fi

export PATH="${SYSROOT}/usr/bin:${SYSROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${SYSROOT}/usr/lib/x86_64-linux-gnu:${SYSROOT}/lib/x86_64-linux-gnu:${SYSROOT}/usr/lib:${SYSROOT}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER="${SYSROOT}/usr/bin/gcc-14"
export RUSTFLAGS="-C link-arg=--sysroot=${SYSROOT}${RUSTFLAGS:+ ${RUSTFLAGS}}"

exec cargo "$@"

