#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
SYSROOT="${HOME}/.cache/rnaview-toolchain/root"

SYS_LINKER="${SYSROOT}/usr/bin/gcc-14"
if [[ ! -x "${SYS_LINKER}" ]]; then
  SYS_LINKER="${SYSROOT}/usr/bin/cc"
fi

# Prefer the local sysroot toolchain when present (devcontainer/No-cc env).
# Fall back to the system toolchain when available (CI/dev machines).
if [[ ! -x "${SYS_LINKER}" ]]; then
  if command -v cc >/dev/null 2>&1; then
    exec cargo "$@"
  fi
  bash "${REPO_ROOT}/tools/build_legacy_rnaview.sh"
fi

if [[ ! -x "${SYS_LINKER}" ]]; then
  echo "missing sysroot linker after setup: ${SYS_LINKER}" >&2
  exit 2
fi

export PATH="${SYSROOT}/usr/bin:${SYSROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${SYSROOT}/usr/lib/x86_64-linux-gnu:${SYSROOT}/lib/x86_64-linux-gnu:${SYSROOT}/usr/lib:${SYSROOT}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER="${SYS_LINKER}"
export RUSTFLAGS="-C link-arg=--sysroot=${SYSROOT}${RUSTFLAGS:+ ${RUSTFLAGS}}"

exec cargo "$@"
