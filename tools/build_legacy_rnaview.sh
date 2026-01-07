#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if command -v make >/dev/null 2>&1 && command -v cc >/dev/null 2>&1; then
  make
  exit 0
fi

TOOLCHAIN_DIR="${HOME}/.cache/rnaview-toolchain"
SYSROOT="${TOOLCHAIN_DIR}/root"

if [[ -x "${SYSROOT}/usr/bin/make" && -x "${SYSROOT}/usr/bin/gcc-14" ]]; then
  export PATH="${SYSROOT}/usr/bin:${SYSROOT}/bin:${PATH}"
  export LD_LIBRARY_PATH="${SYSROOT}/usr/lib/x86_64-linux-gnu:${SYSROOT}/lib/x86_64-linux-gnu:${SYSROOT}/usr/lib:${SYSROOT}/lib"

  make \
    CC="${SYSROOT}/usr/bin/gcc-14" \
    CFLAGS="--sysroot=${SYSROOT} -isystem ${SYSROOT}/usr/include -isystem ${SYSROOT}/usr/include/x86_64-linux-gnu -Iinclude" \
    LDFLAGS="--sysroot=${SYSROOT}"

  echo "built: ${REPO_ROOT}/bin/rnaview"
  exit 0
fi

APTROOT="${HOME}/.apt"
DEBDIR="${TOOLCHAIN_DIR}/debs"

mkdir -p "${APTROOT}/lists/partial" "${APTROOT}/cache/archives/partial"
: > "${APTROOT}/status"

apt-get \
  -o Dir::State="${APTROOT}" \
  -o Dir::State::lists="${APTROOT}/lists" \
  -o Dir::State::status="${APTROOT}/status" \
  -o Dir::Cache="${APTROOT}/cache" \
  -o Dir::Cache::archives="${APTROOT}/cache/archives" \
  -o Debug::NoLocking=1 \
  update >/dev/null

rm -rf "${TOOLCHAIN_DIR}"
mkdir -p "${DEBDIR}" "${SYSROOT}"

pkgs=(
  make
  gcc-14 gcc-14-x86-64-linux-gnu gcc-14-base
  cpp-14 cpp-14-x86-64-linux-gnu
  libgcc-14-dev libcc1-0
  libgmp10 libisl23 libmpc3 libmpfr6
  zlib1g libzstd1
  libstdc++6 libgcc-s1
  libgomp1 libitm1 libatomic1 libasan8 liblsan0 libtsan2 libubsan1 libhwasan0 libquadmath0
  binutils-x86-64-linux-gnu binutils-common libbinutils libctf0 libctf-nobfd0 libgprofng0 libsframe1 libjansson4
  libc6 libc6-dev libc-dev-bin linux-libc-dev libcrypt-dev libcrypt1 rpcsvc-proto
)

cd "${DEBDIR}"
apt-get \
  -o Dir::State="${APTROOT}" \
  -o Dir::State::lists="${APTROOT}/lists" \
  -o Dir::State::status="${APTROOT}/status" \
  -o Dir::Cache="${APTROOT}/cache" \
  -o Dir::Cache::archives="${APTROOT}/cache/archives" \
  -o Debug::NoLocking=1 \
  download "${pkgs[@]}" >/dev/null

for deb in ./*.deb; do
  dpkg-deb -x "$deb" "${SYSROOT}"
done

# Provide merged-usr symlinks expected by linker scripts.
if [[ ! -e "${SYSROOT}/lib" ]]; then
  ln -s usr/lib "${SYSROOT}/lib"
fi
if [[ ! -e "${SYSROOT}/lib64" ]]; then
  ln -s usr/lib64 "${SYSROOT}/lib64"
fi

# Convenience symlinks expected by some toolchains.
ln -sf x86_64-linux-gnu-ld "${SYSROOT}/usr/bin/ld" || true
ln -sf x86_64-linux-gnu-as "${SYSROOT}/usr/bin/as" || true
ln -sf x86_64-linux-gnu-ar "${SYSROOT}/usr/bin/ar" || true
ln -sf x86_64-linux-gnu-ranlib "${SYSROOT}/usr/bin/ranlib" || true
ln -sf gcc-14 "${SYSROOT}/usr/bin/cc" || true

export PATH="${SYSROOT}/usr/bin:${SYSROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${SYSROOT}/usr/lib/x86_64-linux-gnu:${SYSROOT}/lib/x86_64-linux-gnu:${SYSROOT}/usr/lib:${SYSROOT}/lib"

cd "${REPO_ROOT}"
make \
  CC="${SYSROOT}/usr/bin/gcc-14" \
  CFLAGS="--sysroot=${SYSROOT} -isystem ${SYSROOT}/usr/include -isystem ${SYSROOT}/usr/include/x86_64-linux-gnu -Iinclude" \
  LDFLAGS="--sysroot=${SYSROOT}"

echo "built: ${REPO_ROOT}/bin/rnaview"
