# Shared helpers for the externals build scripts. Source, do not execute.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/versions.env"

mkdir -p "$SRC_DIR" "$BLD_DIR" "$INSTALL_DIR"

# Put our Python 3.11 ahead of the system 3.10 for every stage and every
# sub-process. VTK was built with VTK_WRAP_PYTHON=ON, so its exported config
# transitively runs FindPython3 (>=3.11) in any consumer that find_package(VTK)s
# it -- including MITK's own ExternalProject sub-builds (ACVD, etc.) which we
# cannot pass -DPython3_* hints into individually. With python3.11 first on
# PATH, FindPython3 resolves to 3.11 without hints. LD_LIBRARY_PATH is needed
# because our interpreter is built --enable-shared.
if [ -d "$INSTALL_DIR/python/bin" ]; then
  export PATH="$INSTALL_DIR/python/bin:$PATH"
  export LD_LIBRARY_PATH="$INSTALL_DIR/python/lib:${LD_LIBRARY_PATH:-}"
fi

# fetch <url> [output-name] : download into $SRC_DIR unless already present.
fetch() {
  local url="$1"
  local out="${2:-$(basename "$url")}"
  if [ -f "$SRC_DIR/$out" ]; then
    echo "[fetch] $out already present, skipping"
  else
    echo "[fetch] $url"
    wget -q --show-progress -O "$SRC_DIR/$out.part" "$url"
    mv "$SRC_DIR/$out.part" "$SRC_DIR/$out"
  fi
}

# untar <tarball> <dest-dir> : extract fresh (removes any previous tree).
untar() {
  local tarball="$1" dest="$2"
  rm -rf "$dest"
  mkdir -p "$dest"
  tar -xf "$SRC_DIR/$tarball" -C "$dest" --strip-components=1
}

banner() {
  echo
  echo "============================================================"
  echo "  $*"
  echo "============================================================"
}
