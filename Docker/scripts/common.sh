# Shared helpers for the externals build scripts. Source, do not execute.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/versions.env"

mkdir -p "$SRC_DIR" "$BLD_DIR" "$INSTALL_DIR"

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
