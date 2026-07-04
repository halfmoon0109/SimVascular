#!/bin/bash
# Stage 2: build SimVascular (this repo) against the externals built by the
# earlier scripts. This mirrors the repo's own run-cmake.sh (same -DSV_*_DIR
# arguments and install/<pkg> layout), driven from /work.
#
# Usage: 80-simvascular.sh [path-to-SimVascular-repo]  (default: /work/SimVascular)
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

REPO="${1:-$WORK/SimVascular}"
BUILD="$WORK/build"
banner "SimVascular from $REPO"

[ -f "$REPO/CMakeLists.txt" ] || { echo "ERROR: no SimVascular repo at $REPO"; exit 1; }

# ML weights dir: run-cmake.sh requires --ml-data-path. The weights are only
# needed at runtime by the ML segmentation tool; an empty dir lets the build
# proceed.
mkdir -p "$INSTALL_DIR/ml"

mkdir -p "$BUILD"
cd "$BUILD"

export CC=/usr/bin/gcc CXX=/usr/bin/g++
# The build runs the externals Python interpreter; it needs its shared lib.
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$INSTALL_DIR/python/lib"

cmake "$REPO" \
  -DSV_ENABLE_DISTRIBUTION:BOOL=ON \
  -DSV_EXTERNALS_DIR:PATH="$INSTALL_DIR" \
  -DSV_FREETYPE_DIR:PATH="$INSTALL_DIR/freetype" \
  -DSV_GDCM_DIR:PATH="$INSTALL_DIR/gdcm" \
  -DSV_HDF5_DIR:PATH="$INSTALL_DIR/hdf5" \
  -DSV_ITK_DIR:PATH="$INSTALL_DIR/itk" \
  -DSV_MITK_DIR:PATH="$INSTALL_DIR/mitk" \
  -DSV_ML_DIR:PATH="$INSTALL_DIR/ml" \
  -DSV_MMG_DIR:PATH="$INSTALL_DIR/mmg" \
  -DSV_OpenCASCADE_DIR:PATH="$INSTALL_DIR/opencascade" \
  -DSV_Qt6_DIR:PATH="$INSTALL_DIR/qt6" \
  -DSV_PYTHON_DIR:PATH="$INSTALL_DIR/python" \
  -DSV_TINYXML2_DIR:PATH="$INSTALL_DIR/tinyxml2" \
  -DSV_VTK_DIR:PATH="$INSTALL_DIR/vtk" \
  -DPython3_INCLUDE_DIRS:PATH="$INSTALL_DIR/python/include/python3.11" \
  -DPython3_LIBRARY_DIRS:PATH="$INSTALL_DIR/python/lib" \
  -DPython_SV_EXTERNALS_DIR:PATH="$INSTALL_DIR/python" \
  -DPython_DIR:PATH="$INSTALL_DIR/python" \
  -DPYTHON_DIR:PATH="$INSTALL_DIR/python" \
  -DPYTHON_EXECUTABLE:PATH="$INSTALL_DIR/python/bin/python3" \
  -DPython3_EXECUTABLE:PATH="$INSTALL_DIR/python/bin/python3" \
  -DCMAKE_PREFIX_PATH:PATH="$INSTALL_DIR/python"

make -j"$NPROC"

banner "SimVascular build finished -- binaries under $BUILD/SimVascular-build"
