#!/bin/bash
# MITK -- the one external with a deliberate version decision (see
# versions.env). Superbuild against OUR Qt6/VTK/ITK/GDCM/Python, mirroring the
# in-repo recipe (compile-cmake-mitk-linux.sh): superbuild ON, external
# ITK/VTK dirs, GDCM on, SimpleITK off, examples/testing off.
#
# NOTE: the MITK superbuild downloads its remaining third-party deps (CTK,
# DCMTK, POCO, ...) during the build -- network access is required.
#
# This is the longest and riskiest step. If it fails, capture the full error
# and stop; possible documented fallbacks (choose deliberately, not silently):
#   1) Let MITK build its own patched VTK/ITK (drop EXTERNAL_*_DIR), then point
#      SimVascular's SV_VTK_DIR/SV_ITK_DIR at MITK's ep/ trees instead.
#   2) Try MITK_GIT_TAG=v2024.06 (newer API) or v2023.04 (older Qt6 API),
#      depending on which sv4gui API mismatch appears.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
banner "MITK $MITK_GIT_TAG"

if [ -d "$INSTALL_DIR/mitk/lib" ]; then
  echo "already installed, skipping"; exit 0
fi

if [ ! -d "$BLD_DIR/mitk/src/.git" ]; then
  rm -rf "$BLD_DIR/mitk/src"
  git clone --branch "$MITK_GIT_TAG" --depth 1 https://github.com/MITK/MITK.git "$BLD_DIR/mitk/src"
fi

VTK_MM=${VTK_VERSION%.*}
GDCM_MM=${GDCM_VERSION%.*}
PY_MM=${PYTHON_VERSION%.*}
PY="$INSTALL_DIR/python/bin/python3"

# Same VTK-transitively-re-finds-Python3 issue as in 50-itk.sh (EXTERNAL_VTK_DIR
# makes MITK plain find_package(VTK) it too): pin Python3 explicitly so it
# doesn't resolve to the system 3.10 interpreter.
cmake -S "$BLD_DIR/mitk/src" -B "$BLD_DIR/mitk/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DBUILD_TESTING=OFF \
  -DMITK_BUILD_EXAMPLES=OFF \
  -DQt6_DIR="$INSTALL_DIR/qt6/lib/cmake/Qt6" \
  -DCMAKE_PREFIX_PATH="$INSTALL_DIR/qt6" \
  -DEXTERNAL_VTK_DIR="$INSTALL_DIR/vtk/lib/cmake/vtk-$VTK_MM" \
  -DEXTERNAL_ITK_DIR="$INSTALL_DIR/itk/lib/cmake/ITK-${ITK_VERSION%.*}" \
  -DMITK_USE_GDCM=ON \
  -DGDCM_DIR="$INSTALL_DIR/gdcm/lib/gdcm-$GDCM_MM" \
  -DMITK_USE_HDF5=ON \
  -DHDF5_DIR="$INSTALL_DIR/hdf5/cmake" \
  -DMITK_USE_Python3=ON \
  -DPython3_EXECUTABLE="$PY" \
  -DPython3_ROOT_DIR="$INSTALL_DIR/python" \
  -DPython3_FIND_STRATEGY=LOCATION \
  -DPython3_INCLUDE_DIR="$INSTALL_DIR/python/include/python$PY_MM" \
  -DPython3_LIBRARY="$INSTALL_DIR/python/lib/libpython$PY_MM.so" \
  -DMITK_USE_SimpleITK=OFF \
  -DMITK_USE_BLUEBERRY=ON
cmake --build "$BLD_DIR/mitk/build" -j"$NPROC"

# Install the inner MITK-build (superbuild wrapper has no install target) into
# the layout FindMITK.cmake expects: a toplevel dir with bin/lib/include.
cmake --install "$BLD_DIR/mitk/build/MITK-build" --prefix "$INSTALL_DIR/mitk"

echo "OK: MITK installed to $INSTALL_DIR/mitk"
