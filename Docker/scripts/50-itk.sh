#!/bin/bash
# ITK 5.4.0 -- shared, system GDCM/HDF5, VtkGlue + Review modules.
# Flags follow the in-repo recipe (compile-cmake-itk-generic.sh) 1:1.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
banner "ITK $ITK_VERSION"

if [ -d "$INSTALL_DIR/itk/lib/cmake/ITK-${ITK_VERSION%.*}" ]; then
  echo "already installed, skipping"; exit 0
fi

fetch "https://github.com/InsightSoftwareConsortium/ITK/releases/download/v$ITK_VERSION/InsightToolkit-$ITK_VERSION.tar.gz"
untar "InsightToolkit-$ITK_VERSION.tar.gz" "$BLD_DIR/itk"

VTK_MM=${VTK_VERSION%.*}
GDCM_MM=${GDCM_VERSION%.*}
PY_MM=${PYTHON_VERSION%.*}

# VTK was built with VTK_WRAP_PYTHON=ON, so its exported config transitively
# find_package(Python3)s for any consumer (here: ITK's VtkGlue module). Without
# hints CMake's FindPython3 defaults to the system interpreter (3.10 on
# Ubuntu 22.04), which is older than the 3.11 VTK was built against -- pin it
# to our own Python explicitly.
cmake -S "$BLD_DIR/itk" -B "$BLD_DIR/itk/build" -G Ninja \
  -DPython3_EXECUTABLE="$INSTALL_DIR/python/bin/python3" \
  -DPython3_ROOT_DIR="$INSTALL_DIR/python" \
  -DPython3_FIND_STRATEGY=LOCATION \
  -DPython3_INCLUDE_DIR="$INSTALL_DIR/python/include/python$PY_MM" \
  -DPython3_LIBRARY="$INSTALL_DIR/python/lib/libpython$PY_MM.so" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DBUILD_TESTING=OFF \
  -DBUILD_EXAMPLES=OFF \
  -DModule_ITKReview=ON \
  -DModule_ITKVtkGlue=ON \
  -DModule_ITKOpenJPEG=ON \
  -DModule_GrowCut=ON \
  -DModule_IsotropicWavelets=ON \
  -DITK_SKIP_PATH_LENGTH_CHECKS=ON \
  -DITK_USE_SYSTEM_GDCM=ON \
  -DGDCM_DIR="$INSTALL_DIR/gdcm/lib/gdcm-$GDCM_MM" \
  -DITK_USE_SYSTEM_HDF5=ON \
  -DHDF5_DIR="$INSTALL_DIR/hdf5/cmake" \
  -DVTK_DIR="$INSTALL_DIR/vtk/lib/cmake/vtk-$VTK_MM" \
  -DQt6_DIR="$INSTALL_DIR/qt6/lib/cmake/Qt6" \
  -DCMAKE_PREFIX_PATH="$INSTALL_DIR/qt6" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/itk"
cmake --build "$BLD_DIR/itk/build" -j"$NPROC"
cmake --install "$BLD_DIR/itk/build"

echo "OK: ITK installed to $INSTALL_DIR/itk"
