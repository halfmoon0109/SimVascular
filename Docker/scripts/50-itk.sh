#!/bin/bash
# ITK 5.4.0 -- shared, system GDCM/HDF5, VtkGlue + Review modules.
# Flags follow the in-repo recipe (compile-cmake-itk-generic.sh) 1:1.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
banner "ITK $ITK_VERSION"

fetch "https://github.com/InsightSoftwareConsortium/ITK/releases/download/v$ITK_VERSION/InsightToolkit-$ITK_VERSION.tar.gz"
untar "InsightToolkit-$ITK_VERSION.tar.gz" "$BLD_DIR/itk"

VTK_MM=${VTK_VERSION%.*}
GDCM_MM=${GDCM_VERSION%.*}

cmake -S "$BLD_DIR/itk" -B "$BLD_DIR/itk/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DBUILD_TESTING=OFF \
  -DBUILD_EXAMPLES=OFF \
  -DModule_ITKReview=ON \
  -DModule_ITKVtkGlue=ON \
  -DITK_USE_SYSTEM_GDCM=ON \
  -DGDCM_DIR="$INSTALL_DIR/gdcm/lib/gdcm-$GDCM_MM" \
  -DITK_USE_SYSTEM_HDF5=ON \
  -DHDF5_DIR="$INSTALL_DIR/hdf5/cmake" \
  -DVTK_DIR="$INSTALL_DIR/vtk/lib/cmake/vtk-$VTK_MM" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/itk"
cmake --build "$BLD_DIR/itk/build" -j"$NPROC"
cmake --install "$BLD_DIR/itk/build"

echo "OK: ITK installed to $INSTALL_DIR/itk"
