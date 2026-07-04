#!/bin/bash
# VTK 9.3.0 -- shared, Qt6 support, Python 3.11 wrapping.
#
# Flag mapping from the in-repo 2022.10 recipe (compile-cmake-vtk-generic.sh)
# to VTK 9.3 module syntax:
#   VTK_Group_Qt=1        -> VTK_GROUP_ENABLE_Qt=YES, VTK_QT_VERSION=6
#   VTK_Group_Imaging=1   -> VTK_GROUP_ENABLE_Imaging=YES
#   VTK_WRAP_PYTHON=1     -> VTK_WRAP_PYTHON=ON (against our Python 3.11)
#   VTK_WRAP_TCL          -> dropped: Code/CMake/Externals/VTK.cmake has the
#                            Tcl components commented out (no TCL external is
#                            registered for the Code build).
#   Module_vtkTestingRendering=1 -> VTK_MODULE_ENABLE_VTK_TestingRendering=YES
# Components required by Code/CMake/Externals/VTK.cmake (FiltersFlowPaths,
# RenderingLabel, ChartsCore, FiltersVerdict, IOPLY, ImagingStencil, tiff, ...)
# are all part of the default full build.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
banner "VTK $VTK_VERSION"

VTK_MM=${VTK_VERSION%.*}
if [ -d "$INSTALL_DIR/vtk/lib/cmake/vtk-$VTK_MM" ]; then
  echo "already installed, skipping"; exit 0
fi
fetch "https://www.vtk.org/files/release/$VTK_MM/VTK-$VTK_VERSION.tar.gz"
untar "VTK-$VTK_VERSION.tar.gz" "$BLD_DIR/vtk"

PY="$INSTALL_DIR/python/bin/python3"

cmake -S "$BLD_DIR/vtk" -B "$BLD_DIR/vtk/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DVTK_BUILD_TESTING=OFF \
  -DVTK_BUILD_EXAMPLES=OFF \
  -DVTK_GROUP_ENABLE_Qt=YES \
  -DVTK_GROUP_ENABLE_Imaging=YES \
  -DVTK_QT_VERSION=6 \
  -DQt6_DIR="$INSTALL_DIR/qt6/lib/cmake/Qt6" \
  -DCMAKE_PREFIX_PATH="$INSTALL_DIR/qt6" \
  -DVTK_MODULE_ENABLE_VTK_GUISupportQt=YES \
  -DVTK_MODULE_ENABLE_VTK_RenderingQt=YES \
  -DVTK_MODULE_ENABLE_VTK_TestingRendering=YES \
  -DVTK_WRAP_PYTHON=ON \
  -DPython3_EXECUTABLE="$PY" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/vtk"
cmake --build "$BLD_DIR/vtk/build" -j"$NPROC"
cmake --install "$BLD_DIR/vtk/build"

test -d "$INSTALL_DIR/vtk/lib/cmake/vtk-$VTK_MM" \
  && echo "OK: $INSTALL_DIR/vtk/lib/cmake/vtk-$VTK_MM (layout expected by Code/CMake/Externals/VTK.cmake:50)"
