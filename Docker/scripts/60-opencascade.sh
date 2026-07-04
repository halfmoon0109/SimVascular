#!/bin/bash
# OpenCASCADE 7.6.0 -- with VTK integration (IVtk), no Draw (so no Tcl/Tk).
# Follows the in-repo recipe (compile-cmake-opencascade-generic.sh): USE_VTK=ON,
# BUILD_MODULE_Draw=OFF, freetype from our externals. The Qt5WebKit lines in
# the old recipe belonged to the Draw/Inspector module which is OFF here.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
banner "OpenCASCADE $OCC_VERSION"

OCC_TAG="V${OCC_VERSION//./_}"
fetch "https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/$OCC_TAG.tar.gz" "occt-$OCC_VERSION.tar.gz"
untar "occt-$OCC_VERSION.tar.gz" "$BLD_DIR/occt"

VTK_MM=${VTK_VERSION%.*}

cmake -S "$BLD_DIR/occt" -B "$BLD_DIR/occt/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_LIBRARY_TYPE=Shared \
  -DBUILD_MODULE_Draw=OFF \
  -DUSE_VTK=ON \
  -D3RDPARTY_VTK_DIR="$INSTALL_DIR/vtk" \
  -D3RDPARTY_VTK_INCLUDE_DIR="$INSTALL_DIR/vtk/include/vtk-$VTK_MM" \
  -D3RDPARTY_VTK_LIBRARY_DIR="$INSTALL_DIR/vtk/lib" \
  -DUSE_FREETYPE=ON \
  -D3RDPARTY_FREETYPE_DIR="$INSTALL_DIR/freetype" \
  -D3RDPARTY_FREETYPE_INCLUDE_DIR_freetype2="$INSTALL_DIR/freetype/include/freetype2" \
  -D3RDPARTY_FREETYPE_INCLUDE_DIR_ft2build="$INSTALL_DIR/freetype/include/freetype2" \
  -D3RDPARTY_FREETYPE_LIBRARY_DIR="$INSTALL_DIR/freetype/lib" \
  -DINSTALL_DIR="$INSTALL_DIR/opencascade" \
  -DINSTALL_DIR_LIB=lib \
  -DINSTALL_DIR_CMAKE=lib/cmake/opencascade
cmake --build "$BLD_DIR/occt/build" -j"$NPROC"
cmake --install "$BLD_DIR/occt/build"

test -d "$INSTALL_DIR/opencascade/lib/cmake/opencascade" \
  && echo "OK: layout expected by Code/CMake/Externals/OpenCASCADE.cmake:49"
