#!/bin/bash
# The small, independent externals: tinyxml2, freetype, mmg, hdf5, gdcm.
# CMake flags follow the in-repo recipes
# (Externals/Make/2022.10/BuildHelpers/CompileScripts/*), versions follow the
# active 2024.05 block of Code/CMake/SimVascularExternalsVersions.cmake.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

# Each section skips if its install dir is already populated, so this script
# (and build-all.sh) can be safely re-run after a failure.

# ---- tinyxml2 8.0.0 ---------------------------------------------------------
banner "tinyxml2 $TINYXML2_VERSION"
if [ -d "$INSTALL_DIR/tinyxml2/lib" ]; then echo "already installed, skipping"; else
fetch "https://github.com/leethomason/tinyxml2/archive/refs/tags/$TINYXML2_VERSION.tar.gz" "tinyxml2-$TINYXML2_VERSION.tar.gz"
untar "tinyxml2-$TINYXML2_VERSION.tar.gz" "$BLD_DIR/tinyxml2"
cmake -S "$BLD_DIR/tinyxml2" -B "$BLD_DIR/tinyxml2/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/tinyxml2"
cmake --build "$BLD_DIR/tinyxml2/build" -j"$NPROC"
cmake --install "$BLD_DIR/tinyxml2/build"
fi

# ---- freetype 2.13.0 --------------------------------------------------------
banner "freetype $FREETYPE_VERSION"
if [ -d "$INSTALL_DIR/freetype/lib" ]; then echo "already installed, skipping"; else
fetch "https://download.savannah.gnu.org/releases/freetype/freetype-$FREETYPE_VERSION.tar.gz"
untar "freetype-$FREETYPE_VERSION.tar.gz" "$BLD_DIR/freetype"
cmake -S "$BLD_DIR/freetype" -B "$BLD_DIR/freetype/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/freetype"
cmake --build "$BLD_DIR/freetype/build" -j"$NPROC"
cmake --install "$BLD_DIR/freetype/build"
fi

# ---- mmg 5.3.9 (static, as registered: shared OFF) --------------------------
# -fcommon: mmg 5.3.9 predates GCC 10's -fno-common default; without it the
# executables fail to link with "multiple definition" errors on Ubuntu 22.04.
banner "mmg $MMG_VERSION"
if [ -d "$INSTALL_DIR/mmg/lib" ]; then echo "already installed, skipping"; else
fetch "https://github.com/MmgTools/mmg/archive/refs/tags/v$MMG_VERSION.tar.gz" "mmg-$MMG_VERSION.tar.gz"
untar "mmg-$MMG_VERSION.tar.gz" "$BLD_DIR/mmg"
cmake -S "$BLD_DIR/mmg" -B "$BLD_DIR/mmg/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DCMAKE_C_FLAGS="-fcommon" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/mmg"
cmake --build "$BLD_DIR/mmg/build" -j"$NPROC"
cmake --install "$BLD_DIR/mmg/build"
fi

# ---- hdf5 1.14.3 (C + C++ + HL, shared; consumed by GDCM/ITK/MITK) ----------
banner "hdf5 $HDF5_VERSION"
if [ -d "$INSTALL_DIR/hdf5/lib" ]; then echo "already installed, skipping"; else
fetch "https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-${HDF5_VERSION//./_}.tar.gz" "hdf5-$HDF5_VERSION.tar.gz"
untar "hdf5-$HDF5_VERSION.tar.gz" "$BLD_DIR/hdf5"
cmake -S "$BLD_DIR/hdf5" -B "$BLD_DIR/hdf5/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DHDF5_BUILD_CPP_LIB=ON \
  -DHDF5_BUILD_HL_LIB=ON \
  -DBUILD_TESTING=OFF \
  -DHDF5_BUILD_EXAMPLES=OFF \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/hdf5"
cmake --build "$BLD_DIR/hdf5/build" -j"$NPROC"
cmake --install "$BLD_DIR/hdf5/build"
fi

# ---- gdcm 3.0.10 (shared; ITK/MITK use it as system GDCM) -------------------
banner "gdcm $GDCM_VERSION"
if [ -d "$INSTALL_DIR/gdcm/lib" ]; then echo "already installed, skipping"; else
fetch "https://github.com/malaterre/GDCM/archive/refs/tags/v$GDCM_VERSION.tar.gz" "gdcm-$GDCM_VERSION.tar.gz"
untar "gdcm-$GDCM_VERSION.tar.gz" "$BLD_DIR/gdcm"
cmake -S "$BLD_DIR/gdcm" -B "$BLD_DIR/gdcm/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DGDCM_BUILD_SHARED_LIBS=ON \
  -DGDCM_BUILD_TESTING=OFF \
  -DGDCM_BUILD_EXAMPLES=OFF \
  -DGDCM_BUILD_DOCBOOK_MANPAGES=OFF \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/gdcm"
cmake --build "$BLD_DIR/gdcm/build" -j"$NPROC"
cmake --install "$BLD_DIR/gdcm/build"
fi

banner "small libs done"
