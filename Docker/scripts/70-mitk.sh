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

# Marker of a completed SimVascular-layout assembly (created by the source
# header mirror near the end of this script).
if [ -d "$INSTALL_DIR/mitk/include/mitk/Modules/Core" ]; then
  echo "already installed, skipping"; exit 0
fi

# Check for CMakeLists.txt, not just .git -- a clone that fails partway
# (e.g. bad tag name) still leaves an initialized .git dir behind.
if [ ! -f "$BLD_DIR/mitk/src/CMakeLists.txt" ]; then
  rm -rf "$BLD_DIR/mitk/src"
  git clone --branch "$MITK_GIT_TAG" --depth 1 https://github.com/MITK/MITK.git "$BLD_DIR/mitk/src"
fi

VTK_MM=${VTK_VERSION%.*}
GDCM_MM=${GDCM_VERSION%.*}
PY_MM=${PYTHON_VERSION%.*}
PY="$INSTALL_DIR/python/bin/python3"

# Compile the MITK superbuild only if it hasn't already produced its core lib
# (so a re-run after fixing only the install/assembly step is fast).
if [ ! -f "$BLD_DIR/mitk/build/MITK-build/lib/libMitkCore.so" ]; then
# Same VTK-transitively-re-finds-Python3 issue as in 50-itk.sh (EXTERNAL_VTK_DIR
# makes MITK plain find_package(VTK) it too): pin Python3 explicitly so it
# doesn't resolve to the system 3.10 interpreter.
cmake -S "$BLD_DIR/mitk/src" -B "$BLD_DIR/mitk/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON \
  -DBUILD_TESTING=OFF \
  -DMITK_BUILD_EXAMPLES=OFF \
  -DMITK_USE_Qt6=ON \
  -DMITK_USE_Qt5=OFF \
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
fi   # end compile guard

# ---------------------------------------------------------------------------
# Assemble the MITK install tree in the SimVascular-specific layout that
# Code/CMake/FindMITK.cmake expects. This is NOT a plain `cmake --install`:
# SimVascular searches ${MITK_DIR}/lib for MITK + its ep deps (DCMTK, Poco,
# CTK) and ${MITK_DIR}/include/mitk/<Module>/... assembled from the MITK
# SOURCE tree plus build-tree generated headers. Ported from the repo's own
# Externals/Make/2022.10/BuildHelpers/CompileScripts/post-install-mitk-linux.sh
# ---------------------------------------------------------------------------
banner "Assembling MITK install tree ($INSTALL_DIR/mitk)"
SRC="$BLD_DIR/mitk/src"
BLD="$BLD_DIR/mitk/build"
DST="$INSTALL_DIR/mitk"
rm -rf "$DST"
mkdir -p "$DST"/{bin,lib,lib/plugins,include,share} \
         "$DST"/include/mitk/{configs,exports,ui_files,Modules} \
         "$DST"/include/ctk

# ---- libraries + binaries: MITK-build + ep (DCMTK/Poco/CppMicroServices) ----
cp -Rf "$BLD/MITK-build/bin/."  "$DST/bin/"           2>/dev/null || true
cp -Rf "$BLD/MITK-build/lib/."  "$DST/lib/"           2>/dev/null || true
cp -Rf "$BLD/ep/bin/."          "$DST/bin/"           2>/dev/null || true
cp -Rf "$BLD/ep/lib/."          "$DST/lib/"           2>/dev/null || true
cp -Rf "$BLD/ep/include/."      "$DST/include/"       2>/dev/null || true
cp -Rf "$BLD/ep/share/."        "$DST/share/"         2>/dev/null || true
cp -Rf "$BLD/ep/plugins/."      "$DST/lib/plugins/"   2>/dev/null || true

# ---- CTK / qRestAPI / PythonQt shared libs (live under CTK-build) ----
CTKBIN="$BLD/ep/src/CTK-build/CTK-build/bin"
cp -f "$CTKBIN"/*CTK*.so*   "$DST/lib/"         2>/dev/null || true
cp -f "$CTKBIN"/liborg*.so* "$DST/lib/plugins/" 2>/dev/null || true
cp -f "$BLD"/ep/src/CTK-build/qRestAPI-build/*qRestAPI*.so* "$DST/lib/" 2>/dev/null || true
cp -f "$BLD"/ep/src/CTK-build/PythonQt-build/*PythonQt*.so* "$DST/lib/" 2>/dev/null || true

# ---- ctk headers ----
for d in Core Scripting/Python/Core Scripting/Python/Widgets \
         Visualization/VTK/Core Widgets; do
  cp -f "$BLD"/ep/src/CTK/Libs/$d/*.h "$DST/include/ctk/" 2>/dev/null || true
done
find "$BLD/ep/src/CTK-build" -name "*Export.h" -exec cp -f {} "$DST/include/ctk/" \; 2>/dev/null || true
[ -d "$BLD/ep/src/CTK/Libs/PluginFramework" ] && \
  cp -Rf "$BLD/ep/src/CTK/Libs/PluginFramework" "$DST/include/ctk/" 2>/dev/null || true

# ---- top-level MITK-build generated headers ----
cp -f "$BLD"/MITK-build/*.h "$DST/include/mitk/" 2>/dev/null || true

# ---- all source headers from Modules + Plugins, structure preserved ----
# FindMITK.cmake globs ${MITK_DIR}/include up to 6 levels deep, so exact
# placement is irrelevant as long as each header exists somewhere under
# include/mitk within that depth. Mirroring the source tree (headers only)
# is the robust way to satisfy its fixed header list without hand-mapping
# each module's internal subdirectory layout (which changed across MITK
# versions, e.g. Modules/ContourModel/DataManagement).
( cd "$SRC" && find Modules Plugins \
    \( -name '*.h' -o -name '*.hpp' -o -name '*.hxx' -o -name '*.tpp' \) \
    -exec cp -f --parents {} "$DST/include/mitk/" \; ) 2>/dev/null || true
find "$BLD/MITK-build/Plugins" -name "*Export.h" -exec cp -f {} "$DST/include/mitk/exports/" \; 2>/dev/null || true

# ---- build-tree generated headers: exports / configs / ui ----
find "$BLD"                       -name "*Exports.h" -exec cp -f {} "$DST/include/mitk/exports/"  \; 2>/dev/null || true
find "$BLD/MITK-build/Modules"    -name "*Export.h"  -exec cp -f {} "$DST/include/mitk/exports/"  \; 2>/dev/null || true
find "$BLD/MITK-build/Modules"    -name "ui_*.h"     -exec cp -f {} "$DST/include/mitk/ui_files/" \; 2>/dev/null || true
find "$BLD/MITK-build"            -name "*Config.h"  -exec cp -f {} "$DST/include/mitk/configs/"   \; 2>/dev/null || true

echo "OK: MITK assembled at $DST"
echo "    libs:    $(find "$DST/lib" -maxdepth 1 -name '*.so*' | wc -l) shared objects"
echo "    plugins: $(find "$DST/lib/plugins" -maxdepth 1 -name '*.so*' | wc -l) plugin libs"
