#!/bin/bash
# Launch the built SimVascular GUI with all runtime environment set up.
# Usage (inside the container, with a display connected -- see Docker/README.md):
#   bash /work/SimVascular/Docker/scripts/run-gui.sh
#
# This resolves the run-time shared libraries of every external (GDCM, VTK,
# ITK, OpenCASCADE, MITK, Qt, ...), points Qt at its plugins (platform + the
# sqldrivers plugin CTK's plugin database needs), and forces software OpenGL
# which WSLg/headless X sessions require.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/versions.env"
E="$INSTALL_DIR"

export LD_LIBRARY_PATH="\
$E/qt6/lib:$E/vtk/lib:$E/itk/lib:$E/gdcm/lib:$E/hdf5/lib:$E/opencascade/lib:\
$E/tinyxml2/lib:$E/freetype/lib:$E/mmg/lib:$E/python/lib:$E/mitk/lib:\
$E/mitk/lib/plugins:${LD_LIBRARY_PATH:-}"

# CTK/BlueBerry plugin search paths (sv4gui_MitkApp reads SV_PLUGIN_PATH and
# feeds each entry to ctkPluginFrameworkLauncher::addSearchPath). Must include
# BOTH SimVascular's own plugins (build tree) and MITK/BlueBerry/CTK plugins
# (our assembled MITK install), else org.mitk.gui.common fails to resolve.
export SV_PLUGIN_PATH="\
$WORK/build/SimVascular-build/lib/plugins:\
$E/mitk/lib/plugins:\
$E/mitk/bin/plugins"

export QT_PLUGIN_PATH="$E/qt6/plugins"
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-xcb}"
export QT_QPA_PLATFORM_PLUGIN_PATH="$E/qt6/plugins/platforms"
export LIBGL_ALWAYS_SOFTWARE="${LIBGL_ALWAYS_SOFTWARE:-1}"
export LC_ALL="${LC_ALL:-C.UTF-8}"
export LANG="${LANG:-C.UTF-8}"

BIN="$WORK/build/SimVascular-build/bin/simvascular"
echo "Launching $BIN"
exec "$BIN" "$@"
