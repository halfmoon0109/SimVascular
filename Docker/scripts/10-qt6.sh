#!/bin/bash
# Qt 6.6.2 -- official binary release via aqtinstall (Qt mirrors).
#
# The old in-repo superbuild also used prebuilt Qt binaries (the repo notes
# state "Building Qt from source doesn't work!"); the only change here is the
# download source, since the original host is gone.
#
# Modules: qt5compat provides Core5Compat; qtwebengine/qtwebview (+their deps
# qtpositioning/qtwebchannel) provide WebEngineCore/WebEngineWidgets/WebView.
# Everything else required by Code/CMake/Externals/Qt6.cmake (Core, Gui,
# Widgets, Concurrent, Designer, Help, OpenGL, PrintSupport, Sql, Svg, Xml,
# UiTools) ships in the base desktop install (qtbase/qttools/qtsvg).
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
banner "Qt $QT_VERSION (aqtinstall)"

QT_ROOT="$INSTALL_DIR/qt6-dist"
if [ -d "$QT_ROOT/$QT_VERSION/gcc_64/lib/cmake/Qt6" ]; then
  echo "Qt already installed."
else
  aqt install-qt linux desktop "$QT_VERSION" gcc_64 \
      -m qt5compat qtwebengine qtwebview qtpositioning qtwebchannel qtshadertools \
      -O "$QT_ROOT"
fi

# run-cmake.sh expects install/qt6/lib/cmake/Qt6
ln -sfn "$QT_ROOT/$QT_VERSION/gcc_64" "$INSTALL_DIR/qt6"
echo "Qt6 at $INSTALL_DIR/qt6 -> $(readlink -f "$INSTALL_DIR/qt6")"
