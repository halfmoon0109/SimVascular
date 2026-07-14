#!/bin/bash
# Qt 6.6.2 -- official binary release via aqtinstall (Qt mirrors).
#
# The old in-repo superbuild also used prebuilt Qt binaries (the repo notes
# state "Building Qt from source doesn't work!"); the only change here is the
# download source, since the original host is gone.
#
# Modules (aqt -m) that are NOT in the base gcc_64 desktop install but are
# required by consumers:
#   qt5compat     -> Core5Compat   (Code/CMake/Externals/Qt6.cmake)
#   qtscxml       -> StateMachine  (MITK v2024.06 MITK_QT6_COMPONENTS)
#   qtwebengine + qtwebchannel + qtpositioning + qtwebview -> WebEngine*/WebView
#   qtshadertools -> dependency of the above
# Everything else MITK/SimVascular ask for (Core, Gui, Widgets, Network, Qml,
# Svg, OpenGL, OpenGLWidgets, Designer, Help, UiTools, LinguistTools, Sql,
# Concurrent, Xml, PrintSupport) ships in the base desktop install.
QT_MODULES="qt5compat qtscxml qtwebengine qtwebview qtpositioning qtwebchannel qtshadertools"

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
banner "Qt $QT_VERSION (aqtinstall)"

QT_ROOT="$INSTALL_DIR/qt6-dist"
QT_GCC="$QT_ROOT/$QT_VERSION/gcc_64"
# Key the resume-guard on a module-provided config (StateMachine), not just the
# base Qt6 dir, so a Qt tree missing the extra modules gets reinstalled.
if [ -f "$QT_GCC/lib/cmake/Qt6StateMachine/Qt6StateMachineConfig.cmake" ]; then
  echo "Qt already installed with required modules."
else
  rm -rf "$QT_ROOT"
  aqt install-qt linux desktop "$QT_VERSION" gcc_64 -m $QT_MODULES -O "$QT_ROOT"
fi

# run-cmake.sh expects install/qt6/lib/cmake/Qt6
ln -sfn "$QT_ROOT/$QT_VERSION/gcc_64" "$INSTALL_DIR/qt6"
echo "Qt6 at $INSTALL_DIR/qt6 -> $(readlink -f "$INSTALL_DIR/qt6")"
