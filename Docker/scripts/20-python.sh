#!/bin/bash
# CPython 3.11 -- shared build into install/python.
# run-cmake.sh expects install/python/{bin/python3, lib, include/python3.11}.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
banner "Python $PYTHON_VERSION"

fetch "https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz"
untar "Python-$PYTHON_VERSION.tgz" "$BLD_DIR/python"

cd "$BLD_DIR/python"
./configure --prefix="$INSTALL_DIR/python" \
            --enable-shared \
            --with-ensurepip=install \
            LDFLAGS="-Wl,-rpath,$INSTALL_DIR/python/lib"
make -j"$NPROC"
make install

# numpy: needed at runtime by the sv Python modules (NUMPY external in
# SimVascularOptions is delivered via pip in every externals generation).
LD_LIBRARY_PATH="$INSTALL_DIR/python/lib" \
  "$INSTALL_DIR/python/bin/python3" -m pip install --upgrade pip numpy

"$INSTALL_DIR/python/bin/python3" --version
