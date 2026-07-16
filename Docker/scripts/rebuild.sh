#!/bin/bash
# Incremental rebuild after editing SimVascular source only (no externals
# change). Unlike build-all.sh (full externals+SimVascular pipeline, run
# rarely), this is the fast day-to-day loop from WORKFLOW.md section 4.1 --
# just `make` in the existing build tree -- wrapped so it also lands a log
# under $WORK/logs/, matching the per-stage logs build-all.sh produces.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

BUILD="$WORK/build"
[ -d "$BUILD/SimVascular-build" ] && BUILD="$BUILD/SimVascular-build"

[ -f "$BUILD/Makefile" ] || { echo "ERROR: no build tree at $BUILD (run build-all.sh first)"; exit 1; }

mkdir -p "$WORK/logs"
banner "Incremental rebuild: $BUILD"

cd "$BUILD"
make -j"$NPROC" 2>&1 | tee "$WORK/logs/80-simvascular.log"
