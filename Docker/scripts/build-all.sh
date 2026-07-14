#!/bin/bash
# Build all externals in dependency order, then SimVascular.
# Each stage logs to $WORK/logs/<stage>.log and STOPS on first failure --
# per project policy, a failing external is diagnosed, never skipped.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/versions.env"
mkdir -p "$WORK/logs"

STAGES=(
  10-qt6.sh
  20-python.sh
  30-small-libs.sh
  40-vtk.sh
  50-itk.sh
  60-opencascade.sh
  70-mitk.sh
  80-simvascular.sh
)

for stage in "${STAGES[@]}"; do
  log="$WORK/logs/${stage%.sh}.log"
  echo ">>> $stage (log: $log)"
  if ! bash "$SCRIPT_DIR/$stage" 2>&1 | tee "$log"; then
    echo
    echo "!!! FAILED at $stage -- see $log"
    echo "!!! Do not skip this component; report the tail of the log."
    exit 1
  fi
done

echo "All stages completed."
