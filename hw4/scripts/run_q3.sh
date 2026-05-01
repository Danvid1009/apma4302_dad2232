#!/usr/bin/env bash
# HW4 Q3: temperature-driven RHS biharmonic; VTK for ParaView.
set -euo pipefail

HW4_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${HW4_ROOT}/python"
export HW4_RESULT_DIR="${HW4_RESULT_DIR:-${HW4_ROOT}/output/q3}"

echo "Running biharm_temperature_rhs.py -> ${HW4_RESULT_DIR}"
if command -v mpiexec >/dev/null 2>&1; then
  mpiexec -n 1 python3 biharm_temperature_rhs.py "$@"
else
  python3 biharm_temperature_rhs.py "$@"
fi

echo "Open ParaView: ${HW4_RESULT_DIR}/biharm_T_rhs.pvd"
