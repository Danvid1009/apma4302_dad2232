#!/usr/bin/env bash
# HW4 Q3: temperature-driven RHS biharmonic; VTK for ParaView.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=inc_firedrake_apptainer.sh
source "${SCRIPT_DIR}/inc_firedrake_apptainer.sh"

HW4_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${HW4_ROOT}/python"
export HW4_RESULT_DIR="${HW4_RESULT_DIR:-${HW4_ROOT}/output/q3}"

echo "Running biharm_temperature_rhs.py -> ${HW4_RESULT_DIR}"
if [[ -n "${HW4_APPTAINER_SIF:-}${FIREDRAKE_TS_SIF:-}" ]]; then
  echo "(using Apptainer/Singularity image for python3)"
fi
hw4_python3 biharm_temperature_rhs.py "$@"

echo "Open ParaView: ${HW4_RESULT_DIR}/biharm_T_rhs.pvd"
