#!/usr/bin/env bash
# HW4 Q4: convection DAE (firedrake-ts). Override RA, N, TMAX, DT via env vars.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=inc_firedrake_apptainer.sh
source "${SCRIPT_DIR}/inc_firedrake_apptainer.sh"

HW4_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${HW4_ROOT}/python"

RA="${RA:-1e2}"
N="${N:-64}"
TMAX="${TMAX:-1e5}"
DT="${DT:-0.1}"
OUT="${HW4_OUT:-${HW4_ROOT}/output/q4_ra${RA}_n${N}}"

echo "Running convection.py Ra=${RA} N=${N} t_max=${TMAX} dt=${DT} -> ${OUT}"
if [[ -n "${HW4_APPTAINER_SIF:-}${FIREDRAKE_TS_SIF:-}" ]]; then
  echo "(using Apptainer/Singularity image for python3)"
fi
hw4_python3 convection.py \
  --ra "${RA}" --n "${N}" --t-max "${TMAX}" --dt "${DT}" \
  --output-dir "${OUT}" \
  --vtk-every "${VTK_EVERY:-0}" \
  "$@"

echo "CSV: ${OUT}/nu_history.csv"
