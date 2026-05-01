#!/usr/bin/env bash
# HW4 Q4: convection DAE (firedrake-ts). Override RA, N, TMAX, DT via env vars.
set -euo pipefail

HW4_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${HW4_ROOT}/python"

RA="${RA:-1e2}"
N="${N:-64}"
TMAX="${TMAX:-1e5}"
DT="${DT:-0.1}"
OUT="${HW4_OUT:-${HW4_ROOT}/output/q4_ra${RA}_n${N}}"

echo "Running convection.py Ra=${RA} N=${N} t_max=${TMAX} dt=${DT} -> ${OUT}"
if command -v mpiexec >/dev/null 2>&1; then
  mpiexec -n 1 python3 convection.py \
    --ra "${RA}" --n "${N}" --t-max "${TMAX}" --dt "${DT}" \
    --output-dir "${OUT}" \
    --vtk-every "${VTK_EVERY:-0}" \
    "$@"
else
  python3 convection.py \
    --ra "${RA}" --n "${N}" --t-max "${TMAX}" --dt "${DT}" \
    --output-dir "${OUT}" \
    --vtk-every "${VTK_EVERY:-0}" \
    "$@"
fi

echo "CSV: ${OUT}/nu_history.csv"
