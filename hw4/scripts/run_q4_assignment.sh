#!/usr/bin/env bash
# HW4 Q4 — run parts (a), (b), and (c) with the handout parameters.
#
# Uses: BDF-2, monolithic direct + MUMPS (via convection.py defaults), Q1 mesh
# (UnitSquareMesh N×N quadrilateral=True, CG 1).
#
# Prerequisites (same as run_q4.sh):
#   - Firedrake + firedrake-ts in PATH, OR
#   - export HW4_APPTAINER_SIF="$HOME/firedrake-ts.sif"
#     export HW4_APPTAINER_BIND="$HOME:$HOME"
#
# Optional overrides (all times are dimensionless t as in the prompt):
#   DT=0.01                 Default initial TS step (--ts-dt hint; adaptive adjusts). Use TS_FIXED_STEP=1
#                           so DT is the actual fixed step (~t_max/dt steps).
#   TMAX_4A=100000          part (a) end time (default 100000 = 10^5)
#   TMAX_4B=100000          part (b) end time per Ra (raise if Nu not steady)
#   TMAX_4C=100000          part (c) end time per mesh (raise toward Blankenbach steady Nu)
#   VTK_EVERY=0             set >0 for ParaView dumps (large output)
#   TS_RTOL / TS_ATOL       forwarded to convection.py (--ts-rtol / --ts-atol). Part (b)
#                           at high Ra is stiff: adaptive TS often uses tiny dt → many
#                           iterations. Loosening (e.g. TS_RTOL=1e-4 TS_ATOL=1e-8) speeds
#                           runs at some accuracy cost; lowering TMAX_4B shortens drafts.
#   TS_MAX_STEPS            e.g. 50000 → PETSc stops after that many TS steps (t may be
#                           below t_max); use for wall-clock caps / smoke tests.
#   TS_FIXED_STEP=1        fixed timestep (PETSc ts_adapt_type=none); pair with larger DT,
#                           e.g. DT=500 → ~200 steps to t_max=1e5 (may fail at high Ra if dt too large).
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=inc_firedrake_apptainer.sh
source "${SCRIPT_DIR}/inc_firedrake_apptainer.sh"

HW4_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DT="${DT:-0.01}"
VTK_EVERY="${VTK_EVERY:-${VTK:-0}}"

TMAX_4A="${TMAX_4A:-100000}"
TMAX_4B="${TMAX_4B:-100000}"
TMAX_4C="${TMAX_4C:-100000}"

run_one() {
  local ra="$1" n="$2" tmax="$3" tag="$4"
  echo ""
  echo ">>> ${tag}  Ra=${ra}  N=${n}  t_max=${tmax}  dt=${DT}"
  RA="${ra}" N="${n}" TMAX="${tmax}" DT="${DT}" VTK_EVERY="${VTK_EVERY}" \
    HW4_OUT="${HW4_ROOT}/output/${tag}" \
    "${SCRIPT_DIR}/run_q4.sh"
}

echo "HW4 Q4 assignment driver — outputs under ${HW4_ROOT}/output/"
if [[ -n "${HW4_APPTAINER_SIF:-${FIREDRAKE_TS_SIF:-}}" ]]; then
  echo "Apptainer SIF: ${HW4_APPTAINER_SIF:-${FIREDRAKE_TS_SIF}}"
fi

echo ""
echo "=== (a) Ra = 10^2, 64×64, t_max = ${TMAX_4A} (show Nu → 1, non-convective steady state) ==="
run_one 1e2 64 "${TMAX_4A}" "q4_4a_Ra1e2_N64"

echo ""
echo "=== (b) Ra = 10^4, 10^5, 10^6, 64×64, t_max = ${TMAX_4B} each ==="
for ra in 1e4 1e5 1e6; do
  run_one "${ra}" 64 "${TMAX_4B}" "q4_4b_Ra${ra}_N64"
done

echo ""
echo "=== (c) Ra = 10^4, N = 16, 32, 64, 128, t_max = ${TMAX_4C} each (compare late Nu to Blankenbach 4.884) ==="
for n in 16 32 64 128; do
  run_one 1e4 "${n}" "${TMAX_4C}" "q4_4c_Ra1e4_N${n}"
done

echo ""
echo "All legs finished. nu_history.csv locations:"
echo "  (a) ${HW4_ROOT}/output/q4_4a_Ra1e2_N64/nu_history.csv"
echo "  (b) ${HW4_ROOT}/output/q4_4b_Ra1e4_N64/nu_history.csv"
echo "      ${HW4_ROOT}/output/q4_4b_Ra1e5_N64/nu_history.csv"
echo "      ${HW4_ROOT}/output/q4_4b_Ra1e6_N64/nu_history.csv"
echo "  (c) ${HW4_ROOT}/output/q4_4c_Ra1e4_N16|N32|N64|N128/nu_history.csv"
echo ""
echo "Plot (after rsync to laptop): hw4/python/plot_q4_nusselt.py --help"
echo "If Nu has not settled, re-run one leg with a larger TMAX_* or smaller DT."
