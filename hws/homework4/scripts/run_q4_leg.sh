#!/usr/bin/env bash
# HW4 Q4 — run a single assignment leg (save outputs as you go on the cluster).
#
# Same env as run_q4_assignment.sh: DT (default 0.01), TMAX_4A/B/C, TS_*, VTK_EVERY, Apptainer, etc.
# Optional: DS_TOP / DS_BOT forwarded to convection.py (facet IDs for y=1 / y=0).
#
# Usage:
#   ./run_q4_leg.sh list              # printed menu + rough progress guide
#   ./run_q4_leg.sh 1                 # leg 1 only
#   Q4_LEG_LOG=hw4/output/logs/q4_leg03.log ./run_q4_leg.sh 3
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HW4_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DT="${DT:-0.01}"
VTK_EVERY="${VTK_EVERY:-${VTK:-0}}"
TMAX_4A="${TMAX_4A:-100000}"
TMAX_4B="${TMAX_4B:-100000}"
TMAX_4C="${TMAX_4C:-100000}"
DS_TOP="${DS_TOP:-}"
DS_BOT="${DS_BOT:-}"

# Reference PETSc TS *accepted-step* counts from a sample run log (Ra, N, t_max=1e5, adaptive TS, dt≈1e−2 hint).
# Your machine/cluster will differ in wall time per step; use counts only to sanity-check progress via printed "N TS dt" lines.
declare -ra LEG_LABEL=(
  "(a) Ra=1e2  N=64   t_max=${TMAX_4A}"
  "(b) Ra=1e4  N=64   t_max=${TMAX_4B}"
  "(b) Ra=1e5  N=64   t_max=${TMAX_4B}"
  "(b) Ra=1e6  N=64   t_max=${TMAX_4B}"
  "(c) Ra=1e4  N=16   t_max=${TMAX_4C}"
  "(c) Ra=1e4  N=32   t_max=${TMAX_4C}"
  "(c) Ra=1e4  N=64   t_max=${TMAX_4C}"
  "(c) Ra=1e4  N=128  t_max=${TMAX_4C}"
)
declare -ra LEG_TAGS=(
  "q4_4a_Ra1e2_N64"
  "q4_4b_Ra1e4_N64"
  "q4_4b_Ra1e5_N64"
  "q4_4b_Ra1e6_N64"
  "q4_4c_Ra1e4_N16"
  "q4_4c_Ra1e4_N32"
  "q4_4c_Ra1e4_N64"
  "q4_4c_Ra1e4_N128"
)
# "~steps" column: friend's log final "N TS dt …" index for completed legs 1–3,5–8; Ra=1e6 often >> others.
declare -ra LEG_STEPS_HINT=(
  "~1.1×10³ TS steps"
  "~1.3×10³"
  "~2.9×10³"
  "≫others; stiff O(10⁴–10⁵) TS steps"
  "~1.2×10³"
  "~1.2×10³"
  "~1.3×10³"
  "~1.3×10³"
)

print_list() {
  echo "HW4 Q4 legs (run one at a time) — outputs: ${HW4_ROOT}/output/<tag>/nu_history.csv"
  echo ""
  printf " %-2s  %-42s  %-22s  %s\n" "ID" "Case" "~TS steps (sample)" "Output tag"
  echo "------------------------------------------------------------------------------------------------"
  local i
  for i in "${!LEG_LABEL[@]}"; do
    local id=$((i + 1))
    printf " %-2s  %-42s  %-22s  %s\n" "${id}" "${LEG_LABEL[i]}" "${LEG_STEPS_HINT[i]}" "${LEG_TAGS[i]}"
  done
  echo ""
  echo "Relative CPU: each PETSc step does a sparse LU (MUMPS); wall time scales strongly with mesh (N³–N⁴)."
  echo "  Order-of-magnitude expectation: legs 5–8 with N=16 cheapest per step; leg 8 (N=128) heaviest among (c)."
  echo "Leg 4 (Ra=1e6) is usually the longest overall job at t_max=1e5—leave extra wall time or tighten TS_* if needed."
  echo ""
  echo "Current env (override before calling): DT=${DT}  TMAX_4A=${TMAX_4A}  TMAX_4B=${TMAX_4B}  TMAX_4C=${TMAX_4C}"
  [[ -n "${DS_BOT}" || -n "${DS_TOP}" ]] && echo "Facet IDs: DS_BOT=${DS_BOT:-<unset>}  DS_TOP=${DS_TOP:-<unset>}"
  echo 'Optional log: Q4_LEG_LOG=hw4/output/logs/q4_leg01.log '"${SCRIPT_DIR}/run_q4_leg.sh 1"
}

run_leg() {
  local id="$1"
  if [[ ! "$id" =~ ^[1-8]$ ]]; then
    echo "Invalid leg '${id}'. Use: $0 list  or  ID in 1..8" >&2
    exit 2
  fi
  local i=$((id - 1))
  local ra n tmax tag
  case "${id}" in
    1) ra=1e2 n=64 tmax="${TMAX_4A}" tag="${LEG_TAGS[i]}" ;;
    2) ra=1e4 n=64 tmax="${TMAX_4B}" tag="${LEG_TAGS[i]}" ;;
    3) ra=1e5 n=64 tmax="${TMAX_4B}" tag="${LEG_TAGS[i]}" ;;
    4) ra=1e6 n=64 tmax="${TMAX_4B}" tag="${LEG_TAGS[i]}" ;;
    5) ra=1e4 n=16 tmax="${TMAX_4C}" tag="${LEG_TAGS[i]}" ;;
    6) ra=1e4 n=32 tmax="${TMAX_4C}" tag="${LEG_TAGS[i]}" ;;
    7) ra=1e4 n=64 tmax="${TMAX_4C}" tag="${LEG_TAGS[i]}" ;;
    8) ra=1e4 n=128 tmax="${TMAX_4C}" tag="${LEG_TAGS[i]}" ;;
  esac

  echo ""
  echo "======== Leg ${id}/${#LEG_TAGS[@]}: ${LEG_LABEL[i]} ========"
  echo "Tag: ${tag}   ~steps (sample): ${LEG_STEPS_HINT[i]}"
  echo "CSV: ${HW4_ROOT}/output/${tag}/nu_history.csv"
  [[ -n "${Q4_LEG_LOG:-}" ]] && echo "Appending stdout/stderr to: ${Q4_LEG_LOG}"

  mkdir -p "${HW4_ROOT}/output/logs"
  export RA="${ra}" N="${n}" TMAX="${tmax}" DT="${DT}" VTK_EVERY="${VTK_EVERY}"

  [[ -n "${HW4_APPTAINER_SIF:-${FIREDRAKE_TS_SIF:-}}" ]] && echo "Using Apptainer SIF: ${HW4_APPTAINER_SIF:-${FIREDRAKE_TS_SIF}}"

  if [[ -n "${Q4_LEG_LOG:-}" ]]; then
    mkdir -p "$(dirname "${Q4_LEG_LOG}")"
    local st
    HW4_OUT="${HW4_ROOT}/output/${tag}" \
      RA="${ra}" N="${n}" TMAX="${tmax}" DT="${DT}" VTK_EVERY="${VTK_EVERY}" \
      "${SCRIPT_DIR}/run_q4.sh" 2>&1 | tee -a "${Q4_LEG_LOG}"
    st="${PIPESTATUS[0]}"
    return "${st}"
  else
    HW4_OUT="${HW4_ROOT}/output/${tag}" "${SCRIPT_DIR}/run_q4.sh"
  fi
}

case "${1:-}" in
  "" | -h | --help | help | list | ls)
    print_list
    if [[ "${1:-}" =~ ^(-h|--help|help)$ ]]; then
      echo ""
      echo "Run one leg:  $0 <1-8>"
    fi
    ;;
  *)
    run_leg "$1"
    ;;
esac
