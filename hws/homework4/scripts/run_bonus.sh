#!/usr/bin/env bash
# Extra credit drivers: (1) mesh×Ra sweep via bonus_convergence_sweep.py
# (3) Crank–Nicolson + SNES via convection_cn.py
# Run inside Firedrake (+ firedrake-ts for sweep). Set HW4_APPTAINER_* on clusters.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=inc_firedrake_apptainer.sh
source "${SCRIPT_DIR}/inc_firedrake_apptainer.sh"
HW4_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${HW4_ROOT}/python"

MODE="${1:-help}"
case "${MODE}" in
  sweep)
    shift || true
    echo "Running bonus mesh×Ra sweep -> output/bonus_convergence/ (override with --output-root)"
    hw4_python3 bonus_convergence_sweep.py "$@"
    ;;
  cn)
    shift || true
    echo "Running Crank–Nicolson / SNES driver convection_cn.py"
    hw4_python3 convection_cn.py "$@"
    ;;
  *)
    echo "Usage:"
    echo "  $0 sweep [--t-max 5000] [--dt 0.1] [--dry-run]   # extra credit (1)"
    echo "  $0 cn [--ra 1e4] [--n 32] [--t-max 500] ...       # extra credit (3)"
    echo ""
    echo "Set HW4_APPTAINER_SIF (+ bind) on clusters; see hw4/RUN.md §3b."
    exit 1
    ;;
esac
