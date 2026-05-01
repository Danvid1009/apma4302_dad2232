#!/usr/bin/env bash
# HW4 Q2: C PETSc biharm (three options files) + Firedrake biharm.py (three presets).
set -euo pipefail

HW4_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="${HW4_LOG_DIR:-${HW4_ROOT}/output/logs}"
mkdir -p "${LOG_DIR}"

echo "== Q2a: C / PETSc =="
if [[ -z "${PETSC_DIR:-}" ]]; then
  echo "ERROR: Set PETSC_DIR (and usually PETSC_ARCH) before running this script."
  exit 1
fi

cd "${HW4_ROOT}/c"
make biharm

run_c() {
  local name="$1"
  local opts="$2"
  echo "--- ./biharm ${opts} ---"
  ./biharm -options_file "${opts}" -log_view 2>&1 | tee "${LOG_DIR}/q2_c_${name}.log"
  grep -E 'SNESSolve|KSP|SNES' "${LOG_DIR}/q2_c_${name}.log" | tail -n 20 || true
}

run_c direct options_file_direct
run_c split_direct options_file_split_direct
run_c split_mg options_file_split_mg

echo "== Q2b: Firedrake biharm.py (requires active Firedrake env; MPI optional) =="
cd "${HW4_ROOT}/python"

for preset in direct split_direct split_mg; do
  export HW4_RESULT_DIR="${HW4_ROOT}/output/q2_firedrake/${preset}"
  mkdir -p "${HW4_RESULT_DIR}"
  echo "--- python biharm.py --preset ${preset} -> ${HW4_RESULT_DIR} ---"
  if command -v mpiexec >/dev/null 2>&1; then
    mpiexec -n 1 python3 biharm.py --preset "${preset}" 2>&1 | tee "${LOG_DIR}/q2_py_${preset}.log"
  else
    python3 biharm.py --preset "${preset}" 2>&1 | tee "${LOG_DIR}/q2_py_${preset}.log"
  fi
done

echo "Logs under ${LOG_DIR}"
echo "VTK under ${HW4_ROOT}/output/q2_firedrake/<preset>/biharm.pvd"
