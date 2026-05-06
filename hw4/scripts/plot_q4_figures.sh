#!/usr/bin/env bash
# Regenerate HW4 Q4 Nusselt figures (demo curves, or from nu_history.csv after run_q4.sh).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIG="${ROOT}/doc/figures"
PY="${ROOT}/.venv-plot/bin/python"
if [[ ! -x "${PY}" ]]; then
  echo "Create plotting venv first:"
  echo "  python3 -m venv ${ROOT}/.venv-plot && ${ROOT}/.venv-plot/bin/pip install numpy matplotlib"
  exit 1
fi
cd "${ROOT}/python"
if [[ "${1:-}" == "--demo" || $# -eq 0 ]]; then
  exec "${PY}" plot_q4_nusselt.py --figures-dir "${FIG}" --demo
fi
exec "${PY}" plot_q4_nusselt.py --figures-dir "${FIG}" "$@"
