#!/usr/bin/env bash
# Regenerate VTK time series tuned for "hero" ParaView stills (warp, glyphs, colormaps).
# Usage (from repo project root — directory containing scripts/ and src/):
#   bash scripts/paraview_pretty_runs.sh
#
# Override defaults:
#   PY=/usr/bin/python3 MPIEXEC=mpiexec NP=4 bash scripts/paraview_pretty_runs.sh

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

: "${PETSC_DIR:=/Users/dan/Desktop/Columbia/HPC_4302/petsc}"
: "${PETSC_ARCH:=apma4302-pkgs-opt}"
: "${PY:=/usr/bin/python3}"
: "${MPIEXEC:=mpiexec}"
: "${NP:=4}"

export PYTHONPATH="${PETSC_DIR}/${PETSC_ARCH}/lib:${ROOT}/src:${PYTHONPATH:-}"
export DYLD_LIBRARY_PATH="${PETSC_DIR}/${PETSC_ARCH}/lib:${DYLD_LIBRARY_PATH:-}"

echo "== Poisson snapshot (high-res surface / contour in ParaView) =="
"$PY" -u scripts/smoke_paraview_export.py --nx 151 --ny 151 --outdir output/paraview/vtp --basename pretty_poisson

echo "== Heat — smooth IC, mid-late frames =="
"$MPIEXEC" -n "$NP" "$PY" -u scripts/run_heat_paraview.py \
  --nx 161 --ny 161 --steps 300 --kappa 0.18 --save-every 30 \
  --outdir output/paraview/vtp --series-name pretty_heat

echo "== Advection–diffusion — twin blob + diagonal drift + velocity in VTP =="
"$MPIEXEC" -n "$NP" "$PY" -u scripts/run_advection_diffusion_paraview.py \
  --nx 161 --ny 161 --steps 360 --dt 6.5e-4 --nu 0.008 \
  --vel-x 1.05 --vel-y 0.72 --save-every 36 --ic twin \
  --outdir output/paraview/vtp --series-name pretty_advection

echo "== Optional: MPI path metaphor (ParaView Tube filter) =="
"$PY" -u scripts/export_mpi_domain_paths_vtp.py \
  --mpi-ranks 8 --paths-per-rank 20 --steps 140 \
  --out output/paraview/vtp/pretty_mpi_domain_paths.vtp

echo "Done. Open in ParaView:"
echo "  output/paraview/vtp/pretty_poisson.pvd"
echo "  output/paraview/vtp/pretty_heat.pvd"
echo "  output/paraview/vtp/pretty_advection.pvd"
echo "  output/paraview/vtp/pretty_mpi_domain_paths.vtp"
echo "See PROJECT_OUTLINE.md Appendix A (Pretty hero kit subsection)."
