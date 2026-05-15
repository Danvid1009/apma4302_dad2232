#!/usr/bin/env bash
# Regenerate all "cool" HPC + Phase2 figures into visuals/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
mkdir -p visuals

# Prefer Apple /usr/bin/python3 on macOS when Homebrew's python3 is first on PATH but lacks your plot stack.
# Override: export MCTPLOT_PYTHON="$(which python3.9)"
PY="${MCTPLOT_PYTHON:-/usr/bin/python3}"

if [[ -f output/strong_scaling.csv ]]; then
  "$PY" scripts/plot_hpc_speedup_vs_ideal.py \
    --csv output/strong_scaling.csv \
    --out visuals/hpc_speedup_vs_ideal.png
fi

if [[ -f output/weak_scaling.csv ]]; then
  "$PY" scripts/plot_hpc_weak_scaling_walltime.py \
    --csv output/weak_scaling.csv \
    --out visuals/hpc_weak_scaling_walltime.png
fi

"$PY" scripts/export_mpi_domain_paths_vtp.py \
  --mpi-ranks 8 --paths-per-rank 14 --steps 100 \
  --out output/paraview/vtp/mpi_domain_paths.vtp

if [[ -f output/strong_scaling.csv ]] && [[ -f output/weak_scaling.csv ]]; then
  "$PY" scripts/plot_hpc_scaling_dashboard.py \
    --strong-csv output/strong_scaling.csv \
    --weak-csv output/weak_scaling.csv \
    --out visuals/hpc_scaling_dashboard.png
fi

if [[ -f output/convergence.csv ]]; then
  "$PY" scripts/plot_hpc_convergence_triptych.py \
    --infile output/convergence.csv \
    --out visuals/hpc_convergence_triptych.png
fi

if [[ -f data/raw/SPY.csv ]]; then
  "$PY" scripts/plot_spy_volatility_regime.py \
    --csv data/raw/SPY.csv \
    --out visuals/phase2_spy_price_and_vol.png
fi

if [[ -f output/barrier_sweep.csv ]]; then
  "$PY" scripts/plot_hpc_barrier_story.py \
    --infile output/barrier_sweep.csv \
    --out visuals/hpc_barrier_geometry_story.png
fi

if [[ -f output/poisson_pc_sweep.csv ]]; then
  "$PY" scripts/plot_poisson_pc_sweep.py \
    --infile output/poisson_pc_sweep.csv \
    --out-prefix visuals/poisson_pc_sweep
fi

echo "Done. See visuals/ and output/paraview/vtp/mpi_domain_paths.vtp"
