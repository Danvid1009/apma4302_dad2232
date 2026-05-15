# Phase 3 Runbook

This runbook generates report-ready convergence and scaling figures.

## 1) Environment

```bash
export PETSC_DIR="/Users/dan/Desktop/Columbia/HPC_4302/petsc"
export PETSC_ARCH="apma4302-pkgs-opt"
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib"
export PATH="$PETSC_DIR/$PETSC_ARCH/bin:$PETSC_DIR/bin:$PATH"
```

## 2) Convergence sweep

```bash
python scripts/run_poisson_convergence.py --sizes 21 41 81 121 161 --out output/poisson_convergence.csv
python scripts/plot_phase3.py --conv output/poisson_convergence.csv --outdir visuals
```

Outputs:
- `output/poisson_convergence.csv`
- `visuals/spatial_convergence_loglog.png`

## 3) Strong scaling sweep

```bash
# Use a Poisson-specific filename so you do not overwrite Monte Carlo scaling CSVs.
OUT=output/poisson_strong_scaling.csv
mpiexec -n 1 python scripts/run_poisson_timed.py --nx 401 --ny 401 --out "$OUT"
mpiexec -n 2 python scripts/run_poisson_timed.py --nx 401 --ny 401 --out "$OUT"
mpiexec -n 4 python scripts/run_poisson_timed.py --nx 401 --ny 401 --out "$OUT"
mpiexec -n 8 python scripts/run_poisson_timed.py --nx 401 --ny 401 --out "$OUT"
python scripts/plot_phase3.py --scaling "$OUT" --outdir visuals
```

Outputs:
- `output/poisson_strong_scaling.csv` (example path above)
- `visuals/strong_scaling_speedup.png`
- `visuals/strong_scaling_efficiency.png`
