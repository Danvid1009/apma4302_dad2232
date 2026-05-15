# Parallel Monte Carlo Option Pricing (PETSc)

**Final submission:** see [`SUBMISSION.md`](SUBMISSION.md) and [`submission/4302_Final_Project.pdf`](submission/4302_Final_Project.pdf). Git branch: `final-project`.

This project implements parallel Monte Carlo pricing for:
- European call
- Asian call
- Barrier call (up-and-out)

It is aligned to `PROJECT_OUTLINE.md` and includes:
- core simulation/pricing modules
- PETSc-based parallel pricing runner
- validation script against Black-Scholes
- ParaView-friendly path export starter

## Quick Start

1. Create environment and install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Optional quick check:

```bash
python -c "import petsc4py; print('petsc4py ok')"
```

2. Run pricing with PETSc/MPI:

```bash
mpirun -n 4 python scripts/run_pricing.py --option european --paths 200000
```

3. Validate European price against Black-Scholes:

```bash
mpirun -n 4 python scripts/validate_european.py --paths 400000
```

4. Export sample paths for ParaView:

```bash
python scripts/export_paraview_paths.py --paths 200 --steps 250 --out output/path_bundle.csv
```

5. Run convergence sweep (error vs \(N\)) and plot fitted slope:

```bash
mpirun -n 4 python scripts/run_convergence.py --paths-list 20000 50000 100000 200000 400000 --out output/convergence.csv
python scripts/plot_convergence.py --infile output/convergence.csv --out output/convergence.png
```

Tip: use `--reps 5` (or higher) on `run_convergence.py` for a more stable fitted slope.

6. Barrier level sweep (variance vs geometry) and error surface \((N,\sigma)\) vs Black–Scholes:

```bash
mpirun -n 4 python scripts/run_barrier_sweep.py --out output/barrier_sweep.csv
python scripts/plot_barrier_sweep.py --infile output/barrier_sweep.csv --out-prefix output/barrier_sweep

mpirun -n 4 python scripts/run_error_surface.py --out output/error_surface_european.csv
python scripts/plot_error_surface.py --infile output/error_surface_european.csv --out output/error_surface_abs_err.png
```

7. **Phase 2 (historical vol, optional Yahoo CSV):** see [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md). Quick demo without download:

```bash
python scripts/phase2_calibrate_and_price.py --csv data/raw/example_synthetic_ohlcv.csv --paths 100000 --out-json output/phase2_calibration_summary.json
```

With real data: `pip install -r requirements-phase2.txt`, run `scripts/fetch_yahoo_ohlcv.py`, then point `--csv` at `data/raw/<ticker>.csv`.

8. **Regenerate “gallery” HPC figures** (dark-theme dashboards into `visuals/`):

```bash
bash scripts/render_hpc_figures.sh
```

Produces `visuals/hpc_scaling_dashboard.png`, `hpc_convergence_triptych.png`, `hpc_speedup_vs_ideal.png` (if `output/strong_scaling.csv` exists), `hpc_weak_scaling_walltime.png` (if weak CSV exists), `phase2_spy_price_and_vol.png` (if SPY CSV exists), `hpc_barrier_geometry_story.png` (if `output/barrier_sweep.csv` exists), and regenerates `output/paraview/vtp/mpi_domain_paths.vtp` for ParaView.

## Project Layout

- `src/montecarlo/`: core models and pricing engine
- `scripts/`: command-line entry points for experiments
- `output/`: generated experiment outputs
- `visuals/`: generated PNG figures and ParaView captures
- `data/`: optional input/output datasets

## Reporting and tracking

- **Full project summary + PETSc checklist:** [`docs/PROJECT_SUMMARY_AND_PETSC_CHECKLIST.md`](docs/PROJECT_SUMMARY_AND_PETSC_CHECKLIST.md) — what is in the repo and what to run when PETSc is available.
- **Experiment registry (start here):** `docs/HPC_EXPERIMENT_REGISTRY.md` — catalog of scripts, class HPC topics (MPI, collectives, scaling), outputs, and maintenance checklist.
- **Append-only run log:** `output/registry/experiment_runs.csv` — log each campaign with `python scripts/log_experiment_run.py --key ... --output ...`
- **Long / background jobs:** `docs/LONG_RUNNING_JOBS.md` — probe Cursor terminal snapshots with `python scripts/terminal_probe.py --auto`
- **Optional cluster / Apptainer (few nodes, Firedrake+TS image):** [`docs/RUN_LATER_CLUSTER_APPTAINER.md`](docs/RUN_LATER_CLUSTER_APPTAINER.md) — Slurm, `sbatch`, binds; **not** required for the Monte Carlo project.
- Phase 1 results summary: `RESULTS_PHASE1.md`
- ParaView / VTK: [`PROJECT_OUTLINE.md` Appendix A](PROJECT_OUTLINE.md#appendix-a-paraview-and-vtk-workflows) (redirect stubs: [`PARAVIEW_GUIDE_PETSC.md`](PARAVIEW_GUIDE_PETSC.md), [`PARAVIEW_GUIDE.md`](PARAVIEW_GUIDE.md))
- **Optional market data (Phase 2):** [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) — Yahoo OHLCV via `yfinance`; no data needed for core MC/HPC work.

## Next Build Steps

- Run scaling matrix, barrier sweep, and error surface **locally** (`mpirun` + scripts); append rows to `output/registry/experiment_runs.csv`
- Optional later: small Slurm / Apptainer job — [`docs/RUN_LATER_CLUSTER_APPTAINER.md`](docs/RUN_LATER_CLUSTER_APPTAINER.md)
- Optional: control variates / QMC beyond antithetic
