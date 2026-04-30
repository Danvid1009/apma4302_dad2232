# Parallel Monte Carlo Option Pricing (PETSc)

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

## Project Layout

- `src/montecarlo/`: core models and pricing engine
- `scripts/`: command-line entry points for experiments
- `output/`: generated experiment outputs
- `visuals/`: generated PNG figures and ParaView captures
- `data/`: optional input/output datasets

## Reporting Helpers

- Phase 1 results summary: `RESULTS_PHASE1.md`
- ParaView workflow steps: `PARAVIEW_GUIDE.md`

## Next Build Steps

- Add variance reduction toggles (antithetic/control variates)
- Add structured-grid exporter for error surfaces
