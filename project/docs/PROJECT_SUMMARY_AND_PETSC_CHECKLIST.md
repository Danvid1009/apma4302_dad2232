# Project summary and PETSc access checklist

Single reference for **what is in this repository** and **what to run once your PETSc Python and MPI are available**.

---

## 1. What this project is

**Parallel Monte Carlo option pricing** (European, Asian, barrier) using **PETSc / MPI** for path batching and **global reductions** (sums for mean and variance, control-variate moments). **Phase 2** adds historical volatility from CSV, roll-forward repricing, optional quote comparison, and HPC-style reporting (registry, dashboards).

**Math spine:** GBM paths → payoffs → discount → `MPI.SUM`-style reductions → price + standard error + CI.

**Report / planning:** [`PROJECT_OUTLINE.md`](../PROJECT_OUTLINE.md) (full plan + syllabus §0.15 + **Appendix A** ParaView/VTK). **Experiment catalog:** [`HPC_EXPERIMENT_REGISTRY.md`](HPC_EXPERIMENT_REGISTRY.md). **Phase 1 → Phase 2 narrative:** [`PHASE1_TO_PHASE2.md`](PHASE1_TO_PHASE2.md).

---

## 2. Repository layout (high level)

| Area | Role |
| --- | --- |
| [`src/montecarlo/`](../src/montecarlo/) | Core library: `gbm.py`, `payoffs.py`, `engine.py` (MC + Asian CV), `black_scholes.py`, `historical_data.py`, `vol_calibration.py`. |
| [`scripts/`](../scripts/) | All drivers: pricing, validation, scaling, convergence, Phase 2, plots, Poisson/course-bridge tools, registry logging. |
| [`docs/`](../docs/) | Runbooks, data sources, experiment registry, cluster deferral, this file. |
| [`data/raw/`](../data/raw/) | OHLCV CSVs (e.g. `SPY.csv`), `quotes_template.csv`, synthetic example CSV. |
| [`output/`](../output/) | CSV/JSON results, `registry/experiment_runs.csv`, ParaView **VTP/PVD** under `output/paraview/vtp/` (single directory). |
| [`visuals/`](../visuals/) | Dashboard-style PNGs from `plot_hpc_*` and related scripts. |
| [`overleaf/`](../overleaf/) | LaTeX draft assets. |
| [`RESULTS_PHASE1.md`](../RESULTS_PHASE1.md) | Human-readable baseline numbers / narrative snapshots. |

---

## 3. Dependencies

| File | Contents |
| --- | --- |
| [`requirements.txt`](../requirements.txt) | `numpy`, `scipy`, **`petsc4py`**, **`mpi4py`**, `matplotlib` — **core stack**. |
| [`requirements-phase2.txt`](../requirements-phase2.txt) | `yfinance`, `pandas` — **only** for Yahoo fetch / optional CSV tooling. |

Install core (in the **same** environment where PETSc was built or matched):

```bash
pip install -r requirements.txt
```

---

## 4. Your PETSc environment (this course layout)

[`scripts/check_petsc_env.py`](../scripts/check_petsc_env.py) expects:

- **`PETSC_DIR`:** `/Users/dan/Desktop/Columbia/HPC_4302/petsc`
- **`PETSC_ARCH`:** `apma4302-pkgs-opt`
- **`PYTHONPATH`:** includes `.../apma4302-pkgs-opt/lib` (so `petsc4py` resolves)
- **`PATH`:** includes `.../apma4302-pkgs-opt/bin` (for `mpiexec` / PETSc tools aligned with that build)

**First command when you have access:**

```bash
export PETSC_DIR=/Users/dan/Desktop/Columbia/HPC_4302/petsc
export PETSC_ARCH=apma4302-pkgs-opt
# Then set PYTHONPATH / PATH per your local install (see check script messages).

python scripts/check_petsc_env.py
```

Use the **Python interpreter that matches this PETSc build** for every `python` / `mpirun python` below (often the one inside the PETSc arch tree or the conda/env you use for APMA4302).

**Imports for both `petsc4py` and `montecarlo`:** from the project root, combine the arch `lib` with `src`, for example:

```bash
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
```

### 4.1 First import, dyld, and plotting (macOS-friendly)

- **Cold start:** the first `from petsc4py import PETSc` in a fresh process can sit for a long time with **no stdout** while dyld loads `libpetsc`, MPI, and friends; subsequent work in the **same** interpreter is usually fast.
- **`DYLD_LIBRARY_PATH` (optional):** if loads stall or libraries resolve oddly, export `DYLD_LIBRARY_PATH="$PETSC_DIR/$PETSC_ARCH/lib:${DYLD_LIBRARY_PATH:-}"` before the first import of the day.
- **Live prints:** use **`python -u`** or **`PYTHONUNBUFFERED=1`** if you want output while shared objects load.
- **Figure-only scripts:** drivers that only use Matplotlib do not need MPI, but they still need a Python that has NumPy/Matplotlib. If Homebrew’s `python3` is first on `PATH` without those wheels, use Apple **`/usr/bin/python3`** or set **`MCTPLOT_PYTHON`** (see [`render_hpc_figures.sh`](../scripts/render_hpc_figures.sh)).

### 4.2 PETSc-heavy course bridge (compact run map)

These are the **ready PETSc-native knobs** for the report (Vec, SNES, Mat/Vec/KSP/PC, TS + ParaView). They are mostly **not** a second “PETSc C rewrite” of the Monte Carlo pricer: the main pricer already runs under `mpiexec` with `petsc4py` reductions. A deeper PETSc rewrite (paths as `Vec`, SDEs in `TS`, etc.) would be **new code**, not another script you already have.

*Below, `python` means the interpreter that matches your `petsc4py` ABI (§4). Project root = directory containing `scripts/` and `src/`.*

#### 1. Parallel `PETSc.Vec` + global sum (LA / Vec week)

**What it adds:** Makes the “global combine” story explicit with a partitioned `Vec` and a reduction, instead of only ad hoc scalar allreduces.

```bash
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
mpiexec -n 4 python -u scripts/demo_petsc_vec_global_sum.py
```

#### 2. SNES implied volatility (Newton / nonlinear week)

**What it adds:** Scalar nonlinear solve \(\mathrm{BS}(\sigma) = V_{\mathrm{mkt}}\) with PETSc **SNES**; optional hand Jacobian (**vega**) vs finite differences.

Pick a synthetic market price equal to Black–Scholes at e.g. \(\sigma = 0.2\) with defaults \(S_0 = K = 100\), \(r = 0.05\), \(T = 1\) (about **10.4506**; compute once in a REPL if you want the full double).

```bash
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
mpiexec -n 1 python -u scripts/implied_vol_snes.py \
  --market-price 10.4506 \
  --analytic-jacobian
```

Optional JSON:

```bash
mpiexec -n 1 python -u scripts/implied_vol_snes.py \
  --market-price 10.4506 --analytic-jacobian \
  --out output/implied_vol_snes.json
```

#### 3. Poisson: `Mat` + `Vec` + `KSP` + `PC` (Krylov / PC / PDE weeks)

**What it adds:** Same PETSc objects as typical course labs; good for PC comparison or scaling narratives alongside Monte Carlo.

Timed solve (append CSV):

```bash
mpiexec -n 4 python -u scripts/run_poisson_timed.py \
  --nx 201 --ny 201 --ksp-type cg --pc-type jacobi \
  --out output/poisson_strong_append.csv --label poisson2d
```

PC sweep (several preconditioners):

```bash
mpiexec -n 4 python -u scripts/run_poisson_pc_sweep.py \
  --nx 65 --ny 65 --ksp-type cg --pc-types jacobi sor \
  --out output/poisson_pc_sweep.csv
```

KSP sweep (several Krylov solvers, fixed PC — good for “CG vs GMRES” discussion):

```bash
mpiexec -n 4 python -u scripts/run_poisson_ksp_sweep.py \
  --nx 65 --ny 65 --pc-type jacobi --ksp-types cg gmres bcgsl \
  --out output/poisson_ksp_sweep.csv
python scripts/plot_poisson_ksp_sweep.py --infile output/poisson_ksp_sweep.csv --out-prefix visuals/poisson_ksp_sweep
```

Convergence in mesh size:

```bash
mpiexec -n 1 python -u scripts/run_poisson_convergence.py \
  --sizes 21 41 81 121 --out output/poisson_convergence.csv
```

Details and plots: [`PHASE3_RUNBOOK.md`](PHASE3_RUNBOOK.md).

#### 4. Poisson solve vs Monte Carlo (“two cultures” / bottleneck story)

**What it adds:** One CSV row comparing **one KSP solve** vs **one European MC** on the same ranks — solver-bound vs sampling-bound behavior.

```bash
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
mpiexec -n 4 python -u scripts/run_sampling_vs_solver_benchmark.py \
  --nx 129 --ny 129 --ksp-type cg --pc-type jacobi \
  --paths 80000 --steps 252 \
  --out output/sampling_vs_solver_benchmark.csv
```

#### 5. Heat / advection–diffusion + ParaView (`TS` / advection weeks)

**What it adds:** `TS`-style time stepping in PETSc for PDEs, same toolchain as visualization lectures.

```bash
mpiexec -n 4 python -u scripts/run_heat_paraview.py --nx 121 --ny 121 --steps 250 --outdir output/paraview/vtp
mpiexec -n 4 python -u scripts/run_advection_diffusion_paraview.py --nx 121 --ny 121 --steps 320 --outdir output/paraview/vtp
```

Defaults use **`cg` + `jacobi`** (parallel-safe on typical course builds). For nonsymmetric stencils you may try **`--ksp-type gmres --pc-type jacobi`**; avoid parallel **`ilu`** unless your PETSc install provides a compatible factorization backend.

Optional **matplotlib** views of solver logs (after the runs exist):

```bash
python scripts/plot_poisson_pc_sweep.py --infile output/poisson_pc_sweep.csv --out-prefix visuals/poisson_pc_sweep
python scripts/plot_diffusion_solver_history.py --infile output/paraview/vtp/heat_solver_history.csv
python scripts/plot_diffusion_solver_history.py --infile output/paraview/vtp/advection_diffusion_solver_history.csv
```

See [`PROJECT_OUTLINE.md` Appendix A](../PROJECT_OUTLINE.md#appendix-a-paraview-and-vtk-workflows) for opening `.pvd` files (stub: [`PARAVIEW_GUIDE_PETSC.md`](../PARAVIEW_GUIDE_PETSC.md)).

#### 6. Extra “slide-ready” HPC figures (Matplotlib + one ParaView metaphor)

After **`output/strong_scaling.csv`** / **`output/weak_scaling.csv`** exist (`run_scaling.py`):

```bash
python scripts/plot_hpc_speedup_vs_ideal.py --csv output/strong_scaling.csv --out visuals/hpc_speedup_vs_ideal.png
python scripts/plot_hpc_weak_scaling_walltime.py --csv output/weak_scaling.csv --out visuals/hpc_weak_scaling_walltime.png
```

**MPI path batches in ParaView** (synthetic rank field; good next to scaling plots):

```bash
python scripts/export_mpi_domain_paths_vtp.py --mpi-ranks 8 --paths-per-rank 16 --steps 120 --out output/paraview/vtp/mpi_domain_paths.vtp
```

Tube + color by `mpi_rank`: see [`PROJECT_OUTLINE.md` Appendix A](../PROJECT_OUTLINE.md#appendix-a-paraview-and-vtk-workflows), subsection **5) “Cool picture” recipes** → **E) Path / tube metaphor**. Or run **`bash scripts/render_hpc_figures.sh`** (uses `MCTPLOT_PYTHON` for plot scripts; still writes the `.vtp`).

---

## 5. When PETSc is available: ordered checklist

### A. Environment and smoke

1. Run **`python scripts/check_petsc_env.py`** → should print `Status: READY` and `petsc4py_import=OK`.
2. From project root, set **`PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"`** (or install the package in editable mode) so **`petsc4py`** and **`montecarlo`** both import.
3. Quick MC smoke (replace `python` with your PETSc Python):

   ```bash
   export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
   mpiexec -n 2 python -u scripts/validate_european.py --paths 50000
   ```

### B. Phase 1 experiments (core report evidence)

Run and archive outputs under `output/` (and copy key numbers into [`RESULTS_PHASE1.md`](../RESULTS_PHASE1.md) if they change).

| Goal | Commands (typical) |
| --- | --- |
| Validation vs Black–Scholes | `mpiexec -n 4 python scripts/validate_european.py --paths 400000` |
| Convergence vs \(N\) | `mpiexec -n 4 python scripts/run_convergence.py --out output/convergence.csv` (+ optional `--antithetic`) |
| Plot convergence | `python scripts/plot_convergence.py --infile output/convergence.csv --out output/convergence.png` |
| Strong / weak scaling | `mpiexec -n 4 python scripts/run_scaling.py ...` then `python scripts/plot_scaling.py` |
| Variance reduction | `python scripts/run_variance_reduction_bench.py` (see script `--help`) |
| Barrier sweep | `mpiexec -n 4 python scripts/run_barrier_sweep.py --out output/barrier_sweep.csv` + `plot_barrier_sweep.py` |
| Error surface | `mpiexec -n 4 python scripts/run_error_surface.py` + `plot_error_surface.py` |
| Asian control variate | `mpiexec -n 4 python scripts/run_asian_cv_control_variate.py` (see `--help`) |
| HPC figure gallery | `bash scripts/render_hpc_figures.sh` |

### C. Phase 2 (data + realism)

| Goal | Notes |
| --- | --- |
| Fetch OHLCV | `pip install -r requirements-phase2.txt` then `python scripts/fetch_yahoo_ohlcv.py ...` → `data/raw/*.csv` |
| Calibrate + price | `PYTHONPATH=src python scripts/phase2_calibrate_and_price.py --csv data/raw/SPY.csv ...` |
| Roll-forward | `PYTHONPATH=src python scripts/run_phase2_roll_forward.py ...` → `output/phase2_roll_forward.csv` |
| Plot roll-forward | `python scripts/plot_phase2_roll_forward.py` — if `visuals/phase2_roll_forward_panel.png` already exists, the previous file is copied to `visuals/archive/` with a UTC timestamp before overwrite. Same for `plot_spy_volatility_regime.py` → `visuals/archive/`. |
| SPY vol figure | `python scripts/plot_spy_volatility_regime.py` (as configured in your workflow) |
| Quotes scaffold | Fill `data/raw/quotes_template.csv`, run `compare_mc_to_quotes.py` |

### D. Course-methods alignment (PETSc Vec, KSP/PC, SNES)

Copy-paste **commands** for these scripts live in **§4.2** (PETSc-heavy course bridge).

| Script | Purpose |
| --- | --- |
| [`demo_petsc_vec_global_sum.py`](../scripts/demo_petsc_vec_global_sum.py) | Parallel **`PETSc.Vec`**, one row per rank, `Vec.sum`. |
| [`implied_vol_snes.py`](../scripts/implied_vol_snes.py) | Scalar **`SNES`**: \(\mathrm{BS}(\sigma) = V_{\mathrm{mkt}}\); optional `--analytic-jacobian`. |
| [`run_poisson_pc_sweep.py`](../scripts/run_poisson_pc_sweep.py) | Poisson **KSP + PC** sweep → `output/poisson_pc_sweep.csv`. |
| [`run_poisson_ksp_sweep.py`](../scripts/run_poisson_ksp_sweep.py) | Poisson **KSP-type** sweep (fixed PC) → `output/poisson_ksp_sweep.csv`. |
| [`plot_poisson_ksp_sweep.py`](../scripts/plot_poisson_ksp_sweep.py) | Bar charts from `output/poisson_ksp_sweep.csv`. |
| [`run_sampling_vs_solver_benchmark.py`](../scripts/run_sampling_vs_solver_benchmark.py) | **One Poisson solve vs one MC** timing → `output/sampling_vs_solver_benchmark.csv`. |
| [`run_poisson_timed.py`](../scripts/run_poisson_timed.py), [`run_poisson_convergence.py`](../scripts/run_poisson_convergence.py) | Structured-grid Poisson (Mat/Vec/KSP); see [`PHASE3_RUNBOOK.md`](PHASE3_RUNBOOK.md). |
| [`plot_poisson_pc_sweep.py`](../scripts/plot_poisson_pc_sweep.py) | Bar charts from `output/poisson_pc_sweep.csv` (time + iterations per PC). |
| [`plot_diffusion_solver_history.py`](../scripts/plot_diffusion_solver_history.py) | Iterations + residual vs pseudo-time from `heat_solver_history.csv` / advection logs. |
| [`run_heat_paraview.py`](../scripts/run_heat_paraview.py), [`run_advection_diffusion_paraview.py`](../scripts/run_advection_diffusion_paraview.py) | `TS` + ParaView exports; see **§4.2** and [`PROJECT_OUTLINE.md` Appendix A](../PROJECT_OUTLINE.md#appendix-a-paraview-and-vtk-workflows). |

### E. Provenance and long jobs

- After meaningful runs, append the registry: **`python scripts/log_experiment_run.py`** (see [`HPC_EXPERIMENT_REGISTRY.md`](HPC_EXPERIMENT_REGISTRY.md)).
- Long shells: [`LONG_RUNNING_JOBS.md`](LONG_RUNNING_JOBS.md), **`scripts/terminal_probe.py`**.
- Optional cluster later: [`RUN_LATER_CLUSTER_APPTAINER.md`](RUN_LATER_CLUSTER_APPTAINER.md), **`scripts/submit_scaling.slurm`**.

### F. ParaView (optional)

- Path bundles: **`export_paraview_paths.py`**, **`export_stock_paths_vtp.py`**.
- ParaView / VTK: [`PROJECT_OUTLINE.md` Appendix A](../PROJECT_OUTLINE.md#appendix-a-paraview-and-vtk-workflows) (stub [`PARAVIEW_GUIDE_PETSC.md`](../PARAVIEW_GUIDE_PETSC.md)).

---

## 6. “Still to do” once PETSc works (suggested priority)

1. **`check_petsc_env.py` passes** with your real `PYTHONPATH` / `PATH`.
2. **Refresh baseline CSVs/PNGs** for the report (convergence, scaling, barrier, error surface as needed).
3. **Run course-bridge scripts** (§5.D) once; add one short subsection to the report citing iteration counts / timing ratio.
4. **Phase 2** on your chosen ticker window; update `RESULTS_PHASE1.md` or a dedicated results doc if numbers move.
5. **`log_experiment_run.py`** for each machine/rank configuration you show in the final PDF.
6. Optional: **cluster** scaling only if course asks; otherwise local `mpiexec` is enough if documented.

---

## 7. Quick reference: module ↔ responsibility

| Module | Responsibility |
| --- | --- |
| `gbm.py` | GBM path simulation; antithetic pairs. |
| `payoffs.py` | European / Asian / barrier payoffs. |
| `engine.py` | `petsc_monte_carlo_price`, reductions, `petsc_asian_call_cv_european_control`. |
| `black_scholes.py` | Closed-form call + **vega** (for SNES Jacobian). |
| `historical_data.py` | CSV load (closes / adj close + dates). |
| `vol_calibration.py` | Annualized realized vol from closes. |

---

## 8. PETSc options from the shell (no code edits)

PETSc reads **`PETSC_OPTIONS`** (space-separated flags). Useful for demos in the report (residual history, why a PC failed).

**KSP / PC (Poisson, heat, advection drivers):**

```bash
export PETSC_OPTIONS="-ksp_monitor -ksp_converged_reason -ksp_rtol 1e-8"
mpiexec -n 4 python -u scripts/run_poisson_timed.py --nx 65 --ny 65 --ksp-type cg --pc-type jacobi
```

**SNES (implied vol script uses prefix `iv_`):**

```bash
export PETSC_OPTIONS="-iv_snes_monitor -iv_snes_converged_reason"
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
mpiexec -n 1 python -u scripts/implied_vol_snes.py --market-price 10.4506 --analytic-jacobian
```

Clear `PETSC_OPTIONS` when you do not want extra logging (`unset PETSC_OPTIONS`).

---

## 9. One-liners that save time

**Black–Scholes price for a synthetic “mid” (implied-vol sanity):**

```bash
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
python -c "from montecarlo.black_scholes import black_scholes_call as c; print(c(100,100,0.05,0.2,1.0))"
```

**Asian / barrier smokes under MPI:**

```bash
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
mpiexec -n 4 python -u scripts/run_pricing.py --option asian --paths 80000
mpiexec -n 4 python -u scripts/run_pricing.py --option barrier --paths 80000 --barrier 130
```

**Batch scaling across options (manifest under `output/experiments/`):**

```bash
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
mpiexec -n 4 python -u scripts/run_experiment_matrix.py --help
```

**ParaView smoke (small Poisson export):**

```bash
mpiexec -n 1 python -u scripts/smoke_paraview_export.py
```

---

## 10. Registry: log what you will cite

After a run, append a row so the PDF matches machine evidence:

```bash
python scripts/log_experiment_run.py \
  --key sampling_vs_solver_p4 \
  --ranks 4 \
  --output output/sampling_vs_solver_benchmark.csv \
  --notes "Poisson 129² CG+Jacobi vs European MC 80k paths" \
  --cmd 'mpiexec -n 4 python -u scripts/run_sampling_vs_solver_benchmark.py --nx 129 --ny 129 --paths 80000'
```

Use a fresh **`--key`** per distinct configuration; keep **`--output`** paths honest.

---

## 11. Optional “deeper PETSc” (not in repo yet)

If you want **more** PETSc surface area later (each is a small project):

| Idea | PETSc objects | Benefit |
| --- | --- | --- |
| Path statistics in **`Vec`** | `Vec` per rank, `VecDot` / `VecSum` | Same MC math, explicit Vec API in `engine.py`. |
| **`PETSC_OPTIONS=-log_view` + `MC_PETSC_LOG_EVENTS`** | `PETSc.Log.Event` in `engine.py` | **Implemented (opt-in):** pathgen vs payoff vs `allreduce` time inside `petsc_monte_carlo_price` / Asian CV; see **§12**. |
| **`TS`** for deterministic \(dS=rS\,dt\) | `TS` + RHS | Clean time-stepping lecture tie-in (no noise). |
| **GAMG / fieldsplit** on Poisson | `PC` options | Stronger preconditioner narrative (if your build supports it). |

---

## 12. PETSc profiling inside Monte Carlo (opt-in)

`src/montecarlo/engine.py` registers **`PETSc.Log.Event`** regions when:

```bash
export MC_PETSC_LOG_EVENTS=1
```

Event names:

| Event | Where |
| --- | --- |
| `MC_paths` | GBM path simulation |
| `MC_payoff` | Payoff + discount |
| `MC_reduce` | Local sums + MPI `allreduce` |
| `MC_cv_paths` | Asian CV path simulation |
| `MC_cv_local` | Build \(X,Y\) and local sums |
| `MC_cv_reduce_xy` | First reduction block |
| `MC_cv_Z` | Form \(Z\) and local sums |
| `MC_cv_reduce_Z` | Final reduction on \(Z\) |

To print a summary at process exit (typical one-shot):

```bash
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
export MC_PETSC_LOG_EVENTS=1
export PETSC_OPTIONS="-log_view"
mpiexec -n 4 python -u scripts/validate_european.py --paths 200000
```

Unset `MC_PETSC_LOG_EVENTS` (or set to `0`) for zero extra work in production batches. If `Event.begin` errors on your build, drop `MC_PETSC_LOG_EVENTS` and rely on **`-log_view`** only for KSP-heavy scripts.

---

## 13. Troubleshooting (quick)

| Symptom | Things to check |
| --- | --- |
| `ModuleNotFoundError: petsc4py` | Use the **course** Python; set **`PYTHONPATH`** to `$PETSC_DIR/$PETSC_ARCH/lib` before `src`. |
| `ModuleNotFoundError: montecarlo` | From repo root, include **`src`** in **`PYTHONPATH`**. |
| `mpiexec` uses wrong Python | Call **`mpiexec /full/path/to/python script.py`** or `export PATH=.../apma4302-pkgs-opt/bin:$PATH`. |
| Rank hangs / mismatch | Same script on all ranks; avoid rank-0-only `input()`; check collective **order** (MC events are collective when logging is on). |
| First run “idle” for minutes | Cold **dyld** load; use **`DYLD_LIBRARY_PATH`**, **`python -u`**, or wait once per fresh process. |
| Poisson / KSP diverges or wrong PC | Try **`--ksp-type cg --pc-type jacobi`** on Laplacian; avoid parallel **`ilu`** unless your PETSc supports it. |

---

*Living doc: extend §5–§13 as you add figures, cluster runs, or new drivers.*
