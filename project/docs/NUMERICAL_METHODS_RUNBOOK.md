# Numerical methods runbook (course bridge)

Single map from **method** → **math idea** → **script** → **typical output**.  
Use with [`docs/PROJECT_SUMMARY_AND_PETSC_CHECKLIST.md`](PROJECT_SUMMARY_AND_PETSC_CHECKLIST.md) §4.2 and [`docs/HPC_EXPERIMENT_REGISTRY.md`](HPC_EXPERIMENT_REGISTRY.md).

**Environment (every PETSc command):**

```bash
export PETSC_DIR=/Users/dan/Desktop/Columbia/HPC_4302/petsc
export PETSC_ARCH=apma4302-pkgs-opt
export PYTHONPATH="$PETSC_DIR/$PETSC_ARCH/lib:src"
export DYLD_LIBRARY_PATH="$PETSC_DIR/$PETSC_ARCH/lib:${DYLD_LIBRARY_PATH:-}"   # macOS optional
```

Use **`python -u`** on first import of the day; **`mpiexec`** from the same MPI stack PETSc was built against.

---

## 1) Monte Carlo and stochastic numerics

| Method | What it is | Script(s) | Output |
| --- | --- | --- | --- |
| **GBM paths** | Euler–Maruyama / exact-step on log-price; antithetic pairs | `src/montecarlo/gbm.py`; drivers use `engine.py` | paths → payoffs |
| **Weak convergence \(O(N^{-1/2})\)** | MC error vs number of paths | `scripts/run_convergence.py`, `plot_convergence.py` | `output/convergence*.csv`, PNG |
| **Variance reduction (antithetic)** | Correlated paths, lower variance at fixed \(N\) | `--antithetic` on MC drivers; `run_variance_reduction_bench.py` | bench CSV |
| **Control variates** | Asian + European control | `scripts/run_asian_cv_control_variate.py` | stdout / CSV per `--help` |
| **Barrier / discontinuous payoffs** | Pathwise MC; variance vs geometry | `run_barrier_sweep.py`, `plot_barrier_sweep.py` | `output/barrier_sweep.csv` |
| **Error surface \((N,\sigma)\)** | MC vs BS grid | `run_error_surface.py`, `plot_error_surface.py` | `output/error_surface_european.csv` |
| **Black–Scholes / Greeks** | Closed-form + **vega** for SNES Jacobian | `src/montecarlo/black_scholes.py` | used by `implied_vol_snes.py` |

---

## 2) Nonlinear equations (Newton / SNES)

| Method | What it is | Script | Output |
| --- | --- | --- | --- |
| **Scalar root / Newton** | Solve \(\mathrm{BS}(\sigma)=V_{\mathrm{mkt}}\) in \(\sigma\) | `scripts/implied_vol_snes.py` | stdout; optional `--out` JSON |
| **Analytic Jacobian** | \(1\times 1\) Jacobian = vega | `--analytic-jacobian` | fewer SNES matvec probes |

```bash
mpiexec -n 1 python -u scripts/implied_vol_snes.py --market-price 10.4506 --analytic-jacobian
```

---

## 3) Linear algebra on partitioned data (`Vec`)

| Method | What it is | Script | Output |
| --- | --- | --- | --- |
| **Global reduction on a `Vec`** | Partitioned vector + `Vec.sum` (same *combine* pattern as inner products in Krylov) | `scripts/demo_petsc_vec_global_sum.py` | stdout |

```bash
mpiexec -n 4 python -u scripts/demo_petsc_vec_global_sum.py
```

---

## 4) Elliptic PDE: Poisson (FD + Krylov + PC)

| Method | What it is | Script(s) | Output |
| --- | --- | --- | --- |
| **5-point Laplacian** | Sparse `Mat` + `Vec` RHS | `run_poisson_timed.py` (build used by sweeps) | append timed CSV |
| **Spatial convergence** | \(L^2\) error vs mesh (manufactured / discrete norm as implemented) | `run_poisson_convergence.py`, `plot_phase3.py` | `output/poisson_convergence.csv`, `visuals/spatial_convergence_loglog.png` |
| **PC comparison** | Fixed KSP, sweep preconditioners | `run_poisson_pc_sweep.py`, `plot_poisson_pc_sweep.py` | `output/poisson_pc_sweep.csv`, PNGs |
| **KSP comparison** (new) | Fixed PC, sweep **KSP** types (CG vs GMRES vs BiCGStab, etc.) | `run_poisson_ksp_sweep.py`, `plot_poisson_ksp_sweep.py` | `output/poisson_ksp_sweep.csv`, PNGs |

```bash
mpiexec -n 4 python -u scripts/run_poisson_ksp_sweep.py --nx 65 --ny 65 --pc-type jacobi \
  --ksp-types cg gmres bcgsl --out output/poisson_ksp_sweep.csv
python scripts/plot_poisson_ksp_sweep.py --infile output/poisson_ksp_sweep.csv --out-prefix visuals/poisson_ksp_sweep
```

**Note:** Poisson with Dirichlet BCs yields a **symmetric definite**-structured system; **CG** is the natural method; **GMRES / BiCGStab** illustrate general Krylov templates and iteration counts may differ. On **coarse** grids every method may converge in a **single** iteration; use larger `--nx/--ny` (e.g. 129 or 257) if you want iteration counts to separate on plots.

---

## 5) Parabolic / advection–diffusion (time marching + implicit solves)

| Method | What it is | Script | Output |
| --- | --- | --- | --- |
| **Backward Euler + implicit step** | Heat: solve \((I+\Delta t L)u^{n+1}=u^n\) each step | `run_heat_paraview.py` | `output/paraview/vtp/heat.pvd`, `heat_solver_history.csv` |
| **Upwind-ish advection + diffusion** | Nonsymmetric sparse system per step; defaults **cg+jacobi** (parallel-safe) | `run_advection_diffusion_paraview.py` | `output/paraview/vtp/advection_diffusion.pvd` (and `advection_diffusion_*.vtp`) |
| **Solver history plots** | KSP iterations / residual vs pseudo-time | `plot_diffusion_solver_history.py`, `make_phase2_figures.py` | PNGs under `visuals/` |

ParaView workflow: [`PROJECT_OUTLINE.md` Appendix A](../PROJECT_OUTLINE.md#appendix-a-paraview-and-vtk-workflows) (stub [`PARAVIEW_GUIDE_PETSC.md`](../PARAVIEW_GUIDE_PETSC.md)).

---

## 6) “Two cultures”: sampling cost vs linear solve cost

| Method | What it is | Script | Output |
| --- | --- | --- | --- |
| **One KSP vs one MC** | Same `mpiexec` width; compare wall times | `run_sampling_vs_solver_benchmark.py` | `output/sampling_vs_solver_benchmark.csv` |

---

## 7) Suggested order for a “numerical methods” results pass

1. `run_poisson_convergence.py` + `plot_phase3.py --conv …`  
2. `run_poisson_ksp_sweep.py` + `plot_poisson_ksp_sweep.py`  
3. `run_poisson_pc_sweep.py` + `plot_poisson_pc_sweep.py`  
4. `implied_vol_snes.py` (synthetic mid)  
5. `run_sampling_vs_solver_benchmark.py`  
6. `run_heat_paraview.py` + `run_advection_diffusion_paraview.py` + ParaView exports  
7. `run_variance_reduction_bench.py` and/or `run_asian_cv_control_variate.py`  
8. Append `output/registry/experiment_runs.csv` via `log_experiment_run.py`

---

## 8) Syllabus one-liners (for prose)

- **FD + stability:** heat and advection–diffusion use implicit steps (linear solves); MC uses explicit path simulation with statistical error decay.  
- **Krylov:** Poisson drivers expose **KSP** iteration counts; KSP sweep and PC sweep make the **outer iteration vs preconditioner** tradeoff visible.  
- **Nonlinear:** implied vol is a **1D SNES** problem with optional **hand Jacobian** (vega).  
- **HPC:** same PETSc world communicator links **parallel MC** (path batches) to **parallel sparse solves**.
