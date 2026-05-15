# HPC experiment registry (living document)

This file is the **single place** to see what experiments exist, which **class-style HPC ideas** they exercise, where outputs go, and what has been run. Keep it in sync with `PROJECT_OUTLINE.md` (plan and math) and `RESULTS_PHASE1.md` (narrative results snapshots).

## How to log a new run (append-only)

Machine-readable history lives in `output/registry/experiment_runs.csv`. After you finish a batch of jobs, append a row:

```bash
python scripts/log_experiment_run.py \
  --key strong_european_p8 \
  --ranks 8 \
  --output output/scaling_results.csv \
  --notes "laptop PETSc; strong scaling baseline" \
  --cmd 'mpirun -n 8 python scripts/run_pricing.py ...'
```

If you prefer, edit `experiment_runs.csv` by hand (keep the header row). The script fills `hostname` and UTC `timestamp` automatically.

## Class HPC topics ↔ this project

| Typical course theme | Where it shows up here |
| --- | --- |
| MPI / SPMD rank parallelism | `mpirun` drivers; `petsc4py` + `mpi4py` in `src/montecarlo/engine.py` |
| Domain decomposition (embarrassingly parallel “domains” = path batches) | `_split_paths` in `engine.py`; each rank owns disjoint paths |
| Collective communication | `allreduce` for global sums of payoff, payoff², and path counts |
| Strong scaling (fixed problem size, grow \(p\)) | `scripts/run_scaling.py --mode strong` |
| Weak scaling (fixed work per rank, grow \(p\)) | `scripts/run_scaling.py --mode weak` |
| Speedup \(S(p)\), parallel efficiency \(E(p)\) | CSV columns from `run_scaling.py`; `scripts/plot_scaling.py` |
| Serial / parallel time mental model (Amdahl intuition) | Reduction + launcher overhead visible when paths/rank shrink at large \(p\) |
| Load balance | Even path counts via integer split; barrier payoff changes variance, not step count in this code |
| Reproducibility / RNG streams | Per-rank seeds `seed + 10007 * rank` in `petsc_monte_carlo_price` |
| Numerical method + HPC coupling | MC error \(O(N^{-1/2})\) vs wall time from parallelism (`run_convergence.py`) |
| Variance reduction as “algorithmic efficiency” | Antithetic pairs in `gbm.py`; bench `run_variance_reduction_bench.py` |
| **PETSc `Vec` (partitioned)** | `scripts/demo_petsc_vec_global_sum.py` — one row per rank, `Vec.sum` |
| **SNES / Newton (scalar)** | `scripts/implied_vol_snes.py` — \(\mathrm{BS}(\sigma)=\) market; optional analytic Jacobian (`black_scholes_call_vega`) |
| **KSP + PC on structured Poisson** | `scripts/run_poisson_timed.py`, `scripts/run_poisson_convergence.py`, `scripts/run_poisson_pc_sweep.py`, `scripts/run_poisson_ksp_sweep.py` |
| **Sampling vs solver wall time** | `scripts/run_sampling_vs_solver_benchmark.py` — one Poisson solve vs one European MC call |
| **PETSc profiling (`Log.Event`)** | `MC_PETSC_LOG_EVENTS=1` → `engine.py` marks `MC_paths` / `MC_payoff` / `MC_reduce` (and CV variants); use with `PETSC_OPTIONS=-log_view` |

### Spring syllabus weeks ↔ this repo (high level)

Use this table when you want the **final write-up** to echo **specific lectures** without overstating PDE overlap.

| Block (approx.) | Course focus | Honest link here |
| --- | --- | --- |
| Late Jan | HPC + PETSc intro | MPI + `petsc4py` world comm, `mpirun` drivers |
| Early Feb | Parallel LA / Krylov | SPMD + **`PETSc.Vec`** demo (`demo_petsc_vec_global_sum.py`); MC **`allreduce`** same *combine* pattern as Krylov dots |
| Mid Feb | Sparse direct | Not used; optional one-line contrast with PDE labs |
| Late Feb–Mar | PDE + Newton / SNES | Poisson drivers + **`implied_vol_snes.py`** (scalar **`SNES`**) |
| Mid Mar | Time stepping (TS) | Discrete path simulation vs lecture TS themes (explicit marching; different stability) |
| Late Mar–Apr | Preconditioners / MG | **`run_poisson_pc_sweep.py`** / **`run_poisson_ksp_sweep.py`** for PC vs KSP iteration/time; contrast MC sampling bottleneck |
| Late Apr | **Scaling & performance** | Primary alignment: strong/weak scaling, speedup, efficiency |

## Experiment catalog (definitions)

Status legend: **done** = baseline implemented and/or numbers in `RESULTS_PHASE1.md`; **ready** = script exists, fill in runs when you execute; **planned** = backlog in `PROJECT_OUTLINE.md` §12.

| Key | Script(s) | What it measures | Primary HPC methods | Default / typical outputs |
| --- | --- | --- | --- | --- |
| `price_smoke` | `run_pricing.py` | Single MC price + CI | MPI path batching, reductions | stdout |
| `validate_bs` | `validate_european.py` | MC vs Black–Scholes | Same + verification | stdout |
| `convergence_n` | `run_convergence.py`, `plot_convergence.py` | Error vs \(N\), log–log slope | Statistical scaling + MPI | `output/convergence.csv`, PNG |
| `scaling_strong` | `run_scaling.py`, `plot_scaling.py` | Strong scaling | Speedup, efficiency | `output/scaling_results.csv`, PNG |
| `scaling_weak` | `run_scaling.py`, `plot_scaling.py` | Weak scaling | Same | same pattern |
| `scaling_matrix` | `run_experiment_matrix.py` | Cross-option (+ optional antithetic) scaling | Batch strong/weak studies | `output/experiments/*.csv`, `experiment_matrix_manifest.json` |
| `scaling_compare_plot` | `plot_scaling_compare.py` | Overlay speedup curves | Comparative performance reporting | PNG path via `--out` |
| `variance_antithetic` | `run_variance_reduction_bench.py` | stderr / CI width plain vs antithetic | Variance vs wall time at fixed \(N\) | `output/variance_reduction_bench.csv` |
| `paraview_paths` | `export_paraview_paths.py` | 3D path bundle export | I/O for visualization | `output/path_bundle.csv`, etc. |
| `cluster_slurm` | `submit_scaling.slurm` | Optional batch scaling on a shared system | Job scripts, queueing | Slurm logs; **defer setup** to [`RUN_LATER_CLUSTER_APPTAINER.md`](RUN_LATER_CLUSTER_APPTAINER.md) |
| `barrier_sweep` | `run_barrier_sweep.py`, `plot_barrier_sweep.py` | Price/stderr/CI vs barrier \(B\) | Variance / payoff geometry vs MC workload | `output/barrier_sweep.csv`, PNG |
| `error_surface` | `run_error_surface.py`, `plot_error_surface.py` | \(|MC-BS|\) on \((N,\sigma)\) grid | Accuracy geography + wall time per grid point | `output/error_surface_european.csv`, PNG |
| `yahoo_ohlcv_fetch` | `fetch_yahoo_ohlcv.py` | Download daily OHLCV | External public data | `data/raw/*.csv` (needs `requirements-phase2.txt`) |
| `phase2_vol_mc` | `phase2_calibrate_and_price.py` | Calibrated \(\sigma\) vs baseline MC + BS | Same MC engine; CSV-driven parameters | stdout, optional `--out-json` |
| `hpc_gallery` | `render_hpc_figures.sh`, `plot_hpc_*` | Dashboard-style scaling / convergence / barrier / SPY vol | Reporting | `visuals/hpc_*.png`, `visuals/phase2_spy_price_and_vol.png` |
| `petsc_vec_sum` | `demo_petsc_vec_global_sum.py` | Partitioned `Vec` + `Vec.sum` vs hand reduction | Parallel **Vec** layout (LA week) | stdout |
| `implied_vol_snes` | `implied_vol_snes.py` | Scalar implied \(\sigma\); SNES iters; FD vs analytic Jacobian | **SNES** / Newton | stdout, optional JSON |
| `sampling_vs_solver` | `run_sampling_vs_solver_benchmark.py` | Poisson KSP time vs MC time on same ranks | **KSP+PC** vs **Monte Carlo** bottleneck | appends `output/sampling_vs_solver_benchmark.csv` |
| `poisson_pc_sweep` | `run_poisson_pc_sweep.py` | Iterations / time vs `--pc-types` | Preconditioned Krylov | appends `output/poisson_pc_sweep.csv` |
| `poisson_ksp_sweep` | `run_poisson_ksp_sweep.py`, `plot_poisson_ksp_sweep.py` | Iterations / time vs **KSP** type (fixed PC) | Krylov family comparison | `output/poisson_ksp_sweep.csv`, PNGs |
| `mc_petsc_log` | `validate_european.py` / any MC driver | PETSc `Log.Event` breakdown path vs payoff vs reduce | Set `MC_PETSC_LOG_EVENTS=1`, `PETSC_OPTIONS=-log_view` | stdout at exit |
| `hpc_speedup_ideal` | `plot_hpc_speedup_vs_ideal.py` | Strong speedup vs ideal line + efficiency | Reads `output/strong_scaling.csv` | `visuals/hpc_speedup_vs_ideal.png` |
| `hpc_weak_walltime` | `plot_hpc_weak_scaling_walltime.py` | Weak scaling wall time flatness | Reads `output/weak_scaling.csv` | `visuals/hpc_weak_scaling_walltime.png` |
| `mpi_domain_vtp` | `export_mpi_domain_paths_vtp.py` | ParaView polylines with `mpi_rank` field | Synthetic SPMD batches | `output/paraview/vtp/mpi_domain_paths.vtp` |
| `registry_append` | `log_experiment_run.py` | Bookkeeping / provenance | Traceability | appends `output/registry/experiment_runs.csv` |

## Run history (high-level)

Detailed rows: **`output/registry/experiment_runs.csv`**.  
Archived narrative snapshots: **`RESULTS_PHASE1.md`**.  
Batch matrix manifests: **`output/experiments/experiment_matrix_manifest.json`** (after you run the matrix driver).

### Baseline rows (seeded with Phase 1 summary)

These mirror the numbers already documented in `RESULTS_PHASE1.md` so the CSV is never “empty” on a fresh clone.

| When (UTC) | Key | Ranks | Outputs | Notes |
| --- | --- | --- | --- | --- |
| (see CSV) | `results_phase1_doc` | 2–8 | `RESULTS_PHASE1.md`, `output/*.csv` | Documented baseline; re-run scripts to refresh numbers |

When you add new clusters, new option types, or antithetic scaling, **append** to `experiment_runs.csv` rather than deleting old rows.

## Maintenance checklist (quick)

- [ ] Prefer **local PETSc** runs for all figures in the main report; use `RUN_LATER_CLUSTER_APPTAINER.md` only if you add a small optional cluster appendix.
- [ ] After a scaling campaign (local or Slurm): new CSV under `output/` or `output/experiments/`, new PNG under `visuals/` or `output/`, one `log_experiment_run.py` line per distinct configuration.
- [ ] Long or background shells: run `python scripts/terminal_probe.py --auto` twice (~60s apart); confirm transcript or file size moves (`docs/LONG_RUNNING_JOBS.md`).
- [ ] After changing code that affects numbers: update `RESULTS_PHASE1.md` section headers or add “Revision” subsection; bump `PROJECT_OUTLINE.md` status bullets.
- [ ] Before submitting report: confirm every figure path in the report appears in `experiment_runs.csv` or this catalog table.
