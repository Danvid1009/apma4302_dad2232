# Project Outline: Parallel Monte Carlo Option Pricing with ParaView Visualization

## 0. Current Project Status (Updated)

This outline now reflects both the project plan and implementation progress.

**Completed foundation**
- Core GBM path simulator implemented
- Payoff modules implemented for European, Asian, and barrier up-and-out options
- Parallel PETSc-based Monte Carlo engine implemented with global reductions for:
  - payoff sum
  - payoff squared sum
  - total path count
- 95% confidence interval computation implemented
- European closed-form Black-Scholes function implemented
- Validation script (MC vs Black-Scholes) implemented
- ParaView path and barrier helper export script implemented
- Scaling automation scripts implemented (strong/weak)
- Slurm batch template added for SSH/HPC runs

**Current code map**
- Core library: `src/montecarlo/`
  - `gbm.py`
  - `payoffs.py`
  - `engine.py`
  - `black_scholes.py`
- Experiment scripts: `scripts/`
  - `run_pricing.py`
  - `validate_european.py`
  - `export_paraview_paths.py`
  - `run_scaling.py`
  - `plot_scaling.py`
  - `submit_scaling.slurm`

**Next priority tasks**
- Generate final figures/tables for report
- Add optional variance-reduction methods (antithetic/control variate)
- Add realism extension (historical calibration + market quote comparison)

**Project execution order**
- Phase 1 (current): methods-first baseline using synthetic/model inputs
- Phase 2 (after baseline): realism extension with market-calibrated inputs and market quote comparison

## 1. Problem Statement

This project prices financial derivatives using Monte Carlo simulation and studies how parallel computing improves performance while maintaining numerical accuracy.

**Target instruments**
- European call option
- Asian option
- Barrier option

**Core research question**
- How does parallel Monte Carlo reduce runtime while preserving pricing accuracy?

## 2. Mathematical Model

Assume the stock price follows geometric Brownian motion (GBM):

\[
dS_t = rS_t\,dt + \sigma S_t\,dW_t
\]

Exact simulation step:

\[
S_{t+\Delta t}
=
S_t
\exp\left[
\left(r-\frac{1}{2}\sigma^2\right)\Delta t
+
\sigma \sqrt{\Delta t} Z
\right], \quad Z\sim \mathcal{N}(0,1)
\]

## 3. Option Payoffs

European call:

\[
P = e^{-rT}\max(S_T-K,0)
\]

Asian call:

\[
P = e^{-rT}\max(\bar S-K,0)
\]

Barrier call (up-and-out style):

\[
P =
\begin{cases}
e^{-rT}\max(S_T-K,0), & S_t < B \text{ for all } t \\
0, & \text{otherwise}
\end{cases}
\]

## 4. Parallel Monte Carlo Method

Each PETSc/MPI rank simulates an independent batch of paths.

For rank \(p\):

\[
\hat V_p = \frac{1}{N_p}\sum_{i=1}^{N_p} P_i
\]

Global estimate:

\[
\hat V =
\frac{1}{N}
\sum_p
\sum_{i=1}^{N_p} P_i
\]

Use PETSc/MPI reduction to combine:
- payoff sums
- payoff squared sums
- total number of paths

Then compute confidence intervals:

\[
\hat V \pm 1.96\frac{\hat\sigma}{\sqrt{N}}
\]

## 5. Validation

For the European call, compare Monte Carlo against the closed-form Black-Scholes value:

\[
C = S_0\Phi(d_1)-Ke^{-rT}\Phi(d_2)
\]

Report:
- absolute error
- relative error
- convergence rate

Expected convergence:

\[
\text{error} = O(N^{-1/2})
\]

**Implementation status**
- Implemented: `scripts/validate_european.py` prints
  - MC price
  - BS price
  - absolute error
  - relative error
  - standard error
  - 95% confidence interval
- Completed baseline run (\(N=400{,}000\), ranks=4):
  - \(V_{MC}=10.42161941\)
  - \(V_{BS}=10.45058357\)
  - absolute error \(=2.8964\times 10^{-2}\)
  - relative error \(=2.7715\times 10^{-3}\)
  - 95% CI \([10.37609515,\;10.46714366]\)
- Convergence sweep implemented:
  - `scripts/run_convergence.py`
  - `scripts/plot_convergence.py`
- Convergence evidence (paths \(2\times 10^4\) to \(8\times 10^5\), reps=5, ranks=4):
  - fitted slope (mean-abs error) \(=-0.547067\)
  - fitted slope (RMSE) \(=-0.529293\)
  - consistent with expected \(O(N^{-1/2})\)

## 6. HPC Experiments

### Strong Scaling

Fix total paths \(N\), increase processor count \(p\).

Measure:

\[
S(p)=\frac{T(1)}{T(p)}, \qquad
E(p)=\frac{S(p)}{p}
\]

**Implementation status**
- Implemented driver: `scripts/run_scaling.py --mode strong`
- Implemented plotting: `scripts/plot_scaling.py`
- Completed baseline run (European, total paths \(=400{,}000\)):
  - \(p=1:\;T=4.8319s,\;S=1.000,\;E=1.000\)
  - \(p=2:\;T=3.1022s,\;S=1.558,\;E=0.779\)
  - \(p=4:\;T=2.8031s,\;S=1.724,\;E=0.431\)
  - \(p=8:\;T=2.4953s,\;S=1.936,\;E=0.242\)

### Weak Scaling

Fix paths per processor and increase \(p\). Measure runtime stability.

**Implementation status**
- Implemented driver: `scripts/run_scaling.py --mode weak`
- Runtime, speedup, and efficiency exported to CSV for reporting
- Completed baseline run (European, \(100{,}000\) paths per rank):
  - \(p=1:\;T=1.3559s\)
  - \(p=2:\;T=1.7754s\)
  - \(p=4:\;T=2.9858s\)
  - \(p=8:\;T=5.4020s\)
  - observed runtime growth indicates local-machine overhead and weak-scaling inefficiency at higher \(p\)

### Variance Reduction (optional extensions)
- Antithetic variates
- Control variates (using Black-Scholes European call)
- Quasi-Monte Carlo

Compare methods using error per unit runtime.

**Implementation status**
- Not yet implemented (planned extension phase)

### Cluster/SSH Execution Workflow

Run experiments either interactively over SSH or by Slurm using a PETSc-enabled Python environment.

- Interactive:
  - `python scripts/run_scaling.py --mode strong ...`
  - `python scripts/run_scaling.py --mode weak ...`
- Batch:
  - `sbatch scripts/submit_scaling.slurm`

## 6.5 Methods-First Execution Plan (Phase 1)

This phase is the main technical core and should be completed before realism extensions.

### Step A: Numerical correctness baseline
- Run `scripts/run_pricing.py` for European, Asian, and barrier options.
- Run `scripts/validate_european.py` and record:
  - MC vs BS absolute/relative error
  - confidence interval width
- Add convergence sweep over increasing \(N\), estimate slope in log-log error plot.
- Status: completed.

### Step B: HPC performance study
- Strong scaling: fixed total paths, vary ranks.
- Weak scaling: fixed paths per rank, vary ranks.
- Generate speedup/efficiency plots from `scripts/plot_scaling.py`.
- Record practical bottlenecks (MPI overhead, RNG cost, memory traffic).
- Status: completed baseline runs and plots on local PETSc setup.

### Step C: Visualization study
- Export path bundles and barrier helper geometry from `scripts/export_paraview_paths.py`.
- Build ParaView figures for:
  - path ensemble geometry
  - barrier crossing behavior
  - (later) error surface once convergence grid is available
- Status: data export completed; ParaView figure assembly in progress.

### Step D: Report-ready outputs
- Save all CSVs/plots in `output/`.
- Build a concise methods/results narrative:
  - numerical accuracy
  - scaling performance
  - visualization interpretation
- Status: partially completed (outputs generated; final narrative/tables pending).

## 7. ParaView Visualizations

### Visualization 1: Monte Carlo Path Bundle

Export simulated stock paths as 3D curves:
\[
(t, S_t, i)
\]
where \(i\) is path index.

Use ParaView tube rendering to display path ensembles.

**Implementation status**
- Implemented exporter: `scripts/export_paraview_paths.py`
- Output format: CSV with `(path_id, time, stock_price)`

### Visualization 2: Barrier Option Knockout Geometry

Display:
- surviving paths
- knocked-out paths
- barrier level as a plane

This gives a geometric view of barrier crossings.

**Implementation status**
- Barrier helper geometry export implemented (barrier plane CSV)
- Path classification coloring (survive/knockout) is next step in post-processing

### Visualization 3: Error Surface

Compute error on parameter grids such as:
- \((N,\sigma)\), or
- \((S_0,\sigma)\)

Export structured grid with:
\[
z = |\hat V - V_{\text{BS}}|
\]
and visualize the convergence/error landscape.

**Implementation status**
- Not yet implemented (planned after convergence/scaling baseline is complete)

## 8. Expected Challenges

- Slow Monte Carlo convergence
- Random number stream management across MPI ranks
- Load imbalance for barrier options with early path termination
- Communication overhead at high processor counts
- ParaView export format and debugging issues

## 9. Deliverables

- Parallel Monte Carlo implementation
- European, Asian, and barrier option pricing results
- Validation against Black-Scholes
- Strong and weak scaling plots
- Confidence interval analysis
- ParaView visualizations
- Discussion of performance bottlenecks and failure modes

**Deliverable tracking**
- Completed: parallel PETSc Monte Carlo implementation and baseline validation infrastructure
- Completed: scaling automation and plotting infrastructure
- Completed: convergence sweep tooling and evidence of \(O(N^{-1/2})\)-consistent behavior
- In progress: final figures/tables and report narrative
- Pending: realism extension (historical calibration + market quote comparison)

## 10. Realism Extension (Phase 2, after Methods Baseline)

After Phase 1 is complete, add a realism section with empirical market inputs.

### 10.1 Calibrate model inputs from historical data
- Calibrate \(\sigma\) from historical log-returns (rolling and/or full-window estimate).
- Optionally calibrate \(r\) from short-rate proxy/T-bill data or assume a constant risk-free input for each test date.
- Re-run pricing with calibrated parameters and compare behavior to baseline synthetic settings.

### 10.2 Compare simulations to market option quotes
- Collect market option quotes (same underlying/date, known strike/maturity).
- Compute pricing error metrics between MC prices and observed quotes.
- Analyze where model mismatch appears (vol smile, jumps, path features not captured by GBM).

### 10.3 Report framing
- Keep Phase 1 as the core methods/HPC contribution.
- Present Phase 2 as a realism stress test:
  - what improves with calibration
  - what remains limited by model assumptions
