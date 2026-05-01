# Phase 1 Results (Methods-First Baseline)

This document summarizes measured results from the PETSc-based parallel Monte Carlo implementation.

## 1) Pricing Smoke Tests (2 ranks, 20,000 paths, 100 steps)

- European call:
  - `price = 10.57326187`
  - `std_error = 0.10513776`
  - `95% CI = [10.36719185, 10.77933189]`
- Asian call:
  - `price = 5.90038797`
  - `std_error = 0.05757424`
  - `95% CI = [5.78754246, 6.01323348]`
- Barrier call (up-and-out, `B=130`):
  - `price = 3.65435585`
  - `std_error = 0.04574728`
  - `95% CI = [3.56469119, 3.74402051]`

## 2) European Validation vs Black-Scholes

Run configuration: 4 ranks, 400,000 paths, 252 steps.

- Monte Carlo price: `10.42161941`
- Black-Scholes price: `10.45058357`
- Absolute error: `2.89641665e-02`
- Relative error: `2.77153580e-03`
- Standard error: `2.32266613e-02`
- 95% CI: `[10.37609515, 10.46714366]`

Interpretation: the analytical value lies inside the MC confidence interval, consistent with correct implementation.

## 3) Strong Scaling (fixed total paths = 400,000)

- `p=1`: `T=4.8319s`, `S=1.000`, `E=1.000`
- `p=2`: `T=3.1022s`, `S=1.558`, `E=0.779`
- `p=4`: `T=2.8031s`, `S=1.724`, `E=0.431`
- `p=8`: `T=2.4953s`, `S=1.936`, `E=0.242`

Figures:
- `visuals/strong_scaling_speedup.png`
- `visuals/strong_scaling_efficiency.png`

## 4) Weak Scaling (fixed paths/rank = 100,000)

- `p=1`: `T=1.3559s`
- `p=2`: `T=1.7754s`
- `p=4`: `T=2.9858s`
- `p=8`: `T=5.4020s`

Figures:
- `visuals/weak_scaling_speedup.png`
- `visuals/weak_scaling_efficiency.png`

Interpretation: runtime increases with processor count in this local environment, indicating communication/runtime overhead dominates at higher `p`.

## 5) Convergence Study (Error vs N)

Configuration: 4 ranks, `N = [2e4, 5e4, 1e5, 2e5, 4e5, 8e5]`, 5 independent repetitions per `N`.

Fitted log-log slopes:
- Mean absolute error slope: `-0.547067`
- RMSE slope: `-0.529293`

Theoretical Monte Carlo rate: `-0.5` (i.e., error `O(N^{-1/2})`).

Conclusion: observed slopes are close to theory and support the expected convergence behavior.

Figure:
- `visuals/convergence.png`

## 6) Generated Artifacts

- Numerical CSVs:
  - `output/strong_scaling.csv`
  - `output/weak_scaling.csv`
  - `output/convergence.csv`
- ParaView CSV data:
  - `output/path_bundle.csv`
  - `output/barrier_plane.csv`
- Plots:
  - all PNGs in `visuals/`
