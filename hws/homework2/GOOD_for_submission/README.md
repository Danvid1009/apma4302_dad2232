# APMA 4302 — Homework 2 Submission

This folder contains the complete submission for Homework 2: written solutions (PDF), BVP code and executable, solution data, and the required plots.

---

## Contents

### Written solutions (all problems)

| File | Description |
|------|-------------|
| **HW2_4302_LaTeX.pdf** | **All written solutions** for HW2. Contains: **Problem 1** (condition number and error bounds, with proofs); **Problem 2** (discrete Laplacian — eigenvectors, eigenvalues, κ(A) ∼ O(m²)); **Problem 3** (BVP — manufactured solution and f(x), description of code and runs, convergence); **Problem 4** (solver performance — iteration counts for Richardson+Jacobi, CG no PC, CG c=0, CG+ICC, bjacobi+ICC, MUMPS, with explanations). |

### Code and program (Problem 3)

| File | Description |
|------|-------------|
| **bvp.c** | Source code for the BVP solver. Solves −u″ + γu = f on [0,1] with Dirichlet boundary conditions. Implements **3(a)** manufactured solution (default: u(x) = sin(kπx) + c/8, f(x) = (k²π²+γ) sin(kπx) + γc/8). Implements **3(b)**: PETSc options (`-bvp_m`, `-bvp_gamma`, `-bvp_k`, `-bvp_c`, and `-options_file`), parallel assembly, `MatZeroRowsColumns` for BCs (symmetry preserved), relative error printed, HDF5 output of u, f, uexact. Modified from tri.c (p4pdes, Ch. 2). |
| **bvp** | Compiled executable. Built with PETSc (HDF5 support). Intended run: `mpiexec -np P ./bvp -options_file options_file` (requires an options file and PETSc environment). |
| **options_file** | Default PETSc/BVP options (m=201, γ=0, k=5, c=3, KSP tolerances). Used by **bvp** and **run_q4.sh**. |
| **run_q4.sh** | Script that runs all Problem 4 solver configs (a)–(f) and prints convergence/iteration info. Usage: set PETSc env (e.g. `source env_pkgs_opt`), then `./run_q4.sh`. Requires **options_file** and **bvp** in the same directory. For (f) uses `-pc_factor_mat_solver_type mumps`. |
| **bvp_solution.h5** | HDF5 output from a successful BVP run. Contains vectors **u** (numerical solution), **uexact** (exact solution), and **f** (right-hand side). Used to produce the single-solution plot. |

### Plots (Problems 3(b) and 3(c))

| File | Description |
|------|-------------|
| **bvp_solution.pdf** | **3(b) Single-run plot.** Numerical solution vs exact solution u(x) = sin(kπx)+c/8, and pointwise error (right axis). Generated from `bvp_solution.h5` for one run (e.g. m=201, γ=0, k=5, c=3). |
| **bvp_convergence.pdf** | **3(c) Convergence plot.** Relative error (2-norm) vs mesh spacing h on a log-log scale for γ=0, k=1,5,10, m=40,80,160,…,1280. Order of convergence ≈ 2, consistent with the second-order finite difference scheme. |

---

## How the pieces fit the assignment

- **Problems 1 and 2:** Written work only → see **HW2_4302_LaTeX.pdf**.
- **Problem 3(a):** Manufactured solution and f(x) → derived in **HW2_4302_LaTeX.pdf**; implemented in **bvp.c** (see comments and `u_exact_fn`, `f_rhs`).
- **Problem 3(b):** BVP code and run → **bvp.c** and **bvp**; options and `-options_file`, parallel assembly, `MatZeroRowsColumns`, relative error, HDF5 output as required. **bvp_solution.pdf** is the plot of the solution, true solution, and error from **bvp_solution.h5**.
- **Problem 3(c):** Convergence study → **bvp_convergence.pdf** (error vs h and order); discussion in **HW2_4302_LaTeX.pdf**.
- **Problem 4:** Solver performance → iteration counts and explanations in **HW2_4302_LaTeX.pdf**. **run_q4.sh** reproduces the runs (Richardson+Jacobi, CG no PC, CG c=0, CG+ICC, bjacobi+ICC, MUMPS 1 and 4 procs).

---

## Running the executable (optional)

The **bvp** binary was built with PETSc (with HDF5). To run it elsewhere you need the same PETSc install and an options file. Example:

```bash
# Set PETSc environment (PETSC_DIR, PETSC_ARCH), then:
mpiexec -np 1 ./bvp -options_file options_file
```

An options file should define at least `-bvp_m`, `-bvp_gamma`, `-bvp_k`, `-bvp_c` (and optionally KSP options). **options_file** is included. The program prints the relative error and writes **bvp_solution.h5** when HDF5 is available.

To reproduce Q4 iteration counts: set PETSc environment, then run `./run_q4.sh`.

---

## Summary

| Item | Purpose |
|------|---------|
| HW2_4302_LaTeX.pdf | All written solutions (Problems 1–4) |
| bvp.c | BVP source (3(a) manufactured solution, 3(b) full spec) |
| bvp | BVP executable |
| options_file | Default BVP/KSP options (for bvp and run_q4.sh) |
| run_q4.sh | Run all Q4 solver configs and print iteration info |
| bvp_solution.h5 | HDF5 output from one BVP run (u, uexact, f) |
| bvp_solution.pdf | 3(b) plot: numerical vs exact and error |
| bvp_convergence.pdf | 3(c) plot: error vs h, order ≈ 2 |
