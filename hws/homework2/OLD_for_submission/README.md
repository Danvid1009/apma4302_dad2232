# HW2 submission

- **bvp.c** — BVP solver (from tri.c, p4pdes Ch. 2), run with `mpiexec -np P ./bvp -options_file options_file`
- **Makefile** — build with `make bvp`
- **options_file** — PETSc options (BVP parameters and KSP options)
- **plot_bvp.py** — plot solution from `bvp_solution.h5` (3b) and convergence study (3c)
- **hw2_solutions.tex** — written solutions
- **run_q4.sh** — run all Q4 solver configurations and print iteration counts

## Q4: Solver performance

Run all solver configs and print iteration counts (table and explanations in `hw2_solutions.tex`):

```bash
source env_pkgs_opt
./run_q4.sh
```

For (f) the script uses `-pc_factor_mat_solver_type mumps` (PETSc 3.14+).

## Q3(b) completed — Plots

3(b) requirements: PETSc options / `-options_file`, parallel assembly, `MatZeroRowsColumns`, relative error printed, HDF5 output (`u`, `f`, `uexact`), and Python script to plot solution, true solution, and error. All satisfied.

**Generated plots (in this directory):**
- **bvp_solution.pdf** — Single run: numerical vs exact solution and pointwise error (from `bvp_solution.h5`).
- **bvp_convergence.pdf** — Q3(c): error vs h (log-log) for γ=0, k=1,5,10, m=40,…,1280; order of convergence ≈ 2.

To regenerate: run BVP (step 3 below), then `source env_plot_bvp` and `python3 plot_bvp.py` for the single-solution plot; use `BVP_NPROC=8 python3 plot_bvp.py convergence ./bvp` for the convergence plot.

---

## Q3(b) walkthrough: single run and plot

Do these steps **in order**. Use the same shell (or re-source the envs when you open a new one).

### 1. PETSc with HDF5

Your PETSc must be built with **HDF5** so the code can write `bvp_solution.h5`. If you built with `--download-hdf5`, you’re set. To check:

```bash
# From PETSc root:
./config/build.log  # or grep -i hdf5 config/BuildSystem/config/packages/hdf5.py
# Or just try building bvp; the code uses PetscViewerHDF5Open.
```

If HDF5 wasn’t used, reconfigure PETSc with `--download-hdf5`, then rebuild (`make all`).

### 2. Set environment and build

From the **`for submission`** directory:

```bash
cd "/path/to/for submission"
source env_pkgs_opt
make bvp
```

You should get an executable `./bvp` with no errors.

### 3. Run the BVP

```bash
mpiexec -np 1 ./bvp -options_file options_file
```

- Use `-np 1` for a single process, or `-np 2`, etc., for parallel.
- The program reads options from `options_file` (e.g. `-bvp_m 201`, `-bvp_k 5`, `-ksp_rtol`, etc.).
- On success you get:
  - One line: `Relative error (2-norm): ...`
  - A file **`bvp_solution.h5`** in the same directory (contains vectors `u`, `uexact`, `f`).

**If you see MPI errors on your Mac** (e.g. `MPI_ERR_COMM: invalid communicator` with Open MPI): the BVP and options are correct; the failure is from the MPI stack. Run the same command on a **cluster** (or with PETSc built with MPICH and that `mpiexec`); it should work there.

### 4. Plot with Python (petsc4py)

**Requires:** `bvp_solution.h5` from step 3, and **petsc4py** (from your PETSc build).

```bash
source env_plot_bvp
python3 plot_bvp.py
```

- `env_plot_bvp` sets `PYTHONPATH` so Python finds `petsc4py` in `$PETSC_DIR/$PETSC_ARCH/lib`.
- The script reads `bvp_solution.h5` with the PETSc HDF5 viewer and plots numerical solution, exact solution, and error.
- If you set `BVP_PLOT_OUT`, the figure is saved to that file; otherwise it uses the default (e.g. `bvp_solution.pdf` if that’s set in the script).

**If `petsc4py` fails to import** (e.g. loader error on some setups): run step 4 on the **same machine/environment where step 3 worked** (e.g. cluster), or use a Python that matches the PETSc/petsc4py build (e.g. Python 3.9 as in the env comment).

### Summary

| Step | Command |
|------|--------|
| 1 | PETSc built with `--download-hdf5` |
| 2 | `source env_pkgs_opt` → `make bvp` |
| 3 | `mpiexec -np 1 ./bvp -options_file options_file` → creates `bvp_solution.h5` |
| 4 | `source env_plot_bvp` → `python3 plot_bvp.py` → plot (and optional PDF) |

---

**Build and run** (use pkgs-opt; make and make check passed in PETSc):
```bash
source env_pkgs_opt   # or: export PETSC_DIR=... PETSC_ARCH=apma4302-pkgs-opt
make bvp
mpiexec -np P ./bvp -options_file options_file
```
