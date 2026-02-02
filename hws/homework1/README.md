# Homework 1: exp(x) Taylor approximation (PETSc + MPI)

Parallel N-term Taylor approximation of exp(x). You need **expx.c** and a **makefile** to build and run.

## Build

- **Option A — use config:** The repo includes **config.mk** with `PETSC_DIR` and `PETSC_ARCH`. The makefile includes it so `make expx` works if your PETSc paths match (edit `config.mk` if needed).
- **Option B — set env locally:** Set PETSc in your shell instead of using config:
  ```bash
  export PETSC_DIR=/path/to/petsc
  export PETSC_ARCH=your-arch
  make expx
  ```
  The makefile uses `?=` so these env vars override values from config.

## Run

```bash
mpiexec -np 4 ./expx -x -10 -N 50
```

Or use **run.sh** (after `chmod +x run.sh`): it sets the env, runs `make expx`, then `mpiexec` with defaults.

## Files

- **expx.c** — source (required)
- **makefile** — build rules (required)
- **config.mk** — optional; sets `PETSC_DIR` and `PETSC_ARCH` for make
- **run.sh** — optional; build + run in one command
