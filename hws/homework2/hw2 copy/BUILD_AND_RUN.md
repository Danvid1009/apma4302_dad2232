# Why "it doesn't work" and how to fix it

## The problem

- **Build error:** "PETSc was configured with Open MPI but now appears to be compiling using a non-Open MPI mpi.h"  
  → Your PETSc was built with **Open MPI**, but your current `mpicc` is **MPICH**. They must match.

- **Run error:** "MPI_ERR_COMM: invalid communicator"  
  → Happens when using **Open MPI** on macOS with this PETSc. Running with **MPICH** avoids it.

So the reliable path is: **use MPICH for everything** (configure PETSc with MPICH, build bvp with MPICH, run with MPICH's mpiexec).

---

## Fix: One-time PETSc setup with MPICH

1. **Use MPICH**
   ```bash
   brew unlink open-mpi 2>/dev/null
   brew link mpich --force
   ```

2. **Configure PETSc** (run from petsc dir; do NOT Ctrl+C — wait ~15–20 min)
   ```bash
   cd /Users/dan/Desktop/Columbia/HPC_4302/petsc
   export PETSC_DIR=$PWD
   export PETSC_ARCH=apma4302-base-opt
   ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
     --with-debugging=0 --with-shared-libraries=1 \
     --COPTFLAGS="-O3 -march=native" --CXXOPTFLAGS="-O3 -march=native" \
     --FOPTFLAGS="-O3 -march=native"
   ```
   Wait until it prints a summary and returns to the prompt.

3. **Build PETSc**
   ```bash
   make all
   ```
   (Same directory; same PETSC_DIR and PETSC_ARCH.)

---

## Build and run the BVP (every time)

```bash
cd "/Users/dan/Desktop/Columbia/HPC_4302/apma4302_dad2232/hws/homework2/hw2 copy"
export PETSC_DIR=/Users/dan/Desktop/Columbia/HPC_4302/petsc
export PETSC_ARCH=apma4302-base-opt
make bvp
mpiexec -np 1 ./bvp -options_file options_file
```

Use the same MPI (MPICH) for this: `which mpicc` should point to Homebrew's mpich, and run with that same environment.

---

## Quick check

- `which mpicc` → should be `/opt/homebrew/bin/mpicc` (MPICH when linked).
- If you see the "Open MPI but now ... non-Open MPI" compile error → PETSc was built with Open MPI; reconfigure PETSc with MPICH (steps above) and run `make all`, then build bvp again.
