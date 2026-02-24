#!/bin/bash
# One-time setup: configure PETSc with MPICH, build PETSc and bvp, then run bvp.
# Run from this directory. Requires: brew link mpich (and open-mpi unlinked).
set -e
PETSC_DIR="${PETSC_DIR:-/Users/dan/Desktop/Columbia/HPC_4302/petsc}"
HW2_DIR="$(cd "$(dirname "$0")" && pwd)"
export PETSC_DIR
export PETSC_ARCH=apma4302-base-debug

echo "=== Step 1: Configure PETSc with MPICH (this may take 5–10 minutes) ==="
cd "$PETSC_DIR"
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
  --with-debugging=1 --with-shared-libraries=1 || { echo "Configure failed."; exit 1; }

echo "=== Step 2: Build PETSc ==="
make all

echo "=== Step 3: Build bvp ==="
cd "$HW2_DIR"
make clean-bvp 2>/dev/null || true
make bvp

echo "=== Step 4: Run bvp ==="
mpiexec -np 1 ./bvp -options_file options_file

echo "=== Done. To plot: python plot_bvp.py ==="
