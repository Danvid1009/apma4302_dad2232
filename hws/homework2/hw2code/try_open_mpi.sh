#!/bin/bash
# Try building and running with Open MPI (your current PETSc is built with Open MPI).
# If the run fails with MPI_ERR_COMM, use the MPICH path instead (see BUILD_AND_RUN.md).
set -e
HW2="$(cd "$(dirname "$0")" && pwd)"
PETSC_DIR="${PETSC_DIR:-/Users/dan/Desktop/Columbia/HPC_4302/petsc}"
PETSC_ARCH="${PETSC_ARCH:-apma4302-base-opt}"

echo "=== Switching to Open MPI ==="
brew unlink mpich 2>/dev/null || true
brew link open-mpi --force 2>/dev/null || true

echo "=== Building bvp ==="
cd "$HW2"
export PETSC_DIR PETSC_ARCH
make clean-bvp 2>/dev/null || true
make bvp

echo "=== Running bvp (with Open MPI workaround flags; may still fail on Mac) ==="
export OMPI_MCA_btl=self,sm,tcp
export OMPI_MCA_btl_tcp_if_include=lo0
mpiexec -np 1 ./bvp -options_file options_file && echo "SUCCESS" || echo "Run failed (try MPICH path in BUILD_AND_RUN.md)"
