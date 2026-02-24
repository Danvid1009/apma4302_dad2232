#!/bin/bash
# One script: build BVP and try to run it. If run fails (Open MPI bug), print exact fix.
set -e
HW2="$(cd "$(dirname "$0")" && pwd)"
PETSC_DIR="${PETSC_DIR:-/Users/dan/Desktop/Columbia/HPC_4302/petsc}"
PETSC_ARCH="${PETSC_ARCH:-apma4302-base-opt}"

echo "=== Using Open MPI (so existing PETSc build matches) ==="
brew unlink mpich 2>/dev/null || true
brew link open-mpi --force 2>/dev/null || true

echo "=== Building bvp ==="
cd "$HW2"
export PETSC_DIR PETSC_ARCH
make clean-bvp 2>/dev/null || true
make bvp
echo "Build OK."

echo "=== Trying to run (workarounds for Mac Open MPI) ==="
for btl in "self,sm" "self,sm,tcp" "self" "tcp,self"; do
  echo "Try: OMPI_MCA_btl=$btl"
  export OMPI_MCA_btl="$btl"
  if mpiexec -np 1 ./bvp -options_file options_file 2>/dev/null; then
    echo "SUCCESS with OMPI_MCA_btl=$btl"
    exit 0
  fi
done

echo ""
echo "--- Run still fails on this Mac (Open MPI bug). ---"
echo "To get a working run you must build PETSc with MPICH and use it:"
echo ""
echo "1. In a terminal (leave it open 15–20 min, do not Ctrl+C):"
echo "   brew link mpich --force"
echo "   cd $PETSC_DIR"
echo "   export PETSC_DIR=\$PWD PETSC_ARCH=apma4302-base-opt"
echo "   ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-debugging=0 --with-shared-libraries=1 --COPTFLAGS=\"-O3 -march=native\" --CXXOPTFLAGS=\"-O3 -march=native\" --FOPTFLAGS=\"-O3 -march=native\""
echo ""
echo "2. When configure finishes, run:"
echo "   make all"
echo ""
echo "3. Then:"
echo "   cd \"$HW2\""
echo "   export PETSC_DIR=$PETSC_DIR PETSC_ARCH=apma4302-base-opt"
echo "   make bvp"
echo "   mpiexec -np 1 ./bvp -options_file options_file"
echo ""
exit 1
