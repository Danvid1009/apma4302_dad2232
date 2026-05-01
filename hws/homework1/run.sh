#!/usr/bin/env bash
# Run expx: set PETSc env, build, then mpiexec.
# Usage: ./run.sh [mpiexec args --] [expx args]
# Example: ./run.sh
# Example: ./run.sh -np 8 -- -x -10 -N 50

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PETSC_DIR="${PETSC_DIR:-/Users/dan/Desktop/Columbia/HPC_4302/petsc}"
export PETSC_ARCH="${PETSC_ARCH:-apma4302-base-debug}"

cd "$SCRIPT_DIR" || exit 1
make expx || exit 1

# Default: 4 processes, -x -10 -N 50. Or: ./run.sh -np 8 ./expx -x 1 -N 20
if [[ $# -eq 0 ]]; then
  mpiexec -np 4 ./expx -x -10 -N 50
else
  mpiexec "$@"
fi
