#!/bin/bash
# Run this after ~20 min to see if PETSc configure finished.
PETSC_DIR="${PETSC_DIR:-/Users/dan/Desktop/Columbia/HPC_4302/petsc}"
LOG="$PETSC_DIR/configure.log"

if [[ ! -f "$LOG" ]]; then
  echo "No configure.log found at $LOG"
  exit 1
fi

if grep -q "Configure stage complete\|Configure complete" "$LOG" 2>/dev/null; then
  echo "Configure FINISHED successfully."
  echo "Run these next:"
  echo "  cd $PETSC_DIR && export PETSC_DIR=$PETSC_DIR PETSC_ARCH=apma4302-base-opt && make all"
  echo "  cd \"$(cd "$(dirname "$0")" && pwd)\" && export PETSC_DIR=$PETSC_DIR PETSC_ARCH=apma4302-base-opt && make bvp && mpiexec -np 1 ./bvp -options_file options_file"
elif grep -qi "Error\|error:" "$LOG" 2>/dev/null; then
  echo "Configure may have failed. Last 30 lines of log:"
  tail -30 "$LOG"
else
  echo "Configure still running or log incomplete. Last 10 lines:"
  tail -10 "$LOG"
  echo ""
  echo "Wait a bit longer and run this script again, or: tail -f $LOG"
fi
