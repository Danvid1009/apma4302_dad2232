#!/bin/bash
# Run PETSc configure with MPICH in the background. You can close the terminal;
# configure keeps running. Come back in 20 min and run the "after configure" steps.
PETSC_DIR="${PETSC_DIR:-/Users/dan/Desktop/Columbia/HPC_4302/petsc}"
PETSC_ARCH="${PETSC_ARCH:-apma4302-base-opt}"
LOG="$PETSC_DIR/configure.log"

echo "=== Ensure MPICH is active ==="
brew unlink open-mpi 2>/dev/null || true
brew link mpich --force 2>/dev/null || true

echo "=== Starting PETSc configure in background (log: $LOG) ==="
cd "$PETSC_DIR"
export PETSC_DIR PETSC_ARCH
nohup ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
  --with-debugging=0 --with-shared-libraries=1 \
  --COPTFLAGS="-O3 -march=native" --CXXOPTFLAGS="-O3 -march=native" \
  --FOPTFLAGS="-O3 -march=native" > "$LOG" 2>&1 &
CONF_PID=$!
echo "Configure PID: $CONF_PID"
echo ""
echo "--- What to do next ---"
echo "1. Watch progress:  tail -f $LOG"
echo "2. In 15–20 minutes, check if it's done:  grep -i 'configure stage complete\\|Configure complete\\|Error' $LOG | tail -5"
echo "3. When you see 'Configure stage complete' (and no fatal Error), run:"
echo "   cd $PETSC_DIR && export PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH && make all"
echo "4. Then build and run bvp:"
echo "   cd \"$(dirname "$0")\" && export PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH && make bvp && mpiexec -np 1 ./bvp -options_file options_file"
echo ""
