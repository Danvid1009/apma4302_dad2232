#!/usr/bin/env bash
# Part 8: sweep -da_refine 2..6 and MPI ranks 1,2,4 (run from p4pdes/c/ch4 with PETSC_* set)
set -euo pipefail
R2D="${R2D:-./reaction2d}"
for r in 2 3 4 5 6; do
  for np in 1 2 4; do
    echo "=== da_refine=$r  mpi=$np ==="
    /usr/bin/time -p mpirun -n "$np" "$R2D" -rct_gamma 100 -rct_p 3 -da_refine "$r" \
      -snes_max_it 20 -snes_atol 1e-12 -snes_rtol 0 \
      -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps \
      -options_left 0 2>&1 | grep -E '^(on |SNES iterations|real )'
  done
done
