#!/bin/bash
# Q4: Solver performance — run each solver and show iteration count.
# From "for submission": source env_pkgs_opt, then ./run_q4.sh

set -e
OPTS="-options_file options_file -ksp_converged_reason"

echo "=== (a) Richardson + Jacobi, 1 proc ==="
mpiexec -np 1 ./bvp $OPTS -ksp_type richardson -pc_type jacobi 2>&1 | grep -E "Linear solve converged|DIVERGED|Relative error" || true

echo "=== (b) CG, no PC, 1 proc ==="
mpiexec -np 1 ./bvp $OPTS -ksp_type cg -pc_type none 2>&1 | grep -E "Linear solve converged|Relative error"

echo "=== (c) CG, no PC, c=0, 1 proc ==="
mpiexec -np 1 ./bvp $OPTS -ksp_type cg -pc_type none -bvp_c 0 2>&1 | grep -E "Linear solve converged|Relative error"

echo "=== (d) CG + ICC, 1 proc ==="
mpiexec -np 1 ./bvp $OPTS -ksp_type cg -pc_type icc 2>&1 | grep -E "Linear solve converged|Relative error"

echo "=== (e) CG + bjacobi(icc), 4 procs ==="
mpiexec -np 4 ./bvp $OPTS -ksp_type cg -pc_type bjacobi -sub_pc_type icc 2>&1 | grep -E "Linear solve converged|Relative error"

echo "=== (f) MUMPS, 1 proc (direct, no iterations) ==="
mpiexec -np 1 ./bvp -options_file options_file -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps 2>&1 | grep "Relative error"

echo "=== (f) MUMPS, 4 procs (direct) ==="
mpiexec -np 4 ./bvp -options_file options_file -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps 2>&1 | grep "Relative error"

echo "Done."
