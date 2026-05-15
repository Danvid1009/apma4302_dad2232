# Restart Implementation Outline

## Objective
Rebuild the project around a PETSc-first workflow where each numerical experiment produces ParaView-ready outputs before additional analysis.

## Phase 1 - PETSc Alignment (must pass first)
1. Export environment:
   - `PETSC_ARCH=apma4302-pkgs-opt`
   - PETSc `bin` in `PATH`
   - PETSc `lib` in `PYTHONPATH`
2. Validate with:
   - `python scripts/check_petsc_env.py`
3. Execute smoke solve + export:
   - `python scripts/smoke_paraview_export.py`
4. Gate:
   - `output/paraview/vtp/poisson.pvd` exists and opens in ParaView.

## Phase 2 - Core Numerical Problems
1. Poisson (elliptic) baseline:
   - validate residual and solution smoothness.
2. Heat equation (parabolic):
   - explicit or implicit in time with PETSc linear solves.
3. Advection-diffusion:
   - add transport term and visualize fronts/plumes.

## Phase 3 - Quantitative + Performance
1. Convergence tests for spatial/temporal refinement.
2. Strong scaling and efficiency tables.
3. Solver configuration sweep:
   - KSP/PC comparisons.

## Phase 4 - Report and Visual Packaging
1. Export final PNG figures under `visuals/`.
2. Fill LaTeX template from generated metrics.
3. Keep ParaView instructions synchronized with real outputs.
