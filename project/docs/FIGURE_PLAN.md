# Figure Plan (PNG Targets and Usage)

All report figures should be written to `visuals/` with stable names.

## Monte Carlo + MPI “story” figures
- `visuals/hpc_speedup_vs_ideal.png` — measured speedup vs linear ideal + efficiency (`plot_hpc_speedup_vs_ideal.py`).
- `visuals/hpc_weak_scaling_walltime.png` — weak scaling wall time near-flat (`plot_hpc_weak_scaling_walltime.py`).
- `visuals/mpi_domain_paths_tube.png` (ParaView export) — path polylines colored by synthetic **`mpi_rank`** (`export_mpi_domain_paths_vtp.py` + Tube filter; see [`PROJECT_OUTLINE.md` Appendix A](../PROJECT_OUTLINE.md#appendix-a-paraview-and-vtk-workflows), subsection **5) → E)**).

## PETSc + ParaView Figures
- `visuals/mesh_domain_3d.png`  
  Domain overview and mesh layout.
- `visuals/poisson_solution_surface.png`  
  Scalar field for Poisson solution.
- `visuals/poisson_contour_slice.png`  
  Slice + contours for spatial structure.
- `visuals/heat_snapshots_t0_tmid_tend.png`  
  Three-time comparison for heat equation.
- `visuals/advection_diffusion_volume.png`  
  Volume rendering / slice blend showing transport front.

## Solver / Numerical Figures
- `visuals/poisson_residual_history.png`  
  Iterative solver residual history.
- `visuals/heat_time_error_loglog.png`  
  Time-step error slope.
- `visuals/spatial_convergence_loglog.png`  
  Spatial refinement convergence.
- `visuals/strong_scaling_speedup.png`  
  Speedup vs ranks.
- `visuals/strong_scaling_efficiency.png`  
  Efficiency vs ranks.
- `visuals/solver_compare_runtime.png`  
  Runtime comparison across KSP/PC choices.

## Report Mapping (where each goes)
- Methods / setup section: `mesh_domain_3d.png`
- Results (Poisson): `poisson_solution_surface.png`, `poisson_contour_slice.png`, `poisson_residual_history.png`
- Results (Heat): `heat_snapshots_t0_tmid_tend.png`, `heat_time_error_loglog.png`
- Results (Advection-diffusion): `advection_diffusion_volume.png`
- Performance: `strong_scaling_speedup.png`, `strong_scaling_efficiency.png`, `solver_compare_runtime.png`
