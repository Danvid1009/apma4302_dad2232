# APMA 4302 — Homework 3 (submission folder)

**Student:** Daniel (UNI: dad2232)

## Main write-up (PDF)

- **`Solutions_Daniel_HW3_4302.pdf`** — compiled solutions (Problems 1–8: theory, implementation summary, experiments, figures).

Source LaTeX: `latex/hw3.tex` (compile from `latex/` with `pdflatex hw3.tex`; figures are under `solutions/figures/`).

## Code (not duplicated here)

The nonlinear solver **`reaction2d`** is implemented in the **p4pdes** tree (not inside this folder):

- `p4pdes/c/ch4/reaction2d.c`
- Build: `cd p4pdes/c/ch4 && make reaction2d`

Course template / reference: `poisson2d/poisson2d.c` (linear Poisson + VTK style).

## PETSc option files (assignment)

| File | Purpose |
|------|--------|
| `solutions/options_file_gamma0` | Part 4: \(\gamma=0\), \(65\times 65\), one Newton, tight residual |
| `solutions/options_file_fd` | Part 5: same setup + finite-difference Jacobian (`-snes_fd_color`) |

Run from `p4pdes/c/ch4`, e.g.  
`mpirun -n 2 ./reaction2d -options_file /path/to/options_file_gamma0`

## Figures & logs

- `solutions/figures/paraview_u.png`, `paraview_uexact.png` — ParaView (Problem 7)
- `solutions/figures/snes_residual_history.png` — SNES residual history (Problem 6)
- `solutions/logs/scaling_hw3_runs.txt` — sample scaling numbers (Problem 8)
- `solutions/plot_snes_residual.py` — regenerate residual plot if needed
- `solutions/run_scaling_study.sh` — optional sweep over `-da_refine` and MPI ranks

## Extra notes

See `solutions/README.md` for the original assignment checklist / layout.

Environment variables for PETSc builds: `PETSC_DIR`, `PETSC_ARCH` (see course instructions).
