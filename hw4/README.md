# Homework 4 — coupled multiphysics (Firedrake + PETSc)

Layout for submission and `git clone` reproduction:

| Path | Contents |
|------|----------|
| **`doc/`** | Assignment handout (`hw4.tex`), solutions (`solutions.tex`), built **`hw4.pdf`**, scanned handout PNGs under **`doc/figures/`**, LaTeX build notes. |
| **`c/`** | PETSc finite-difference biharmonic driver (`biharm.c`, `poissonfunctions.*`), `makefile`, `options_file_*`. |
| **`python/`** | Firedrake drivers: `biharm.py`, `biharm_temperature_rhs.py`, `convection.py`, `convection_cn.py`, `heat.py`, `bonus_convergence_sweep.py`, etc. |
| **`scripts/`** | `run_q2.sh`, `run_q3.sh`, `run_q4.sh`, `run_bonus.sh`, `inc_firedrake_apptainer.sh`. |
| **`RUN.md`** | Environment setup, cluster Apptainer notes, and recommended run order. |
| **`output/`** | Logs, VTK, CSV from runs (see `.gitignore`; regenerate with the scripts). |

**Compile the write-up:** see [`doc/README.md`](doc/README.md).

**Run numerics:** start at [`RUN.md`](RUN.md).

**Which scripts match the PDF:** see [`doc/SCRIPTS.md`](doc/SCRIPTS.md).
