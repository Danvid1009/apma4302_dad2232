# Homework 4 — coupled multiphysics (Firedrake + PETSc)

Outline of where everything lives and how the pieces connect to the assignment.

---

## Where everything is

### Top level

| Item | Purpose |
|------|---------|
| **`RUN.md`** | Full setup: PETSc/Firedrake, cluster Apptainer env vars, **recommended run order** (Q2 → Q3 → Q4 → bonus). |
| **`.gitignore`** | Ignores regenerated `output/`, built binary `c/biharm`, LaTeX aux files, `__pycache__/`, local `.venv/`. |
| **`README.md`** | This file — map of the tree. |

### `doc/` — write-up (LaTeX + PDF + figures)

| Item | Purpose |
|------|---------|
| **`hw4.tex`** | Original assignment handout (problem statement). |
| **`solutions.tex`** | Your solutions document (references code paths below). |
| **`hw4.pdf`** | Last committed PDF build of the write-up (rebuild from `solutions.tex` if you change content). |
| **`doc/README.md`** | How to run `pdflatex` on `solutions.tex`. |
| **`doc/SCRIPTS.md`** | Table mapping each problem / figure mention in the PDF to **repo files** and how to run them. |
| **`doc/figures/`** | Handout scans (`hw4_*.png`) + optional exported PDFs for placeholders (see `figures/README.md`). |

### `c/` — PETSc (Problem 2a)

| Item | Purpose |
|------|---------|
| **`biharm.c`**, **`poissonfunctions.c`**, **`poissonfunctions.h`** | DMDA finite-difference biharmonic + manufactured solution. |
| **`makefile`** | Build `biharm` (needs `PETSC_DIR` / `PETSC_ARCH`). |
| **`options_file_*`** | Three solver presets (direct, split+direct, split+MG). |
| **`*.log`** (if present) | Captured runs from local experiments. |

### `python/` — Firedrake (Problems 2b–4 + extra credit)

| Item | Purpose |
|------|---------|
| **`biharm.py`** | FE biharmonic; same three solver presets as Q2a. |
| **`biharm_temperature_rhs.py`** | Q3: temperature-driven RHS; VTK for ParaView. |
| **`biharm_temp.py`** | Q3 (write-up variant): same physics as above, `T`/`f` via `interpolate` + `biharm_temp.pvd`. |
| **`plot_q4_nusselt.py`** | Q4: build Nu vs time / Nu vs mesh PNGs from `nu_history.csv` (or `--demo` for layout drafts). |
| **`convection.py`** | Q4: full coupled convection DAE (`firedrake-ts`), Nusselt CSV. |
| **`heat.py`** | Reference heat / TS example used to build `convection.py`. |
| **`convection_cn.py`** | Extra (3): Crank–Nicolson + SNES per step. |
| **`bonus_convergence_sweep.py`** | Extra (1): mesh × Ra sweep. |
| Other `*.py` | Supporting experiments (e.g. convection variants) as you add them. |

### `scripts/` — one-command runners

| Script | Runs |
|--------|------|
| **`run_q2.sh`** | C `biharm` (three options) + Firedrake `biharm.py` (three presets); logs under `output/logs/`. |
| **`run_q3.sh`** | `biharm_temperature_rhs.py` → VTK under `output/q3/`. |
| **`run_q4.sh`** | `convection.py` with env vars `RA`, `N`, `TMAX`, `DT`, `HW4_OUT`, etc. |
| **`plot_q4_figures.sh`** | `plot_q4_nusselt.py` → `doc/figures/hw4_q4_*.png` (default `--demo`; pass script args for real CSVs). |
| **`run_bonus.sh`** | `sweep` → `bonus_convergence_sweep.py`; `cn` → `convection_cn.py`. |
| **`inc_firedrake_apptainer.sh`** | Shared helper for running `python3` inside a Firedrake Apptainer image on clusters. |

### `output/` — generated artifacts (not required in clone)

Created by the scripts: logs, VTK/PVD, `nu_history.csv`, bonus CSVs, etc. Regenerate after `git pull` using **`RUN.md`**.

---

## Quick links

- **Build PDF:** [`doc/README.md`](doc/README.md)  
- **Run all numerics:** [`RUN.md`](RUN.md)  
- **PDF ↔ code map + gap-fill checklist:** [`doc/SCRIPTS.md`](doc/SCRIPTS.md) (**“Checklist — fill gaps”** at bottom)

---

## NOT YET COMPLETE

The following are still open relative to a “finished” homework submission (see `doc/solutions.tex` for `\missing{...}` and `\placeholdfig` blocks):

- **Problem 2 (Firedrake):** §2(b) filled from `doc/data/fe_biharm/`; §2(c) from `doc/data/q2c_split_mg/` (optional bar chart remains). Refresh from cluster logs if you prefer your own numbers.  
- **Problem 3:** Export ParaView figure to PDF and replace the placeholder (`figures/hw4_q3_T_psi_velocity.pdf` or edit `solutions.tex`).  
- **Problem 4:** Steady-state / Nusselt discussion for (a)–(c); plots for `Nu(t)`, `Nu` vs mesh vs Blankenbach; mesh sweep numbers.  
- **Extra credit:** Populate bonus sweep table; optional heatmaps / BDF vs CN overlay figures; PETSc C+TS port (2) remains unimplemented unless you add it.  
- **PDF rebuild:** After edits, re-run `pdflatex` in `doc/` and refresh `hw4.pdf` if you want the committed PDF to match the latest `solutions.tex`.
