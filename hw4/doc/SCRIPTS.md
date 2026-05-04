# Scripts and code referenced by `solutions.tex`

Paths are relative to the **repository root** (`apma4302_dad2232/`), unless noted.

| Problem / topic | File in PDF | Repo path | How to run (summary) |
|-----------------|-------------|-----------|----------------------|
| Q2a — PETSc C biharm | `\texttt{c/biharm.c}` | `hw4/c/biharm.c` | `cd hw4/c && make && ./biharm -options_file …` or `hw4/scripts/run_q2.sh` |
| Q2a — options | (implicit) | `hw4/c/options_file_direct`, `options_file_split_direct`, `options_file_split_mg` | Passed to `./biharm` |
| Q2b — Firedrake biharm | `\texttt{biharm.py}` | `hw4/python/biharm.py` | `hw4/scripts/run_q2.sh` or `python3 biharm.py --preset …` |
| Q3 — temperature RHS | `\texttt{python/biharm\_temperature\_rhs.py}` | `hw4/python/biharm_temperature_rhs.py` | `hw4/scripts/run_q3.sh` |
| Q3 — VTK output | `\texttt{hw4/output/q3/biharm\_T\_rhs.pvd}` | created under `hw4/output/q3/` | After `run_q3.sh` |
| Q4 — full convection DAE | `\texttt{python/convection.py}` | `hw4/python/convection.py` | `hw4/scripts/run_q4.sh` (env `RA`, `N`, `TMAX`, …) |
| Q4 — Nusselt CSV | `\texttt{nu\_history.csv}` | under each run’s `--output-dir` / `HW4_OUT` | Written by `convection.py` |
| Extra (1) — sweep | `\texttt{python/bonus\_convergence\_sweep.py}` | `hw4/python/bonus_convergence_sweep.py` | `hw4/scripts/run_bonus.sh sweep …` |
| Extra (1) — wrapper | `\texttt{hw4/scripts/run\_bonus.sh}` | `hw4/scripts/run_bonus.sh` | `chmod +x hw4/scripts/run_bonus.sh` |
| Extra (3) — CN / SNES | `\texttt{python/convection\_cn.py}` | `hw4/python/convection_cn.py` | `hw4/scripts/run_bonus.sh cn …` |

Heat-only DAE reference code (handout): `hw4/python/heat.py` (used to build `convection.py`; not always cited by name in `solutions.tex`).

For full commands and cluster Apptainer setup, use **`../RUN.md`**.
