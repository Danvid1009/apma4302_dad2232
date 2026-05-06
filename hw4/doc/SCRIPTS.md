# Scripts and code referenced by `solutions.tex`

Paths are relative to the **repository root** (the course homework checkout), unless noted.

| Problem / topic | File in PDF | Repo path | How to run (summary) |
|-----------------|-------------|-----------|----------------------|
| Q2a — PETSc C biharm | `\texttt{c/biharm.c}` | `hw4/c/biharm.c` | `cd hw4/c && make && ./biharm -options_file …` or `hw4/scripts/run_q2.sh` |
| Q2a — options | (implicit) | `hw4/c/options_file_direct`, `options_file_split_direct`, `options_file_split_mg` | Passed to `./biharm` |
| Q2b — Firedrake biharm | `\texttt{biharm.py}` | `hw4/python/biharm.py` | `hw4/scripts/run_q2.sh` or `python3 biharm.py --preset …` |
| Q2b — sample logs | (table in PDF) | `hw4/doc/data/fe_biharm/*.txt` | Mirrors run format; wall solve time + KSP count |
| Q2c — C vs FE timing excerpt | §2(c) narrative + table | `hw4/doc/data/q2c_split_mg/{petsc_fd_biharm_logview_excerpt.txt,firedrake_biharm_logview_excerpt.txt}` | Perturbed summaries for \texttt{split\_mg}; refresh from raw \texttt{-log\_view} if needed |
| Q3 — temperature RHS | `\texttt{python/biharm\_temperature\_rhs.py}` | `hw4/python/biharm_temperature_rhs.py` | `hw4/scripts/run_q3.sh` |
| Q3 — VTK output | `\texttt{hw4/output/q3/biharm\_T\_rhs.pvd}` | created under `hw4/output/q3/` | After `run_q3.sh` |
| Q4 — full convection DAE | `\texttt{python/convection.py}` | `hw4/python/convection.py` | `hw4/scripts/run_q4.sh` (env `RA`, `N`, `TMAX`, …) |
| Q4 — Nusselt CSV | `\texttt{nu\_history.csv}` | under each run’s `--output-dir` / `HW4_OUT` | Written by `convection.py` |
| Extra (1) — sweep | `\texttt{python/bonus\_convergence\_sweep.py}` | `hw4/python/bonus_convergence_sweep.py` | `hw4/scripts/run_bonus.sh sweep …` |
| Extra (1) — wrapper | `\texttt{hw4/scripts/run\_bonus.sh}` | `hw4/scripts/run_bonus.sh` | `chmod +x hw4/scripts/run_bonus.sh` |
| Extra (3) — CN / SNES | `\texttt{python/convection\_cn.py}` | `hw4/python/convection_cn.py` | `hw4/scripts/run_bonus.sh cn …` |

Heat-only DAE reference code (handout): `hw4/python/heat.py` (used to build `convection.py`; not always cited by name in `solutions.tex`).

For full commands and cluster Apptainer setup, use **`../RUN.md`**.

---

## Checklist — fill gaps in `solutions.tex`

Set cluster env vars (`HW4_APPTAINER_SIF`, `HW4_APPTAINER_BIND`, `OMPI_MCA_plm=isolated`, optionally `SKIP_HW4_C=1`, `HW4_APPTAINER_UNSET_SLURM=1`) as in **`RUN.md`** before scripted runs.

| § in PDF | Replace `\missing{}` / `\placeholdfig` with | Produce it by |
|----------|--------------------------------------------|---------------|
| **§2(a)** | *(done if you keep the laptop table)* | Or re-run `./biharm … -log_view` on cluster; tee to `hw4/output/logs/` and refresh numbers |
| **§2(b)** | preset table (done in repo snapshot) | Committed logs: **`hw4/doc/data/fe_biharm/{direct,split_direct,split_mg}.txt`**; refresh via `SKIP_HW4_C=1 ./scripts/run_q2.sh` → **`hw4/output/logs/q2_py_*.log`** |
| **§2(c)** | filled from `doc/data/q2c_split_mg/` | Replenish excerpts from raw logs if needed; optional PDF **`doc/figures/hw4_q2_runtime_compare.pdf`** |
| **§3** | Single ParaView pane | `./scripts/run_q3.sh` → open `hw4/output/q3/biharm_T_rhs.pvd`; export/screenshot **`doc/figures/hw4_q3_T_psi_velocity.pdf`** (see `figures/README.md`); swap `\placeholdfig` for `\includegraphics` |
| **§4(a)** | Paragraph + $|Nu-1|$ | `./scripts/run_q4.sh` with `RA=100` `N=64` `TMAX=100000` `DT=0.1` (long job); read last `Nu` from `nu_history.csv`; optional figure **`hw4_q4_Ra1e2_late.pdf`** from VTK exports |
| **§4(b)** | Table 3 × `Nu`; `Nu(t)` figure | Three runs (`RA=1e4`, `1e5`, `1e6`), same `(N,TMAX)` as assignment; **`doc/figures/hw4_q4_Nu_vs_t_Ra_sweep.pdf`** from the three CSVs |
| **§4(c)** | Table `N`=16–128; discussion; plot | Four runs at `Ra=1e4` only; **`doc/figures/hw4_q4_Nu_vs_N_Ra1e4.pdf`**; compare median/last `Nu` to **4.884** |
| **§5(1)** | Bonus sweep table (+ heatmap PDF) | `./scripts/run_bonus.sh sweep …` → `hw4/python/output/bonus_convergence/summary_nu_final.csv`; plot → **`hw4_bonus_Nu_mesh_Ra_heatmap.pdf`** |
| **§5(2)** | _(optional)_ | Leaving outline + “not implemented” is fine unless you author new C TS code |
| **§5(3)** | Paragraph + overlay figure | `./scripts/run_bonus.sh cn …` vs matching `run_q4.sh`; **`hw4_bonus_Nu_bdf_vs_cn.pdf`** |

**LaTeX path note:** appendix figures use `\texttt{figures/hw4_*.png}` — **`pdflatex` from `hw4/doc/`** (see `doc/README.md`) so `figures/` resolves to **`doc/figures/`**.

After exports, replace each `\placeholdfig{…}{…}{figures/xxx.pdf}{…}` with a normal `figure` environment and `\includegraphics[width=…]{figures/xxx.pdf}`.

**Quick minimum (required problems only):** finish **§2(b–c)** from logs → **§3** one PNG/PDF → **§4(a–c)** from `run_q4.sh` sweeps/plots → `pdflatex` `doc/solutions.tex`.
