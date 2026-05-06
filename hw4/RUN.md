# Homework 4 — pull, environment, run order

Use this on a **login node, compute node, SSH session, or in a Docker/app container** as long as the dependencies for each part are available there.

**Layout:** write-up lives under **`hw4/doc/`** (`solutions.tex`, `hw4.pdf`, `figures/`). Code under **`hw4/c/`**, **`hw4/python/`**, runners under **`hw4/scripts/`**. See `hw4/README.md` and `hw4/doc/SCRIPTS.md`.

## 1) When the repository is “ready to push”

You (or whoever has write access) should **commit and push** after:

- `hw4/python/biharm.py`, `biharm_temperature_rhs.py`, `convection.py` are present.
- `hw4/scripts/run_q2.sh`, `run_q3.sh`, `run_q4.sh` are present and executable (`chmod +x hw4/scripts/*.sh`).
- This file `hw4/RUN.md` is present.

Then anyone can `git pull` and follow the sections below.

---

## 2) Clone or pull the repository (SSH / container)

Replace URL and path with yours.

```bash
# first time
git clone <YOUR_GIT_REMOTE_URL>
cd <your-repo-directory>

# later updates
cd <your-repo-directory>
git pull
```

All homework paths below assume you are at the **repository root** or use absolute paths.

---

## 3) Software you need (by part)

| Part | What you need |
|------|----------------|
| **Q2 — C `biharm`** | PETSc configured with **MUMPS** (and usual PETSc build). Set `PETSC_DIR` and `PETSC_ARCH`. `make` in `hw4/c`. |
| **Q2 — Python `biharm.py`** | **Firedrake** (not necessarily `firedrake-ts`; steady solve only). **MUMPS** visible to PETSc inside that environment. |
| **Q3 — `biharm_temperature_rhs.py`** | Same as Q2 Firedrake (steady solve + VTK). |
| **Q4 — `convection.py`** | **Firedrake + firedrake-ts**, **MUMPS** for the monolithic LU per step. Enough wall time for large `--t-max`. |

Activate your course / lab Firedrake environment the same way you do for other homework (e.g. `source firedrake/bin/activate` or module load + venv — **use whatever your container documents**).

### Cluster: Apptainer / Singularity (`firedrake-ts.sif`)

Host `python3` usually has **no** `petsc4py` / Firedrake. The `run_q2.sh` / `run_q3.sh` / `run_q4.sh` scripts can run **`python3` inside the image** when you set:

```bash
export HW4_APPTAINER_SIF=/path/to/firedrake-ts.sif
export HW4_APPTAINER_BIND=/path/to/your/work:/path/to/your/work
```

Optional:

- `FIREDRAKE_TS_SIF` — used if `HW4_APPTAINER_SIF` is empty (same path).
- `HW4_APPTAINER_RUNNER` — default: `apptainer`, else `singularity` if only that exists.
- `HW4_APPTAINER_EXTRA_ARGS` — e.g. `--nv` passed before `--bind`.
- `HW4_APPTAINER_UNSET_SLURM=1` — unset all `SLURM_*` in the script before `apptainer exec` if Open MPI still tries Slurm PMI and crashes.

The scripts set `OMPI_MCA_plm=isolated` by default when using the SIF (helps under interactive `srun`).

---

## 4) Recommended run order

From the repo root:

```bash
chmod +x hw4/scripts/run_q2.sh hw4/scripts/run_q3.sh hw4/scripts/run_q4.sh
```

### Q2 — C + Firedrake biharmonic drivers

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=your-arch   # if your PETSc install uses PETSC_ARCH

hw4/scripts/run_q2.sh
```

**Python-only on a cluster** (no PETSc for C): set the Apptainer variables in §3b, then:

```bash
SKIP_HW4_C=1 hw4/scripts/run_q2.sh
```

- C logs: `hw4/output/logs/q2_c_*.log`
- Python logs: `hw4/output/logs/q2_py_*.log`
- Note: the script runs three Firedrake presets in a row; each overwrites `output/q2_firedrake/biharm.pvd` unless you copy/rename between runs.

**Manual C timing (matches assignment wording):**

```bash
cd hw4/c
./biharm -options_file options_file_direct -log_view | grep SNESSolve
```

**Manual Firedrake (one preset):**

```bash
cd hw4/python
python3 biharm.py --preset direct          # or split_direct | split_mg
```

### Q3 — temperature RHS + ParaView figure

```bash
hw4/scripts/run_q3.sh
# optional mesh controls:
# hw4/scripts/run_q3.sh --n-base 2 --levels 8 --A 0.1
```

Open `hw4/output/q3/biharm_T_rhs.pvd` in ParaView (copy to your laptop if the run was on a cluster).

### Q4 — convection + Nusselt number CSV

Quick sanity (short time, no VTK):

```bash
RA=1e2 N=32 TMAX=50 DT=0.1 VTK_EVERY=0 hw4/scripts/run_q4.sh
```

Assignment-style long run (example; may need **batch** job and hours of wall time):

```bash
RA=1e2 N=64 TMAX=100000 DT=0.1 VTK_EVERY=0 hw4/scripts/run_q4.sh
```

With occasional VTK:

```bash
RA=1e4 N=64 TMAX=5000 DT=0.1 VTK_EVERY=500 HW4_OUT=hw4/output/q4_try hw4/scripts/run_q4.sh
```

Outputs:

- `nu_history.csv` in the run directory (`HW4_OUT` or default under `hw4/output/...`)
- Optional `convection.pvd` if `--vtk-every` / `VTK_EVERY` > 0

**Mesh / Ra sweeps:** call `run_q4.sh` multiple times with different `RA`, `N`, `HW4_OUT`, or pass extra CLI args through to `convection.py`:

```bash
cd hw4/python
python3 convection.py --help
```

### Extra credit (bonus)

From `hw4/` with Firedrake (+ `firedrake-ts` for the sweep, since it calls `convection.py`):

```bash
chmod +x scripts/run_bonus.sh

# (1) Mesh × Ra sweep → python/output/bonus_convergence/summary_nu_final.csv
# Dry-run prints commands only:
./scripts/run_bonus.sh sweep --dry-run
# Real run (expensive: 12 DAE jobs by default). Increase --t-max for steadier Nu.
./scripts/run_bonus.sh sweep --t-max 8000 --nu-every 200

# (3) Crank–Nicolson + SNES per step (plain Firedrake; no firedrake-ts)
./scripts/run_bonus.sh cn --ra 1e4 --n 32 --t-max 500 --dt 0.1
```

Item **(2)** (PETSc C + `DMDA` + `TS` + fieldsplit) is not scripted here; it would be a new C driver alongside `biharm.c`.

---

## 5) If something fails

- **C build:** `PETSC_DIR` not set or wrong; missing MUMPS.
- **`import firedrake` fails:** Firedrake environment not activated in that shell/container.
- **`import firedrake_ts` fails (Q4 only):** install/activate the **firedrake-ts** add-on for that Firedrake build.
- **Nu looks wrong / NaN:** facet IDs for top/bottom may differ by mesh generator; try `python3 convection.py --ds-top 4 --ds-bot 2` (or swap) after checking your Firedrake `UnitSquareMesh` boundary marker convention.

---

## 6) For a colleague in SSH + Docker

1. Ensure the **same** Firedrake (+ `firedrake-ts` for Q4) image or bind-mounted install is documented for the container.
2. `git pull` inside the container workspace.
3. Activate Firedrake, set `PETSC_DIR` for C if running Q2 C, then run `hw4/scripts/run_q2.sh` → `run_q3.sh` → `run_q4.sh` in that order (or only the parts they need).
