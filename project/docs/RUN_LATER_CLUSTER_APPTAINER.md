# Run later: cluster, Slurm, Apptainer (optional)

**Default workflow for this Monte Carlo project:** run everything **locally** with your PETSc-enabled Python (`mpirun`, `petsc4py`, scripts under `scripts/`). That is enough for correctness, convergence, barrier sweeps, error surfaces, antithetic benches, and moderate strong/weak scaling on your own machine.

This file is **only** for when you want extra HPC “production” flavor (batch jobs, shared systems, or a container you already use elsewhere).

---

## Your environment (notes you gave)

- You have **Firedrake** and **PETSc TS** (time stepping / the TS stack) inside an **Apptainer** image.
- You have access to **very few nodes** on the shared system.

**Implication:** treat cluster runs as **small, single-node (or very low node-count) jobs**. The report narrative can still discuss MPI collectives and scaling; you are measuring **algorithm + local parallelism**, not a huge multi-node campaign.

---

## What to defer here (optional checklist)

- [ ] Edit `scripts/submit_scaling.slurm` for your site (modules, walltime, `--ntasks-per-node` capped to what you are allocated).
- [ ] Run `run_experiment_matrix.py` **once** with a **short** `--ranks` list (e.g. `1 2 4`) so one job fits a small allocation.
- [ ] Save `slurm-*.out` / CSVs under `output/` and append rows to `output/registry/experiment_runs.csv` via `scripts/log_experiment_run.py`.

---

## Slurm template (already in repo)

- `scripts/submit_scaling.slurm` — example: single node, 8 tasks, calls `run_scaling.py` then `plot_scaling.py`.
- **Few nodes:** keep `#SBATCH --nodes=1` and reduce `--ntasks-per-node` to match what you can actually get (e.g. 4). Match `--ranks` inside `run_scaling.py` to the same cap so `mpirun` is not oversubscribed.

After edits:

```bash
sbatch scripts/submit_scaling.slurm
```

Tail progress: `tail -f output/slurm-<jobid>.out` (path may differ if your site rewrites `#SBATCH --output`).

---

## Apptainer (if you run this project inside the same image as Firedrake)

Typical pattern (you must adjust image path and binds):

```bash
apptainer exec --bind /path/to/apma4302_dad2232/project:/work \
  /path/to/your.sif \
  bash -lc 'cd /work/project && mpirun -n 4 python scripts/run_pricing.py --paths 100000'
```

**Notes:**

- Use the **Python and MPI** inside the container so `petsc4py` matches the linked PETSc.
- Bind-mount the **git checkout** (this `project` directory or repo root) so outputs land on the host.
- This repo does **not** require Firedrake for the Monte Carlo drivers; the image is just a convenient MPI+PETSc runtime if that is where your toolchain lives.

---

## SSH interactive (optional)

Same commands as `README.md`, from a login node or a `salloc` shell, with the same Python/MPI modules (or container) you use for batch jobs.

---

## Registry

When you eventually run any of the above, log artifacts with:

```bash
python scripts/log_experiment_run.py --key slurm_scaling_may --notes "1 node 4 ranks" --output output/strong_scaling.csv
```

Keep **`docs/HPC_EXPERIMENT_REGISTRY.md`** and **`RESULTS_PHASE1.md`** focused on **where the numbers came from** (local vs cluster); one sentence is enough.
