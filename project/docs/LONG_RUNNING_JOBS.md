# Long-running jobs and “is this stuck?” checks

Use this when a shell was **backgrounded** (Cursor agent, `&`, Slurm, or `nohup`) and you want evidence of **forward progress**, not just a silent prompt.

## 1) Cursor integrated terminals (this repo)

Cursor mirrors each terminal to a text file under your machine’s project metadata folder (name includes the workspace path). Typical layout:

`~/.cursor/projects/<workspace-slug>/terminals/*.txt`

Each file has a small YAML-ish header (`pid`, `cwd`, `command`, `started_at`, sometimes `running_for_ms`) and then the live transcript.

**Quick manual check**

```bash
ls -lt ~/.cursor/projects/*/terminals/*.txt | head
tail -n 30 ~/.cursor/projects/<your-slug>/terminals/692916.txt
```

If the transcript **grows** (new lines, new timestamps) or the command finishes and you see a footer with `exit_code`, it is not stuck. If the transcript is frozen for a long time *and* the workload should be chatty, investigate.

**Automated probe (recommended)**

From the repository root (any cwd works if you pass `--dir`):

```bash
python scripts/terminal_probe.py --auto
```

This finds terminal files whose recorded `cwd` matches this Git repo, prints whether the header `pid` still exists, file age/size, and the tail of the log. It also scans for an `exit_code` footer so you can tell **running vs finished** even when the UI still shows “busy.”

**Interpreting “stuck” vs “slow”**

- **`pid_alive=True`** only means the OS still has a process with that id. A hung PETSc/MPI or blocked I/O job can look alive forever with **no new transcript lines**.
- Prefer **two probes a minute apart**: if `mtime`/`size` never change *and* the tail never grows, treat it as stuck (or blocked on network/filesystem) even when `pid_alive=True`.
- If `pid_alive=False` but there is still no `exit_code` footer, the shell snapshot may be stale; re-open the terminal tab or check `ps` manually.
- If you confirm a **zombie** agent command (PID alive, no progress, known-bad command), you can end it with `kill <pid>` from another shell; only do that when you are sure it is not someone else’s session.

Override the folder explicitly if needed:

```bash
export CURSOR_TERMINALS_DIR="$HOME/.cursor/projects/Users-dan-Desktop-Columbia-HPC-4302-apma4302-dad2232-project/terminals"
python scripts/terminal_probe.py --dir "$CURSOR_TERMINALS_DIR"
```

## 2) PETSc / MPI jobs you start yourself

- **Unbuffered Python** so prints appear immediately: `python -u scripts/run_scaling.py ...` or set `PYTHONUNBUFFERED=1`.
- **Log everything**: `mpirun ... 2>&1 | tee run-$(date +%Y%m%d-%H%M%S).log` so you can `tail -f` the log in another window.
- **Occasional progress** for very long inner loops: print every \(k\) paths or ranks (keep it rate-limited so you do not drown I/O on clusters).

## 3) Cluster (Slurm) — optional

Bulk of this project is intended to run **locally** with PETSc/MPI. Slurm is only if you choose a shared system; see **[`RUN_LATER_CLUSTER_APPTAINER.md`](RUN_LATER_CLUSTER_APPTAINER.md)** for templates and few-node tips.

- `squeue -u $USER` — job still pending/running?
- `tail -f slurm-<jobid>.out` — output advancing?
- If walltime expired you will see `CANCELLED` / `TIMEOUT` in the accounting record; keep that `.out` path in `experiment_runs.csv` via `scripts/log_experiment_run.py`.

## 4) Optional wrapper: wall-clock heartbeat

If a third-party command is completely silent, wrap it so you still see time advancing:

```bash
bash scripts/with_heartbeat.sh 60 mpirun -n 8 python -u scripts/run_pricing.py --paths 2000000
```

The first argument is the interval in seconds between `date` lines on stderr; the rest is the real command.
