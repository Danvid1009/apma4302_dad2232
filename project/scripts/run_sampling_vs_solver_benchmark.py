#!/usr/bin/env python3
"""
Compare wall time: **one Poisson KSP solve** (Mat/Vec + Krylov + PC) vs **one
European MC pricer** on the same MPI communicator.

Course tie-in: contrasts **solver-dominated PDE** mini-app with **sampling-dominated**
Monte Carlo (same PETSc runtime, different numerical bottleneck). Uses
``build_system`` from ``run_poisson_timed.py`` (structured 5-point Laplacian).

Example::

    mpiexec -n 4 python scripts/run_sampling_vs_solver_benchmark.py \\
      --nx 129 --ny 129 --paths 80000 --steps 252
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import pathlib
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from petsc4py import PETSc

from montecarlo.engine import petsc_monte_carlo_price
from montecarlo.payoffs import european_call_payoff


def _load_poisson_builder():
    poisson_path = Path(__file__).resolve().parent / "run_poisson_timed.py"
    spec = importlib.util.spec_from_file_location("run_poisson_timed", poisson_path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod.build_system


def parse_args():
    p = argparse.ArgumentParser(description="Poisson KSP vs Monte Carlo timing on same ranks.")
    p.add_argument("--nx", type=int, default=129)
    p.add_argument("--ny", type=int, default=129)
    p.add_argument("--ksp-type", type=str, default="cg")
    p.add_argument("--pc-type", type=str, default="jacobi")
    p.add_argument("--paths", type=int, default=50000)
    p.add_argument("--steps", type=int, default=252)
    p.add_argument("--s0", type=float, default=100.0)
    p.add_argument("--strike", type=float, default=100.0)
    p.add_argument("--r", type=float, default=0.05)
    p.add_argument("--sigma", type=float, default=0.2)
    p.add_argument("--t", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=str, default="output/sampling_vs_solver_benchmark.csv")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    build_system = _load_poisson_builder()
    comm = PETSc.COMM_WORLD
    rank = comm.getRank()
    size = comm.getSize()

    A, b, x = build_system(args.nx, args.ny, comm)
    ksp = PETSc.KSP().create(comm)
    ksp.setOperators(A)
    ksp.setType(args.ksp_type)
    ksp.getPC().setType(args.pc_type)
    ksp.setFromOptions()

    comm.barrier()
    t0 = time.perf_counter()
    ksp.solve(b, x)
    comm.barrier()
    t_poisson = time.perf_counter() - t0
    poisson_iters = ksp.getIterationNumber()

    payoff_builder = lambda strike: lambda paths: european_call_payoff(paths, strike)
    comm.barrier()
    t1 = time.perf_counter()
    res = petsc_monte_carlo_price(
        args.s0,
        args.strike,
        args.r,
        args.sigma,
        args.t,
        args.steps,
        args.paths,
        payoff_builder=payoff_builder,
        seed=args.seed,
    )
    comm.barrier()
    t_mc = time.perf_counter() - t1

    ksp.destroy()
    A.destroy()
    b.destroy()
    x.destroy()

    if rank == 0:
        out = pathlib.Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        row = {
            "ranks": size,
            "nx": args.nx,
            "ny": args.ny,
            "ksp_type": args.ksp_type,
            "pc_type": args.pc_type,
            "poisson_time_s": f"{t_poisson:.8e}",
            "poisson_iters": poisson_iters,
            "mc_paths": args.paths,
            "mc_steps": args.steps,
            "mc_time_s": f"{t_mc:.8e}",
            "mc_price": f"{res['price']:.10f}",
            "ratio_mc_over_poisson": f"{(t_mc / max(t_poisson, 1e-30)):.4f}",
        }
        write_header = not out.exists()
        with out.open("a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(row.keys()))
            if write_header:
                w.writeheader()
            w.writerow(row)
        print(
            f"ranks={size} poisson_s={t_poisson:.4f} iters={poisson_iters} "
            f"mc_s={t_mc:.4f} ratio_mc/poisson={float(row['ratio_mc_over_poisson']):.2f} "
            f"mc_price={res['price']:.6f} wrote {out}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
