#!/usr/bin/env python3
"""
Sweep **preconditioner types** on the same Poisson system (fixed grid).

Course tie-in: preconditioned Krylov methods + PC choices (Jacobi vs GAMG,
etc.). Appends one CSV row per ``--pc-type``. Some PCs may be unavailable on
minimal PETSc builds; failures print a warning and skip.

Example::

    mpiexec -n 4 python scripts/run_poisson_pc_sweep.py --nx 65 --ny 65 \\
      --pc-types jacobi sor gamg
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import pathlib
import time
from pathlib import Path

from petsc4py import PETSc


def _load_poisson_builder():
    poisson_path = Path(__file__).resolve().parent / "run_poisson_timed.py"
    spec = importlib.util.spec_from_file_location("run_poisson_timed", poisson_path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod.build_system


def parse_args():
    p = argparse.ArgumentParser(description="Poisson KSP: sweep PC types, append CSV.")
    p.add_argument("--nx", type=int, default=65)
    p.add_argument("--ny", type=int, default=65)
    p.add_argument("--ksp-type", type=str, default="cg")
    p.add_argument(
        "--pc-types",
        type=str,
        nargs="+",
        default=["jacobi", "sor"],
        help="PETSc PC types to try in order.",
    )
    p.add_argument("--out", type=str, default="output/poisson_pc_sweep.csv")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    build_system = _load_poisson_builder()
    comm = PETSc.COMM_WORLD
    rank = comm.getRank()
    size = comm.getSize()

    out = pathlib.Path(args.out)
    if rank == 0:
        out.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "ranks",
        "nx",
        "ny",
        "ksp_type",
        "pc_type",
        "solve_time_s",
        "iterations",
        "residual_norm",
        "ok",
        "error",
    ]

    for pc in args.pc_types:
        A, b, x = build_system(args.nx, args.ny, comm)
        ksp = PETSc.KSP().create(comm)
        ksp.setOperators(A)
        ksp.setType(args.ksp_type)
        pc_obj = ksp.getPC()
        ok = True
        err_msg = ""
        try:
            pc_obj.setType(pc)
        except Exception as e:
            ok = False
            err_msg = f"setType failed: {e}"
        if ok:
            try:
                ksp.setFromOptions()
            except Exception as e:
                ok = False
                err_msg = f"setFromOptions failed: {e}"
        elapsed = 0.0
        iters = -1
        resn = 0.0
        if ok:
            x.set(0.0)
            comm.barrier()
            t0 = time.perf_counter()
            try:
                ksp.solve(b, x)
            except Exception as e:
                ok = False
                err_msg = f"solve failed: {e}"
            comm.barrier()
            elapsed = time.perf_counter() - t0
            iters = ksp.getIterationNumber()
            resn = ksp.getResidualNorm() if ok else 0.0

        ksp.destroy()
        A.destroy()
        b.destroy()
        x.destroy()

        if rank == 0:
            row = {
                "ranks": size,
                "nx": args.nx,
                "ny": args.ny,
                "ksp_type": args.ksp_type,
                "pc_type": pc,
                "solve_time_s": f"{elapsed:.8e}",
                "iterations": iters,
                "residual_norm": f"{resn:.12e}",
                "ok": ok,
                "error": err_msg,
            }
            write_header = not out.exists()
            with out.open("a", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fieldnames)
                if write_header:
                    w.writeheader()
                w.writerow(row)
            status = "ok" if ok else f"FAIL {err_msg}"
            print(f"pc={pc} time={elapsed:.4f}s iters={iters} {status}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
