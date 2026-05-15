#!/usr/bin/env python3
"""
Sweep **KSP types** on the same 2D Poisson system (fixed grid, fixed PC).

Course tie-in: Krylov subspace methods (CG for SPD-ish Laplacian vs GMRES /
BiCGStab as general templates). Use a **cheap PC** (default ``jacobi``) so the
comparison stresses the *outer* Krylov iteration counts and timings.

Example::

    mpiexec -n 4 python scripts/run_poisson_ksp_sweep.py --nx 65 --ny 65 \\
      --pc-type jacobi --ksp-types cg gmres bcgsl
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
    p = argparse.ArgumentParser(description="Poisson KSP: sweep KSP types, append CSV.")
    p.add_argument("--nx", type=int, default=65)
    p.add_argument("--ny", type=int, default=65)
    p.add_argument("--pc-type", type=str, default="jacobi", help="Fixed PC for all KSP trials.")
    p.add_argument(
        "--ksp-types",
        type=str,
        nargs="+",
        default=["cg", "gmres", "bcgsl"],
        help="PETSc KSP types to try in order.",
    )
    p.add_argument("--out", type=str, default="output/poisson_ksp_sweep.csv")
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

    for kname in args.ksp_types:
        A, b, x = build_system(args.nx, args.ny, comm)
        ksp = PETSc.KSP().create(comm)
        ksp.setOperators(A)
        ok = True
        err_msg = ""
        try:
            ksp.setType(kname)
        except Exception as e:
            ok = False
            err_msg = f"setType failed: {e}"
        if ok:
            try:
                ksp.getPC().setType(args.pc_type)
            except Exception as e:
                ok = False
                err_msg = f"PC setType failed: {e}"
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
                "ksp_type": kname,
                "pc_type": args.pc_type,
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
            print(f"ksp={kname} time={elapsed:.4f}s iters={iters} {status}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
