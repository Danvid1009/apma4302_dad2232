#!/usr/bin/env python3
"""Compute Poisson spatial convergence and write CSV."""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

import numpy as np
from petsc4py import PETSc

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from vtk_export import gather_vec_to_grid  # noqa: E402


def build_system(nx: int, ny: int, comm: PETSc.Comm) -> tuple[PETSc.Mat, PETSc.Vec, PETSc.Vec]:
    n = nx * ny
    A = PETSc.Mat().createAIJ([n, n], nnz=5, comm=comm)
    A.setUp()
    b = PETSc.Vec().createMPI(n, comm=comm)
    x = PETSc.Vec().createMPI(n, comm=comm)

    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)
    hx2 = 1.0 / (hx * hx)
    hy2 = 1.0 / (hy * hy)

    r0, r1 = A.getOwnershipRange()
    for row in range(r0, r1):
        j = row // nx
        i = row % nx
        if i == 0 or j == 0 or i == nx - 1 or j == ny - 1:
            A.setValue(row, row, 1.0)
            b.setValue(row, 0.0)
            continue
        A.setValue(row, row, 2.0 * (hx2 + hy2))
        A.setValue(row, row - 1, -hx2)
        A.setValue(row, row + 1, -hx2)
        A.setValue(row, row - nx, -hy2)
        A.setValue(row, row + nx, -hy2)
        xx = i * hx
        yy = j * hy
        b.setValue(row, 2.0 * np.pi**2 * np.sin(np.pi * xx) * np.sin(np.pi * yy))
    A.assemble()
    b.assemble()
    x.set(0.0)
    return A, b, x


def solve_once(nx: int, ny: int, ksp_type: str, pc_type: str) -> tuple[int, float, float, float]:
    comm = PETSc.COMM_WORLD
    A, b, x = build_system(nx, ny, comm)
    ksp = PETSc.KSP().create(comm)
    ksp.setOperators(A)
    ksp.setType(ksp_type)
    ksp.getPC().setType(pc_type)
    ksp.setFromOptions()

    t0 = time.perf_counter()
    ksp.solve(b, x)
    elapsed = time.perf_counter() - t0

    grid = gather_vec_to_grid(x, nx, ny)
    if comm.getRank() != 0:
        return 0, 0.0, 0.0, elapsed

    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)
    xs = np.linspace(0.0, 1.0, nx)
    ys = np.linspace(0.0, 1.0, ny)
    xx, yy = np.meshgrid(xs, ys)
    exact = np.sin(np.pi * xx) * np.sin(np.pi * yy)
    diff = grid - exact
    err_l2 = np.sqrt(np.sum(diff * diff) * hx * hy)
    err_inf = np.max(np.abs(diff))
    return ksp.getIterationNumber(), float(err_l2), float(err_inf), elapsed


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Poisson convergence sweep.")
    parser.add_argument("--sizes", nargs="+", type=int, default=[21, 41, 81, 121, 161])
    parser.add_argument("--ksp-type", default="cg")
    parser.add_argument("--pc-type", default="jacobi")
    parser.add_argument("--out", default="output/poisson_convergence.csv")
    args = parser.parse_args()

    comm = PETSc.COMM_WORLD
    rank = comm.getRank()

    rows: list[tuple[int, float, float, float, int]] = []
    for n in args.sizes:
        its, l2, infn, sec = solve_once(n, n, args.ksp_type, args.pc_type)
        if rank == 0:
            rows.append((n, l2, infn, sec, its))
            print(f"n={n} l2={l2:.6e} inf={infn:.6e} t={sec:.4f}s iters={its}")

    if rank == 0:
        out = pathlib.Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8") as fh:
            fh.write("n,error_l2,error_inf,solve_time_s,iterations\n")
            for n, l2, infn, sec, its in rows:
                fh.write(f"{n},{l2:.12e},{infn:.12e},{sec:.8e},{its}\n")
        print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
