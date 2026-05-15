#!/usr/bin/env python3
"""Run one timed Poisson solve and append row for scaling study."""

from __future__ import annotations

import argparse
import pathlib
import time

import numpy as np
from petsc4py import PETSc


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


def main() -> int:
    parser = argparse.ArgumentParser(description="Timed Poisson solve for scaling CSV.")
    parser.add_argument("--nx", type=int, default=401)
    parser.add_argument("--ny", type=int, default=401)
    parser.add_argument("--ksp-type", default="cg")
    parser.add_argument("--pc-type", default="jacobi")
    parser.add_argument("--out", default="output/strong_scaling.csv")
    parser.add_argument("--label", default="poisson2d")
    args = parser.parse_args()

    comm = PETSc.COMM_WORLD
    rank = comm.getRank()
    ranks = comm.getSize()

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
    elapsed = time.perf_counter() - t0

    if rank == 0:
        out = pathlib.Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        header = not out.exists()
        with out.open("a", encoding="utf-8") as fh:
            if header:
                fh.write("label,nx,ny,ranks,solve_time_s,iterations,residual_norm,ksp_type,pc_type\n")
            fh.write(
                f"{args.label},{args.nx},{args.ny},{ranks},{elapsed:.8e},{ksp.getIterationNumber()},"
                f"{ksp.getResidualNorm():.12e},{ksp.getType()},{ksp.getPC().getType()}\n"
            )
        print(f"Appended row to {out}")
        print(f"ranks={ranks} time={elapsed:.4f}s iters={ksp.getIterationNumber()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
