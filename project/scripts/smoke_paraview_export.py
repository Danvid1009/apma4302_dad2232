#!/usr/bin/env python3
"""PETSc Poisson smoke test with direct ParaView VTP/PVD export."""

from __future__ import annotations

import argparse
import pathlib
import sys
from typing import Tuple

import numpy as np
from petsc4py import PETSc

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from vtk_export import write_grid_surface_vtp, write_pvd  # noqa: E402


def build_system(nx: int, ny: int, comm: PETSc.Comm) -> Tuple[PETSc.Mat, PETSc.Vec, PETSc.Vec]:
    n = nx * ny
    A = PETSc.Mat().createAIJ([n, n], nnz=5, comm=comm)
    A.setUp()
    b = PETSc.Vec().createMPI(n, comm=comm)
    x = PETSc.Vec().createMPI(n, comm=comm)

    istart, iend = A.getOwnershipRange()
    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)
    hx2 = 1.0 / (hx * hx)
    hy2 = 1.0 / (hy * hy)

    for row in range(istart, iend):
        j = row // nx
        i = row % nx
        on_boundary = i == 0 or j == 0 or i == nx - 1 or j == ny - 1

        if on_boundary:
            A.setValue(row, row, 1.0)
            b.setValue(row, 0.0)
            continue

        center = 2.0 * (hx2 + hy2)
        west = row - 1
        east = row + 1
        south = row - nx
        north = row + nx

        A.setValue(row, row, center)
        A.setValue(row, west, -hx2)
        A.setValue(row, east, -hx2)
        A.setValue(row, south, -hy2)
        A.setValue(row, north, -hy2)

        xcoord = i * hx
        ycoord = j * hy
        rhs = 2.0 * np.pi**2 * np.sin(np.pi * xcoord) * np.sin(np.pi * ycoord)
        b.setValue(row, rhs)

    A.assemble()
    b.assemble()
    x.set(0.0)
    return A, b, x


def solve_system(A: PETSc.Mat, b: PETSc.Vec, x: PETSc.Vec, ksp_type: str, pc_type: str) -> PETSc.KSP:
    ksp = PETSc.KSP().create(A.getComm())
    ksp.setOperators(A)
    ksp.setType(ksp_type)
    ksp.getPC().setType(pc_type)
    ksp.setFromOptions()
    ksp.solve(b, x)
    return ksp


def gather_solution(x: PETSc.Vec, nx: int, ny: int) -> np.ndarray:
    comm = x.getComm()
    rank = comm.getRank()
    arr_local = x.getArray(readonly=True).copy()
    gathered = comm.tompi4py().gather(arr_local, root=0)
    if rank != 0:
        return np.empty((0, 0), dtype=np.float64)

    full = np.concatenate(gathered)
    return full.reshape((ny, nx))


def main() -> int:
    parser = argparse.ArgumentParser(description="PETSc Poisson smoke test with ParaView output.")
    parser.add_argument("--nx", type=int, default=101, help="Grid points in x")
    parser.add_argument("--ny", type=int, default=101, help="Grid points in y")
    parser.add_argument("--ksp-type", default="cg", help="PETSc KSP type (default: cg)")
    parser.add_argument("--pc-type", default="jacobi", help="PETSc PC type (default: jacobi)")
    parser.add_argument(
        "--outdir",
        default="output/paraview/vtp",
        help="Output directory for VTP/PVD files (default: single shared ParaView folder)",
    )
    parser.add_argument(
        "--basename",
        default="poisson",
        help="Base name for outputs: {basename}_0000.vtp, {basename}.pvd, {basename}_solver_log.txt",
    )
    args = parser.parse_args()

    comm = PETSc.COMM_WORLD
    rank = comm.getRank()

    if args.nx < 3 or args.ny < 3:
        raise ValueError("nx and ny must be at least 3.")

    A, b, x = build_system(args.nx, args.ny, comm)
    ksp = solve_system(A, b, x, args.ksp_type, args.pc_type)

    grid = gather_solution(x, args.nx, args.ny)
    if rank == 0:
        outdir = pathlib.Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        base = args.basename
        vtp_name = f"{base}_0000.vtp"
        vtp_path = outdir / vtp_name
        pvd_path = outdir / f"{base}.pvd"

        dx = 1.0 / (args.nx - 1)
        dy = 1.0 / (args.ny - 1)
        write_grid_surface_vtp(vtp_path, grid, dx, dy, "u")
        write_pvd(pvd_path, [(0.0, vtp_name)])

        log_path = outdir / f"{base}_solver_log.txt"
        with log_path.open("w", encoding="utf-8") as fh:
            fh.write(f"ksp_type={ksp.getType()}\n")
            fh.write(f"pc_type={ksp.getPC().getType()}\n")
            fh.write(f"iterations={ksp.getIterationNumber()}\n")
            fh.write(f"residual_norm={ksp.getResidualNorm():.12e}\n")
            fh.write(f"converged_reason={ksp.getConvergedReason()}\n")

        print(f"Wrote {vtp_path}")
        print(f"Wrote {pvd_path}")
        print(f"Wrote {log_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
