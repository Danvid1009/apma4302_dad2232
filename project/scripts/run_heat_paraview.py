#!/usr/bin/env python3
"""Implicit-Euler heat equation solve with ParaView VTP/PVD export."""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
from petsc4py import PETSc

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from vtk_export import gather_vec_to_grid, write_grid_surface_vtp, write_pvd  # noqa: E402


def idx(i: int, j: int, nx: int) -> int:
    return j * nx + i


def build_heat_matrix(nx: int, ny: int, dt: float, kappa: float, comm: PETSc.Comm) -> PETSc.Mat:
    n = nx * ny
    A = PETSc.Mat().createAIJ([n, n], nnz=5, comm=comm)
    A.setUp()

    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)
    cx = dt * kappa / (hx * hx)
    cy = dt * kappa / (hy * hy)

    r0, r1 = A.getOwnershipRange()
    for row in range(r0, r1):
        j = row // nx
        i = row % nx
        if i == 0 or j == 0 or i == nx - 1 or j == ny - 1:
            A.setValue(row, row, 1.0)
            continue

        A.setValue(row, row, 1.0 + 2.0 * cx + 2.0 * cy)
        A.setValue(row, idx(i - 1, j, nx), -cx)
        A.setValue(row, idx(i + 1, j, nx), -cx)
        A.setValue(row, idx(i, j - 1, nx), -cy)
        A.setValue(row, idx(i, j + 1, nx), -cy)
    A.assemble()
    return A


def initialize_state(vec: PETSc.Vec, nx: int, ny: int) -> None:
    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)
    r0, r1 = vec.getOwnershipRange()
    for row in range(r0, r1):
        j = row // nx
        i = row % nx
        x = i * hx
        y = j * hy
        # Smooth bump plus sinusoidal mode for colorful later slices.
        val = np.sin(np.pi * x) * np.sin(np.pi * y) + 0.35 * np.exp(-80.0 * ((x - 0.3) ** 2 + (y - 0.7) ** 2))
        if i == 0 or j == 0 or i == nx - 1 or j == ny - 1:
            val = 0.0
        vec.setValue(row, val)
    vec.assemble()


def enforce_boundary_zero(vec: PETSc.Vec, nx: int, ny: int) -> None:
    r0, r1 = vec.getOwnershipRange()
    for row in range(r0, r1):
        j = row // nx
        i = row % nx
        if i == 0 or j == 0 or i == nx - 1 or j == ny - 1:
            vec.setValue(row, 0.0)
    vec.assemble()


def main() -> int:
    parser = argparse.ArgumentParser(description="Heat equation ParaView exporter.")
    parser.add_argument("--nx", type=int, default=121)
    parser.add_argument("--ny", type=int, default=121)
    parser.add_argument("--dt", type=float, default=5e-4)
    parser.add_argument("--steps", type=int, default=250)
    parser.add_argument("--kappa", type=float, default=0.15)
    parser.add_argument("--outdir", default="output/paraview/vtp")
    parser.add_argument(
        "--series-name",
        default="heat",
        help="File name stem: {series_name}_NNNN.vtp, {series_name}.pvd, {series_name}_solver_history.csv",
    )
    parser.add_argument("--ksp-type", default="cg")
    parser.add_argument("--pc-type", default="jacobi")
    parser.add_argument("--save-every", type=int, default=25)
    args = parser.parse_args()

    comm = PETSc.COMM_WORLD
    rank = comm.getRank()
    n = args.nx * args.ny

    A = build_heat_matrix(args.nx, args.ny, args.dt, args.kappa, comm)
    u = PETSc.Vec().createMPI(n, comm=comm)
    rhs = PETSc.Vec().createMPI(n, comm=comm)
    initialize_state(u, args.nx, args.ny)

    ksp = PETSc.KSP().create(comm)
    ksp.setOperators(A)
    ksp.setType(args.ksp_type)
    ksp.getPC().setType(args.pc_type)
    ksp.setFromOptions()

    outdir = pathlib.Path(args.outdir)
    dx = 1.0 / (args.nx - 1)
    dy = 1.0 / (args.ny - 1)
    pvd_entries: list[tuple[float, str]] = []
    residual_rows: list[tuple[float, int, float]] = []

    def save_snapshot(step: int, time: float) -> None:
        grid = gather_vec_to_grid(u, args.nx, args.ny)
        if rank != 0:
            return
        fname = f"{args.series_name}_{step:04d}.vtp"
        write_grid_surface_vtp(outdir / fname, grid, dx, dy, "temperature")
        pvd_entries.append((time, fname))

    save_snapshot(0, 0.0)

    for step in range(1, args.steps + 1):
        rhs.copy(u)
        enforce_boundary_zero(rhs, args.nx, args.ny)
        ksp.solve(rhs, u)
        residual_rows.append((step * args.dt, ksp.getIterationNumber(), float(ksp.getResidualNorm())))
        if step % args.save_every == 0 or step == args.steps:
            save_snapshot(step, step * args.dt)

    if rank == 0:
        pvd_path = outdir / f"{args.series_name}.pvd"
        hist_path = outdir / f"{args.series_name}_solver_history.csv"
        write_pvd(pvd_path, pvd_entries)
        with hist_path.open("w", encoding="utf-8") as fh:
            fh.write("time,iterations,residual_norm\n")
            for t, itn, rn in residual_rows:
                fh.write(f"{t:.8e},{itn},{rn:.12e}\n")
        print(f"Wrote {pvd_path}")
        print(f"Wrote {len(pvd_entries)} VTP surface snapshots")
        print(f"Wrote {hist_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
