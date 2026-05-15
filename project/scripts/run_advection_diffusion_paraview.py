#!/usr/bin/env python3
"""Advection-diffusion solve with PETSc and ParaView VTP/PVD export."""

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


def build_system_matrix(
    nx: int,
    ny: int,
    dt: float,
    nu: float,
    vel_x: float,
    vel_y: float,
    comm: PETSc.Comm,
) -> PETSc.Mat:
    n = nx * ny
    A = PETSc.Mat().createAIJ([n, n], nnz=5, comm=comm)
    A.setUp()

    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)
    cx = dt * nu / (hx * hx)
    cy = dt * nu / (hy * hy)
    ax = dt * vel_x / hx
    ay = dt * vel_y / hy

    r0, r1 = A.getOwnershipRange()
    for row in range(r0, r1):
        j = row // nx
        i = row % nx
        if i == 0 or j == 0 or i == nx - 1 or j == ny - 1:
            A.setValue(row, row, 1.0)
            continue

        # Backward Euler with upwind-like advection stencil for positive velocities.
        center = 1.0 + 2.0 * cx + 2.0 * cy + ax + ay
        west = -cx - ax
        east = -cx
        south = -cy - ay
        north = -cy

        A.setValue(row, row, center)
        A.setValue(row, idx(i - 1, j, nx), west)
        A.setValue(row, idx(i + 1, j, nx), east)
        A.setValue(row, idx(i, j - 1, nx), south)
        A.setValue(row, idx(i, j + 1, nx), north)
    A.assemble()
    return A


def set_initial_condition(vec: PETSc.Vec, nx: int, ny: int, ic: str) -> None:
    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)
    r0, r1 = vec.getOwnershipRange()
    for row in range(r0, r1):
        j = row // nx
        i = row % nx
        x = i * hx
        y = j * hy
        if ic == "twin":
            val = np.exp(-120.0 * ((x - 0.2) ** 2 + (y - 0.2) ** 2)) + 0.65 * np.exp(
                -95.0 * ((x - 0.72) ** 2 + (y - 0.38) ** 2)
            )
        else:
            val = np.exp(-120.0 * ((x - 0.2) ** 2 + (y - 0.2) ** 2))
        if i == 0 or j == 0 or i == nx - 1 or j == ny - 1:
            val = 0.0
        vec.setValue(row, val)
    vec.assemble()


def enforce_dirichlet(vec: PETSc.Vec, nx: int, ny: int) -> None:
    r0, r1 = vec.getOwnershipRange()
    for row in range(r0, r1):
        j = row // nx
        i = row % nx
        if i == 0 or j == 0 or i == nx - 1 or j == ny - 1:
            vec.setValue(row, 0.0)
    vec.assemble()


def main() -> int:
    parser = argparse.ArgumentParser(description="Advection-diffusion ParaView exporter.")
    parser.add_argument("--nx", type=int, default=121)
    parser.add_argument("--ny", type=int, default=121)
    parser.add_argument("--dt", type=float, default=8e-4)
    parser.add_argument("--steps", type=int, default=320)
    parser.add_argument("--nu", type=float, default=0.01)
    parser.add_argument("--vel-x", type=float, default=1.0)
    parser.add_argument("--vel-y", type=float, default=0.5)
    parser.add_argument("--outdir", default="output/paraview/vtp")
    parser.add_argument(
        "--series-name",
        default="advection_diffusion",
        help="File name stem for snapshots, PVD, and solver history CSV",
    )
    # Defaults safe under mpiexec on many PETSc builds (parallel ILU factor often unavailable).
    parser.add_argument("--ksp-type", default="cg")
    parser.add_argument("--pc-type", default="jacobi")
    parser.add_argument("--save-every", type=int, default=32)
    parser.add_argument(
        "--ic",
        choices=("single", "twin"),
        default="single",
        help="Initial concentration: single Gaussian blob, or twin blobs (richer ParaView frames).",
    )
    parser.add_argument(
        "--no-velocity-vti",
        action="store_true",
        help="Omit constant velocity vector field from VTP (smaller files; disables Glyph-on-velocity in ParaView).",
    )
    args = parser.parse_args()

    comm = PETSc.COMM_WORLD
    rank = comm.getRank()
    n = args.nx * args.ny
    outdir = pathlib.Path(args.outdir)
    dx = 1.0 / (args.nx - 1)
    dy = 1.0 / (args.ny - 1)

    A = build_system_matrix(args.nx, args.ny, args.dt, args.nu, args.vel_x, args.vel_y, comm)
    u = PETSc.Vec().createMPI(n, comm=comm)
    rhs = PETSc.Vec().createMPI(n, comm=comm)
    set_initial_condition(u, args.nx, args.ny, args.ic)

    ksp = PETSc.KSP().create(comm)
    ksp.setOperators(A)
    ksp.setType(args.ksp_type)
    ksp.getPC().setType(args.pc_type)
    ksp.setFromOptions()

    pvd_entries: list[tuple[float, str]] = []
    history_rows: list[tuple[float, int, float]] = []

    def save_snapshot(step: int, time: float) -> None:
        grid = gather_vec_to_grid(u, args.nx, args.ny)
        if rank != 0:
            return
        fname = f"{args.series_name}_{step:04d}.vtp"
        ny, nx = grid.shape
        vel = None
        if not args.no_velocity_vti:
            vel = np.empty((ny, nx, 3), dtype=np.float64)
            vel[..., 0] = args.vel_x
            vel[..., 1] = args.vel_y
            vel[..., 2] = 0.0
        write_grid_surface_vtp(
            outdir / fname,
            grid,
            dx,
            dy,
            "concentration",
            vector_field=vel,
            vector_name="velocity",
        )
        pvd_entries.append((time, fname))

    save_snapshot(0, 0.0)

    for step in range(1, args.steps + 1):
        rhs.copy(u)
        enforce_dirichlet(rhs, args.nx, args.ny)
        ksp.solve(rhs, u)
        history_rows.append((step * args.dt, ksp.getIterationNumber(), float(ksp.getResidualNorm())))
        if step % args.save_every == 0 or step == args.steps:
            save_snapshot(step, step * args.dt)

    if rank == 0:
        pvd_path = outdir / f"{args.series_name}.pvd"
        hist_path = outdir / f"{args.series_name}_solver_history.csv"
        write_pvd(pvd_path, pvd_entries)
        with hist_path.open("w", encoding="utf-8") as fh:
            fh.write("time,iterations,residual_norm\n")
            for t, itn, rn in history_rows:
                fh.write(f"{t:.8e},{itn},{rn:.12e}\n")
        print(f"Wrote {pvd_path}")
        print(f"Wrote {len(pvd_entries)} VTP surface snapshots")
        print(f"Wrote {hist_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
