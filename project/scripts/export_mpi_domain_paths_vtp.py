#!/usr/bin/env python3
"""
Export GBM paths as **VTP polylines** with a synthetic **MPI rank** field for
ParaView coloring — visual “domain decomposition” of path batches (same idea as
``_split_paths`` in ``engine.py``).

ParaView: color by ``mpi_rank``, tube filter, screenshot for the report.

Example::

    python scripts/export_mpi_domain_paths_vtp.py \\
      --mpi-ranks 8 --paths-per-rank 16 --steps 120 \\
      --out output/paraview/vtp/mpi_domain_paths.vtp
"""

from __future__ import annotations

import argparse
import pathlib
import xml.etree.ElementTree as ET

import numpy as np


def simulate_gbm_paths(
    n_paths: int,
    n_steps: int,
    s0: float,
    r: float,
    sigma: float,
    t_final: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    dt = t_final / n_steps
    t = np.linspace(0.0, t_final, n_steps + 1)
    z = rng.standard_normal((n_paths, n_steps))
    incr = np.exp((r - 0.5 * sigma * sigma) * dt + sigma * np.sqrt(dt) * z)
    s = np.empty((n_paths, n_steps + 1), dtype=np.float64)
    s[:, 0] = s0
    s[:, 1:] = s0 * np.cumprod(incr, axis=1)
    return t, s


def write_vtp_mpi_paths(
    out_path: pathlib.Path,
    t: np.ndarray,
    all_paths: list[tuple[int, np.ndarray]],
    z_rank_scale: float,
    z_path_wobble: float,
) -> None:
    """all_paths: list of (mpi_rank, s_local) with s_local shape (n_local, n_steps+1)."""
    n_cols = int(t.size)
    n_paths = sum(s.shape[0] for _, s in all_paths)
    n_points = n_paths * n_cols

    points = np.zeros((n_points, 3), dtype=np.float64)
    path_id = np.zeros(n_points, dtype=np.int32)
    mpi_rank = np.zeros(n_points, dtype=np.int32)
    connectivity: list[int] = []
    offsets: list[int] = []
    offset = 0
    gid = 0

    for rank, s in all_paths:
        n_loc, nc = s.shape
        assert nc == n_cols
        for pid in range(n_loc):
            base = gid * n_cols
            path_vals = s[pid]
            z_base = float(rank) * z_rank_scale + float(pid) * z_path_wobble
            for j in range(n_cols):
                idx = base + j
                points[idx, 0] = float(t[j])
                points[idx, 1] = float(path_vals[j])
                points[idx, 2] = z_base + 0.04 * z_rank_scale * np.sin(4.0 * np.pi * t[j] / max(t[-1], 1e-12))
                path_id[idx] = gid
                mpi_rank[idx] = int(rank)
                connectivity.append(idx)
            offset += n_cols
            offsets.append(offset)
            gid += 1

    root = ET.Element("VTKFile", type="PolyData", version="0.1", byte_order="LittleEndian")
    poly = ET.SubElement(root, "PolyData")
    piece = ET.SubElement(
        poly,
        "Piece",
        NumberOfPoints=str(n_points),
        NumberOfVerts="0",
        NumberOfLines=str(n_paths),
        NumberOfStrips="0",
        NumberOfPolys="0",
    )

    point_data = ET.SubElement(piece, "PointData", Scalars="mpi_rank")
    arr_r = ET.SubElement(point_data, "DataArray", type="Int32", Name="mpi_rank", format="ascii")
    arr_r.text = " ".join(map(str, mpi_rank.tolist()))
    arr_pid = ET.SubElement(point_data, "DataArray", type="Int32", Name="global_path_id", format="ascii")
    arr_pid.text = " ".join(map(str, path_id.tolist()))

    pts = ET.SubElement(piece, "Points")
    arr_pts = ET.SubElement(pts, "DataArray", type="Float64", NumberOfComponents="3", format="ascii")
    arr_pts.text = " ".join(f"{v:.12e}" for v in points.reshape(-1))

    lines = ET.SubElement(piece, "Lines")
    arr_conn = ET.SubElement(lines, "DataArray", type="Int32", Name="connectivity", format="ascii")
    arr_conn.text = " ".join(map(str, connectivity))
    arr_off = ET.SubElement(lines, "DataArray", type="Int32", Name="offsets", format="ascii")
    arr_off.text = " ".join(map(str, offsets))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(root).write(out_path, encoding="utf-8", xml_declaration=True)


def parse_args():
    p = argparse.ArgumentParser(description="VTP path bundle colored by synthetic MPI rank.")
    p.add_argument("--mpi-ranks", type=int, default=8, help="Number of synthetic ranks (path batches).")
    p.add_argument("--paths-per-rank", type=int, default=16)
    p.add_argument("--steps", type=int, default=120)
    p.add_argument("--s0", type=float, default=100.0)
    p.add_argument("--r", type=float, default=0.05)
    p.add_argument("--sigma", type=float, default=0.2)
    p.add_argument("--t-final", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=4242)
    p.add_argument("--z-rank-scale", type=float, default=0.35, help="Separation in z between ranks.")
    p.add_argument("--z-path-scale", type=float, default=0.004, help="Small offset between paths within a rank.")
    p.add_argument("--out", type=str, default="output/paraview/vtp/mpi_domain_paths.vtp")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    P = max(1, int(args.mpi_ranks))
    nloc = max(1, int(args.paths_per_rank))
    t_first, _ = simulate_gbm_paths(1, args.steps, args.s0, args.r, args.sigma, args.t_final, args.seed)
    t = t_first

    chunks: list[tuple[int, np.ndarray]] = []
    for rank in range(P):
        seed_r = int(args.seed) + 100_007 * rank
        _, s = simulate_gbm_paths(nloc, args.steps, args.s0, args.r, args.sigma, args.t_final, seed_r)
        chunks.append((rank, s))

    out = pathlib.Path(args.out)
    write_vtp_mpi_paths(out, t, chunks, args.z_rank_scale, args.z_path_scale)
    n_total = P * nloc
    print(f"Wrote {out} ({P} synthetic ranks × {nloc} paths = {n_total} polylines)")
    print("ParaView: open .vtp → Tube → Color by mpi_rank → Cool to Warm → screenshot.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
