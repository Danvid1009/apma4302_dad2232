#!/usr/bin/env python3
"""Export GBM stock paths as VTP polylines for ParaView Tube rendering."""

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
    """Return time grid and simulated GBM paths with shape (n_paths, n_steps + 1)."""
    rng = np.random.default_rng(seed)
    dt = t_final / n_steps
    t = np.linspace(0.0, t_final, n_steps + 1)

    z = rng.standard_normal((n_paths, n_steps))
    incr = np.exp((r - 0.5 * sigma * sigma) * dt + sigma * np.sqrt(dt) * z)
    s = np.empty((n_paths, n_steps + 1), dtype=np.float64)
    s[:, 0] = s0
    s[:, 1:] = s0 * np.cumprod(incr, axis=1)
    return t, s


def write_vtp_polylines(
    out_path: pathlib.Path,
    t: np.ndarray,
    s: np.ndarray,
    z_scale: float,
    barrier: float,
    layout: str,
) -> None:
    """Write path bundle as VTK PolyData with line connectivity."""
    n_paths, n_cols = s.shape
    n_points = n_paths * n_cols

    points = np.zeros((n_points, 3), dtype=np.float64)
    path_id = np.zeros(n_points, dtype=np.int32)
    knocked = np.zeros(n_points, dtype=np.int32)
    terminal_price = np.zeros(n_points, dtype=np.float64)
    max_drawdown = np.zeros(n_points, dtype=np.float64)
    barrier_gap = np.zeros(n_points, dtype=np.float64)
    feature_score = np.zeros(n_points, dtype=np.float64)
    connectivity: list[int] = []
    offsets: list[int] = []
    offset = 0

    # Per-path features for a nonlinear embedding (novel geometry vs plain path index).
    s_terminal = s[:, -1]
    running_max = np.maximum.accumulate(s, axis=1)
    drawdown = np.max((running_max - s) / np.maximum(running_max, 1e-12), axis=1)
    max_price = np.max(s, axis=1)
    gap = np.maximum(0.0, barrier - max_price)

    # Normalize robustly into [0, 1] for feature combination.
    def norm(v: np.ndarray) -> np.ndarray:
        lo = float(np.min(v))
        hi = float(np.max(v))
        if hi - lo < 1e-12:
            return np.zeros_like(v)
        return (v - lo) / (hi - lo)

    s_term_n = norm(s_terminal)
    drawdown_n = norm(drawdown)
    gap_n = norm(gap)
    score = 0.55 * s_term_n + 0.30 * (1.0 - drawdown_n) + 0.15 * (1.0 - gap_n)

    for pid in range(n_paths):
        base = pid * n_cols
        path_vals = s[pid]
        crossed = int(np.any(path_vals >= barrier))
        z_base_index = float(pid) * z_scale
        z_base_feature = float(score[pid]) * max(z_scale * max(n_paths - 1, 1), 1e-6)
        for j in range(n_cols):
            idx = base + j
            points[idx, 0] = float(t[j])
            points[idx, 1] = float(path_vals[j])
            if layout == "feature":
                # Add a light time-varying wobble so lines separate without looking artificial.
                wobble = 0.08 * z_scale * np.sin(6.0 * np.pi * t[j] + 0.2 * pid)
                points[idx, 2] = z_base_feature + wobble
            else:
                points[idx, 2] = z_base_index
            path_id[idx] = pid
            knocked[idx] = crossed
            terminal_price[idx] = s_terminal[pid]
            max_drawdown[idx] = drawdown[pid]
            barrier_gap[idx] = gap[pid]
            feature_score[idx] = score[pid]
            connectivity.append(idx)
        offset += n_cols
        offsets.append(offset)

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

    point_data = ET.SubElement(piece, "PointData", Scalars="path_id")
    arr_pid = ET.SubElement(point_data, "DataArray", type="Int32", Name="path_id", format="ascii")
    arr_pid.text = " ".join(map(str, path_id.tolist()))

    arr_knock = ET.SubElement(point_data, "DataArray", type="Int32", Name="knocked_out", format="ascii")
    arr_knock.text = " ".join(map(str, knocked.tolist()))

    arr_term = ET.SubElement(point_data, "DataArray", type="Float64", Name="terminal_price", format="ascii")
    arr_term.text = " ".join(f"{v:.12e}" for v in terminal_price.tolist())

    arr_dd = ET.SubElement(point_data, "DataArray", type="Float64", Name="max_drawdown", format="ascii")
    arr_dd.text = " ".join(f"{v:.12e}" for v in max_drawdown.tolist())

    arr_gap = ET.SubElement(point_data, "DataArray", type="Float64", Name="barrier_gap", format="ascii")
    arr_gap.text = " ".join(f"{v:.12e}" for v in barrier_gap.tolist())

    arr_score = ET.SubElement(point_data, "DataArray", type="Float64", Name="feature_score", format="ascii")
    arr_score.text = " ".join(f"{v:.12e}" for v in feature_score.tolist())

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


def main() -> int:
    parser = argparse.ArgumentParser(description="Export GBM paths as VTP polylines.")
    parser.add_argument("--paths", type=int, default=200)
    parser.add_argument("--steps", type=int, default=252)
    parser.add_argument("--s0", type=float, default=100.0)
    parser.add_argument("--strike", type=float, default=100.0)
    parser.add_argument("--barrier", type=float, default=130.0)
    parser.add_argument("--r", type=float, default=0.05)
    parser.add_argument("--sigma", type=float, default=0.2)
    parser.add_argument("--t-final", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--z-scale", type=float, default=0.02)
    parser.add_argument(
        "--layout",
        choices=["index", "feature"],
        default="feature",
        help="Path embedding in z: index (classic) or feature (nonlinear).",
    )
    parser.add_argument("--out", default="output/paraview/vtp/stock_paths.vtp")
    args = parser.parse_args()

    t, s = simulate_gbm_paths(
        n_paths=args.paths,
        n_steps=args.steps,
        s0=args.s0,
        r=args.r,
        sigma=args.sigma,
        t_final=args.t_final,
        seed=args.seed,
    )
    out = pathlib.Path(args.out)
    write_vtp_polylines(out, t, s, args.z_scale, args.barrier, args.layout)
    print(f"Wrote {out}")
    print(
        f"Model: GBM paths={args.paths} steps={args.steps} "
        f"s0={args.s0} r={args.r} sigma={args.sigma} T={args.t_final} layout={args.layout}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
