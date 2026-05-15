import argparse
import csv
import pathlib
import sys
import xml.etree.ElementTree as ET

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.gbm import simulate_gbm_paths

# Default polyline exports share one directory with other ParaView geometry.
DEFAULT_MC_VTP_DIR = ROOT / "output" / "paraview" / "vtp"


def parse_args():
    parser = argparse.ArgumentParser(description="Export Monte Carlo paths for ParaView.")
    parser.add_argument("--s0", type=float, default=100.0)
    parser.add_argument("--r", type=float, default=0.05)
    parser.add_argument("--sigma", type=float, default=0.2)
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=252)
    parser.add_argument("--paths", type=int, default=200)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--barrier", type=float, default=130.0)
    parser.add_argument("--out", type=str, default="output/path_bundle.csv")
    parser.add_argument("--barrier-out", type=str, default="output/barrier_plane.csv")
    parser.add_argument(
        "--vtp",
        nargs="?",
        const="AUTO",
        default=None,
        metavar="PATH",
        help="Also write polyline .vtp for ParaView (no Table To Points). "
        "Use flag alone → output/paraview/vtp/path_bundle.vtp (barrier_plane.vtp alongside); "
        "or pass an explicit .vtp path.",
    )
    parser.add_argument(
        "--vtp-z-scale",
        type=float,
        default=0.02,
        help="Only for --vtp-axes physical: Z offset per path_id (default: 0.02).",
    )
    parser.add_argument(
        "--vtp-axes",
        choices=("normalized", "physical"),
        default="normalized",
        help="VTP vertex layout: normalized unit cube [0,1]^3 (default, even ParaView framing) "
        "or physical (time, price, path_id×z-scale).",
    )
    return parser.parse_args()


def export_paths(out_path: pathlib.Path, paths: np.ndarray, t: float) -> None:
    n_paths, n_cols = paths.shape
    times = np.linspace(0.0, t, n_cols)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["path_id", "time", "stock_price"])
        for i in range(n_paths):
            for j, tt in enumerate(times):
                writer.writerow([i, float(tt), float(paths[i, j])])


def export_barrier_plane(out_path: pathlib.Path, barrier: float, t: float, n_paths: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["path_id", "time", "barrier"])
        for i in range(n_paths):
            writer.writerow([i, 0.0, barrier])
            writer.writerow([i, t, barrier])


def write_path_bundle_vtp(
    out_path: pathlib.Path,
    times: np.ndarray,
    paths: np.ndarray,
    *,
    axes: str,
    z_path_scale: float,
) -> None:
    """VTK PolyData with one vtkLine per path.

    ``axes=='normalized'`` (default): map (time, price, path index) into ~[0,1]^3 so ParaView
    does not squash the bundle. Point arrays ``time``, ``stock_price``, ``path_id`` keep
    physical values for coloring / tooltips.

    ``axes=='physical'``: legacy layout X=time, Y=price, Z=path_id * z_path_scale.
    """
    n_paths, n_cols = paths.shape
    if times.size != n_cols:
        raise ValueError(f"times.size {times.size} != paths.shape[1] {n_cols}")
    t_span = float(times[-1] - times[0]) + 1e-15
    s_min = float(paths.min())
    s_max = float(paths.max())
    s_span = s_max - s_min + 1e-15
    z_den = max(n_paths - 1, 1)

    n_points = n_paths * n_cols
    points = np.zeros((n_points, 3), dtype=np.float64)
    path_id = np.zeros(n_points, dtype=np.int32)
    stock_price = np.zeros(n_points, dtype=np.float64)
    time_arr = np.zeros(n_points, dtype=np.float64)
    connectivity: list[int] = []
    offsets: list[int] = []
    for i in range(n_paths):
        if axes == "normalized":
            z0 = float(i) / z_den
        else:
            z0 = float(i) * z_path_scale
        for j in range(n_cols):
            idx = i * n_cols + j
            tt = float(times[j])
            sp = float(paths[i, j])
            time_arr[idx] = tt
            stock_price[idx] = sp
            path_id[idx] = i
            if axes == "normalized":
                points[idx, 0] = (tt - float(times[0])) / t_span
                points[idx, 1] = (sp - s_min) / s_span
                points[idx, 2] = z0
            else:
                points[idx, 0] = tt
                points[idx, 1] = sp
                points[idx, 2] = z0
            connectivity.append(idx)
        offsets.append((i + 1) * n_cols)

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
    point_data = ET.SubElement(piece, "PointData", Scalars="stock_price")
    arr_sp = ET.SubElement(point_data, "DataArray", type="Float64", Name="stock_price", format="ascii")
    arr_sp.text = " ".join(f"{v:.12e}" for v in stock_price)
    arr_pid = ET.SubElement(point_data, "DataArray", type="Int32", Name="path_id", format="ascii")
    arr_pid.text = " ".join(map(str, path_id.tolist()))
    arr_time = ET.SubElement(point_data, "DataArray", type="Float64", Name="time", format="ascii")
    arr_time.text = " ".join(f"{v:.12e}" for v in time_arr)
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


def write_barrier_lines_vtp(
    out_path: pathlib.Path,
    barrier: float,
    t0: float,
    t1: float,
    n_paths: int,
    s_min: float,
    s_max: float,
    *,
    axes: str,
    z_path_scale: float,
) -> None:
    """One 2-point line per path at constant barrier (t=0 and t=T), same frame as ``write_path_bundle_vtp``."""
    t_span = float(t1 - t0) + 1e-15
    s_span = float(s_max - s_min) + 1e-15
    z_den = max(n_paths - 1, 1)
    n_points = 2 * n_paths
    points = np.zeros((n_points, 3), dtype=np.float64)
    path_id = np.zeros(n_points, dtype=np.int32)
    barrier_arr = np.zeros(n_points, dtype=np.float64)
    time_arr = np.zeros(n_points, dtype=np.float64)
    connectivity: list[int] = []
    offsets: list[int] = []
    for i in range(n_paths):
        z0 = float(i) / z_den if axes == "normalized" else float(i) * z_path_scale
        base = 2 * i
        for k, tt in enumerate((t0, t1)):
            idx = base + k
            time_arr[idx] = float(tt)
            path_id[idx] = i
            barrier_arr[idx] = float(barrier)
            if axes == "normalized":
                points[idx, 0] = (float(tt) - t0) / t_span
                points[idx, 1] = (float(barrier) - s_min) / s_span
                points[idx, 2] = z0
            else:
                points[idx, 0] = float(tt)
                points[idx, 1] = float(barrier)
                points[idx, 2] = z0
            connectivity.append(idx)
        offsets.append(2 * (i + 1))

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
    point_data = ET.SubElement(piece, "PointData", Scalars="barrier")
    arr_b = ET.SubElement(point_data, "DataArray", type="Float64", Name="barrier", format="ascii")
    arr_b.text = " ".join(f"{v:.12e}" for v in barrier_arr)
    arr_pid = ET.SubElement(point_data, "DataArray", type="Int32", Name="path_id", format="ascii")
    arr_pid.text = " ".join(map(str, path_id.tolist()))
    arr_time = ET.SubElement(point_data, "DataArray", type="Float64", Name="time", format="ascii")
    arr_time.text = " ".join(f"{v:.12e}" for v in time_arr)
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
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    paths = simulate_gbm_paths(args.s0, args.r, args.sigma, args.t, args.steps, args.paths, rng)

    export_paths(pathlib.Path(args.out), paths, args.t)
    export_barrier_plane(pathlib.Path(args.barrier_out), args.barrier, args.t, args.paths)

    print(f"wrote paths to {args.out}")
    print(f"wrote barrier plane helper to {args.barrier_out}")

    if args.vtp is not None:
        if args.vtp == "AUTO":
            DEFAULT_MC_VTP_DIR.mkdir(parents=True, exist_ok=True)
            vtp_path = DEFAULT_MC_VTP_DIR / "path_bundle.vtp"
        else:
            vtp_path = pathlib.Path(args.vtp)
        times = np.linspace(0.0, args.t, paths.shape[1])
        s_min = float(paths.min())
        s_max = float(paths.max())
        write_path_bundle_vtp(
            vtp_path,
            times,
            paths,
            axes=args.vtp_axes,
            z_path_scale=args.vtp_z_scale,
        )
        barrier_vtp = vtp_path.parent / "barrier_plane.vtp"
        write_barrier_lines_vtp(
            barrier_vtp,
            args.barrier,
            float(times[0]),
            float(times[-1]),
            paths.shape[0],
            s_min,
            s_max,
            axes=args.vtp_axes,
            z_path_scale=args.vtp_z_scale,
        )
        print(f"wrote barrier helper polylines to {barrier_vtp} (ParaView: open this first, then path_bundle.vtp)")
        print(
            f"wrote polyline bundle to {vtp_path} (axes={args.vtp_axes}; "
            "ParaView: Line / Tube; color by stock_price or time)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
