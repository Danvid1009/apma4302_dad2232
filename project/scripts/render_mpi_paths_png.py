#!/usr/bin/env python3
"""Matplotlib stand-in for ParaView Tube + color-by-mpi_rank (no VTK install required).

Reads ``export_mpi_domain_paths_vtp.py`` PolyData XML and draws thick 3D polylines
with a cool–warm map on ``mpi_rank``. Intended for ``visuals/paraview_mpi_*.png``.
"""

from __future__ import annotations

import argparse
import pathlib
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  # registers 3d projection


def _named_data_array(parent: ET.Element, name: str) -> ET.Element:
    for c in list(parent):
        if c.tag.endswith("DataArray") and c.get("Name") == name:
            return c
    raise ValueError(f"DataArray Name={name!r} not found under {parent.tag}")


def _points_data_array(piece: ET.Element) -> ET.Element:
    pts_el = piece.find("Points")
    if pts_el is None:
        raise ValueError("no Points")
    for c in list(pts_el):
        if c.tag.endswith("DataArray"):
            return c
    raise ValueError("no Points/DataArray")


def _array_text(elem: ET.Element) -> str:
    t = elem.text
    if not t or not t.strip():
        raise ValueError("empty DataArray text")
    return t.strip()


def load_mpi_vtp(path: pathlib.Path) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    tree = ET.parse(path)
    root = tree.getroot()
    piece = root.find(".//Piece")
    if piece is None:
        raise ValueError("no PolyData Piece in VTP")

    pts_el = _points_data_array(piece)
    ncomp = int(pts_el.get("NumberOfComponents", "3"))
    pts_flat = np.fromstring(_array_text(pts_el), sep=" ", dtype=np.float64)
    pts = pts_flat.reshape(-1, ncomp)
    if pts.shape[1] < 3:
        p2 = np.zeros((pts.shape[0], 3), dtype=np.float64)
        p2[:, : pts.shape[1]] = pts
        pts = p2

    pd = piece.find("PointData")
    if pd is None:
        raise ValueError("no PointData")
    mpi_rank = np.fromstring(_array_text(_named_data_array(pd, "mpi_rank")), sep=" ", dtype=np.int32)

    lines_el = piece.find("Lines")
    if lines_el is None:
        raise ValueError("no Lines")
    conn = np.fromstring(_array_text(_named_data_array(lines_el, "connectivity")), sep=" ", dtype=np.int32)
    offs = np.fromstring(_array_text(_named_data_array(lines_el, "offsets")), sep=" ", dtype=np.int32)

    lines: list[np.ndarray] = []
    lo = 0
    for hi in offs:
        idx = conn[lo:hi]
        lines.append(pts[idx])
        lo = int(hi)
    return pts, mpi_rank, lines


def render(path: pathlib.Path, out: pathlib.Path, *, lw: float, dpi: int) -> None:
    _, mpi_rank, lines = load_mpi_vtp(path)
    rmin = int(mpi_rank.min())
    rmax = int(mpi_rank.max())
    norm = plt.Normalize(rmin, max(rmax, rmin + 1))
    cmap = cm.coolwarm

    fig = plt.figure(figsize=(10.0, 7.0), facecolor="#2e2e2e")
    ax = fig.add_subplot(111, projection="3d", facecolor="#2e2e2e")
    ax.set_facecolor("#2e2e2e")
    for poly in lines:
        rid = int(mpi_rank[int(poly.shape[0] // 2)])
        color = cmap(norm(rid))
        ax.plot(poly[:, 0], poly[:, 1], poly[:, 2], color=color, linewidth=lw, alpha=0.92)

    ax.set_xlabel("time", color="0.85")
    ax.set_ylabel("price", color="0.85")
    ax.set_zlabel("z (rank layout)", color="0.85")
    ax.tick_params(colors="0.75")
    ax.xaxis.pane.fill = False  # type: ignore[attr-defined]
    ax.yaxis.pane.fill = False  # type: ignore[attr-defined]
    ax.zaxis.pane.fill = False  # type: ignore[attr-defined]
    ax.grid(False)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.04)
    cb.set_label("mpi_rank", color="0.9")
    cb.ax.yaxis.set_tick_params(color="0.85")
    plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color="0.85")

    ax.view_init(elev=22, azim=-58)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="#2e2e2e")
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    p = argparse.ArgumentParser(description="Render mpi_domain_paths.vtp to PNG (matplotlib).")
    p.add_argument("--vtp", type=pathlib.Path, required=True)
    p.add_argument("--out", type=pathlib.Path, required=True)
    p.add_argument("--linewidth", type=float, default=2.4, help="Polyline width (Tube-ish).")
    p.add_argument("--dpi", type=int, default=200)
    args = p.parse_args()
    render(args.vtp, args.out, lw=args.linewidth, dpi=args.dpi)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
