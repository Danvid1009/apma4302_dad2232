"""Helpers for PETSc solution gathering and ParaView VTI/PVD/VTP export."""

from __future__ import annotations

import pathlib
import xml.etree.ElementTree as ET
from typing import Iterable

import numpy as np
from petsc4py import PETSc


def gather_vec_to_grid(vec: PETSc.Vec, nx: int, ny: int) -> np.ndarray:
    """Gather a distributed PETSc Vec onto rank 0 and reshape as (ny, nx)."""
    comm = vec.getComm()
    rank = comm.getRank()
    local = vec.getArray(readonly=True).copy()
    gathered = comm.tompi4py().gather(local, root=0)
    if rank != 0:
        return np.empty((0, 0), dtype=np.float64)
    return np.concatenate(gathered).reshape((ny, nx))


def write_vti(
    path: pathlib.Path,
    field: np.ndarray,
    dx: float,
    dy: float,
    scalar_name: str,
    *,
    vector_field: np.ndarray | None = None,
    vector_name: str = "velocity",
) -> None:
    """Write legacy VTK XML ImageData (.vti) with optional per-point 3-vector array.

    ``vector_field`` must have shape ``(ny, nx, 3)`` matching ``field`` and use the same
    row-major flatten order as ``field`` (last index ``x`` varies fastest). When present,
    ``PointData`` is tagged with ``Vectors`` so ParaView can drive **Glyph** filters.
    """
    ny, nx = field.shape
    if vector_field is not None:
        if vector_field.shape != (ny, nx, 3):
            raise ValueError(f"vector_field must be (ny,nx,3)={(ny, nx, 3)}, got {vector_field.shape}")
    attrs: dict[str, str] = {"Scalars": scalar_name}
    if vector_field is not None:
        attrs["Vectors"] = vector_name
    root = ET.Element("VTKFile", type="ImageData", version="0.1", byte_order="LittleEndian")
    image = ET.SubElement(
        root,
        "ImageData",
        WholeExtent=f"0 {nx - 1} 0 {ny - 1} 0 0",
        Origin="0 0 0",
        Spacing=f"{dx} {dy} 1.0",
    )
    piece = ET.SubElement(image, "Piece", Extent=f"0 {nx - 1} 0 {ny - 1} 0 0")
    point_data = ET.SubElement(piece, "PointData", attrs)
    data = ET.SubElement(point_data, "DataArray", type="Float64", Name=scalar_name, format="ascii")
    data.text = " ".join(f"{val:.12e}" for val in field.reshape(-1))
    if vector_field is not None:
        vflat = vector_field.reshape(-1)
        vdata = ET.SubElement(
            point_data,
            "DataArray",
            type="Float64",
            Name=vector_name,
            NumberOfComponents="3",
            format="ascii",
        )
        vdata.text = " ".join(f"{val:.12e}" for val in vflat)
    ET.SubElement(piece, "CellData")
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def write_grid_surface_vtp(
    path: pathlib.Path,
    field: np.ndarray,
    dx: float,
    dy: float,
    scalar_name: str,
    *,
    vector_field: np.ndarray | None = None,
    vector_name: str = "velocity",
) -> None:
    """Write a 2D structured scalar field as VTK PolyData quad surface in the z=0 plane.

    One quad per grid cell; point ordering matches ``field.reshape(-1)`` (x varies fastest).
    Optional ``vector_field`` shape ``(ny, nx, 3)`` is stored on points for Glyph filters.
    """
    ny, nx = field.shape
    if nx < 2 or ny < 2:
        raise ValueError("write_grid_surface_vtp requires nx>=2 and ny>=2")
    if vector_field is not None and vector_field.shape != (ny, nx, 3):
        raise ValueError(f"vector_field must be (ny,nx,3), got {vector_field.shape}")

    n_pts = nx * ny
    pts = np.zeros((n_pts, 3), dtype=np.float64)
    pid = 0
    for j in range(ny):
        for i in range(nx):
            pts[pid, 0] = float(i) * dx
            pts[pid, 1] = float(j) * dy
            pts[pid, 2] = 0.0
            pid += 1

    n_polys = (nx - 1) * (ny - 1)
    conn: list[int] = []
    offs: list[int] = []
    acc = 0
    for j in range(ny - 1):
        for i in range(nx - 1):
            i0 = j * nx + i
            conn.extend((i0, i0 + 1, i0 + nx + 1, i0 + nx))
            acc += 4
            offs.append(acc)

    attrs: dict[str, str] = {"Scalars": scalar_name}
    if vector_field is not None:
        attrs["Vectors"] = vector_name

    root = ET.Element("VTKFile", type="PolyData", version="0.1", byte_order="LittleEndian")
    poly = ET.SubElement(root, "PolyData")
    piece = ET.SubElement(
        poly,
        "Piece",
        NumberOfPoints=str(n_pts),
        NumberOfVerts="0",
        NumberOfLines="0",
        NumberOfStrips="0",
        NumberOfPolys=str(n_polys),
    )
    point_data = ET.SubElement(piece, "PointData", attrs)
    sflat = field.reshape(-1)
    arr_s = ET.SubElement(point_data, "DataArray", type="Float64", Name=scalar_name, format="ascii")
    arr_s.text = " ".join(f"{v:.12e}" for v in sflat)
    if vector_field is not None:
        vflat = vector_field.reshape(-1)
        arr_v = ET.SubElement(
            point_data,
            "DataArray",
            type="Float64",
            Name=vector_name,
            NumberOfComponents="3",
            format="ascii",
        )
        arr_v.text = " ".join(f"{v:.12e}" for v in vflat)

    pts_el = ET.SubElement(piece, "Points")
    arr_pts = ET.SubElement(pts_el, "DataArray", type="Float64", NumberOfComponents="3", format="ascii")
    arr_pts.text = " ".join(f"{v:.12e}" for v in pts.reshape(-1))

    polys = ET.SubElement(piece, "Polys")
    arr_c = ET.SubElement(polys, "DataArray", type="Int32", Name="connectivity", format="ascii")
    arr_c.text = " ".join(map(str, conn))
    arr_o = ET.SubElement(polys, "DataArray", type="Int32", Name="offsets", format="ascii")
    arr_o.text = " ".join(map(str, offs))

    path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def write_pvd(path: pathlib.Path, entries: Iterable[tuple[float, str]]) -> None:
    root = ET.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
    coll = ET.SubElement(root, "Collection")
    for timestep, filename in entries:
        ET.SubElement(coll, "DataSet", timestep=f"{timestep:.8g}", group="", part="0", file=filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)
