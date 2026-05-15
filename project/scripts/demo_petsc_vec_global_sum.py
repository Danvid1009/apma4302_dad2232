#!/usr/bin/env python3
"""
PETSc ``Vec`` + parallel layout (each rank owns one global row).

Course tie-in: **Parallel Vec** objects (early-semester LA weeks). Each rank
stores one scalar (e.g. local partial sum of payoffs); ``Vec.sum`` matches a
single ``MPI.Allreduce`` of the same values.

Example::

    mpiexec -n 4 python scripts/demo_petsc_vec_global_sum.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from petsc4py import PETSc


def parse_args():
    p = argparse.ArgumentParser(description="Demo PETSc Vec partitioned sum.")
    p.add_argument("--reps", type=int, default=1, help="Repeat sum for timing (microbench).")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    comm = PETSc.COMM_WORLD
    rank = comm.getRank()
    size = comm.getSize()

    # Each rank contributes (rank+1) as a toy local partial sum.
    local_value = float(rank + 1)
    expected = float(sum(range(1, size + 1)))

    v = PETSc.Vec().create(comm=comm)
    v.setSizes((1, size))
    v.setFromOptions()
    rlo, rhi = v.getOwnershipRange()
    with v.localForm() as loc:
        a = loc.asarray()
        if a.size:
            a[:] = local_value

    total = 0.0
    for _ in range(max(1, args.reps)):
        total = v.sum()

    if rank == 0:
        print(f"ranks={size} vec.sum={total:.12f} expected={expected:.12f} ok={abs(total - expected) < 1e-10}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
