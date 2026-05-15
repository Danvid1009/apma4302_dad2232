#!/usr/bin/env python3
"""
Weak scaling: wall-clock should stay ~flat as total work grows with p.

Reads ``run_scaling.py`` CSV with ``mode=weak``: ranks, elapsed_s, total_paths.
"""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Weak scaling wall time + total paths.")
    p.add_argument("--csv", type=str, default="output/weak_scaling.csv")
    p.add_argument("--out", type=str, default="visuals/hpc_weak_scaling_walltime.png")
    return p.parse_args()


def read_weak(path: pathlib.Path):
    ranks, el, tp = [], [], []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if row.get("mode", "").lower() != "weak":
                continue
            ranks.append(int(row["ranks"]))
            el.append(float(row["elapsed_s"]))
            tp.append(int(float(row["total_paths"])))
    order = np.argsort(ranks)
    return np.array(ranks, dtype=int)[order], np.array(el)[order], np.array(tp)[order]


def style():
    plt.rcParams.update(
        {
            "figure.facecolor": "#0f1115",
            "axes.facecolor": "#161b22",
            "axes.edgecolor": "#30363d",
            "axes.labelcolor": "#c9d1d9",
            "text.color": "#c9d1d9",
            "xtick.color": "#8b949e",
            "ytick.color": "#8b949e",
            "grid.color": "#30363d",
            "grid.alpha": 0.55,
        }
    )


def main():
    args = parse_args()
    style()
    path = pathlib.Path(args.csv)
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    if not path.exists():
        print(f"missing {path}; run run_scaling.py --mode weak first")
        return 1

    r, elapsed, total = read_weak(path)
    if r.size == 0:
        print("no weak rows in csv")
        return 1

    fig, ax0 = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    fig.suptitle(
        "Weak scaling — wall time vs ranks (flat ≈ perfect weak scaling)",
        color="#f0f6fc",
        fontsize=14,
        fontweight="bold",
    )
    ax0.set_title(
        f"Total paths {int(total[0])} → {int(total[-1])} (paths/rank held ~constant)",
        color="#8b949e",
        fontsize=10,
        pad=10,
    )

    ax0.plot(r, elapsed, "o-", color="#3fb950", lw=2.4, markersize=10, label="Elapsed (s)")
    if elapsed.size:
        ref = float(elapsed[0])
        ax0.axhline(ref, color="#8b949e", ls="--", lw=1.8, label=f"p={int(r[0])} reference ({ref:.3f}s)")
    ax0.set_xlabel("MPI ranks")
    ax0.set_ylabel("Wall time (s)")
    ax0.set_xticks(r)
    ax0.legend(loc="best", framealpha=0.15)
    ax0.grid(True)

    fig.savefig(out, dpi=220, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
