#!/usr/bin/env python3
"""
Strong scaling: actual speedup vs ideal (y = p), plus parallel efficiency.

Reads CSV from ``run_scaling.py`` (``mode=strong``): ranks, speedup, efficiency.
"""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Speedup vs ideal + efficiency (strong scaling).")
    p.add_argument("--csv", type=str, default="output/strong_scaling.csv")
    p.add_argument("--out", type=str, default="visuals/hpc_speedup_vs_ideal.png")
    return p.parse_args()


def read_strong(path: pathlib.Path):
    ranks, sp, ef = [], [], []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if row.get("mode", "").lower() != "strong":
                continue
            ranks.append(int(row["ranks"]))
            sp.append(float(row["speedup"]))
            ef.append(float(row["efficiency"]))
    order = np.argsort(ranks)
    return np.array(ranks, dtype=int)[order], np.array(sp)[order], np.array(ef)[order]


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
        print(f"missing {path}; run run_scaling.py --mode strong first")
        return 1

    r, su, ef = read_strong(path)
    if r.size == 0:
        print("no strong rows in csv")
        return 1

    ideal = r.astype(np.float64) / float(r[0])

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12.2, 5.2), constrained_layout=True)
    fig.suptitle("Strong scaling — speedup vs ideal & parallel efficiency", color="#f0f6fc", fontsize=14, fontweight="bold")

    ax0.plot(r, ideal, "--", color="#8b949e", lw=2.0, label="Ideal (linear speedup)")
    ax0.plot(r, su, "o-", color="#58a6ff", lw=2.2, markersize=9, label="Measured speedup")
    ax0.set_xlabel("MPI ranks")
    ax0.set_ylabel("Speedup (vs 1 rank)")
    ax0.set_xticks(r)
    ax0.legend(loc="upper left", framealpha=0.15)
    ax0.grid(True)

    ax1.axhline(1.0, color="#8b949e", ls="--", lw=1.5, label="100% efficiency")
    ax1.plot(r, ef, "s-", color="#f0883e", lw=2.2, markersize=8, label="Parallel efficiency")
    ax1.set_xlabel("MPI ranks")
    ax1.set_ylabel("Efficiency = speedup / p")
    ax1.set_xticks(r)
    ax1.set_ylim(0.0, min(1.15, float(np.max(ef)) * 1.15 + 0.05))
    ax1.legend(loc="upper right", framealpha=0.15)
    ax1.grid(True)

    fig.savefig(out, dpi=220, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
