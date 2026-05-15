"""
Barrier sweep as an HPC/statistics story: price vs barrier + relative stderr (uncertainty per geometry).
"""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Barrier sweep visualization.")
    p.add_argument("--infile", type=str, default="output/barrier_sweep.csv")
    p.add_argument("--out", type=str, default="visuals/hpc_barrier_geometry_story.png")
    return p.parse_args()


def main():
    args = parse_args()
    path = pathlib.Path(args.infile)
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    b, price, stderr, wall = [], [], [], []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            b.append(float(row["barrier"]))
            price.append(float(row["price"]))
            stderr.append(float(row["std_error"]))
            wall.append(float(row["wall_s"]))

    b = np.array(b)
    price = np.array(price)
    stderr = np.array(stderr)
    wall = np.array(wall)
    rel = stderr / np.maximum(price, 1e-12)

    plt.rcParams.update(
        {
            "figure.facecolor": "#0d1117",
            "axes.facecolor": "#161b22",
            "axes.edgecolor": "#30363d",
            "axes.labelcolor": "#c9d1d9",
            "text.color": "#c9d1d9",
            "xtick.color": "#8b949e",
            "ytick.color": "#8b949e",
            "grid.color": "#30363d",
            "grid.alpha": 0.5,
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5), constrained_layout=True)
    fig.suptitle(
        "Barrier up-and-out call — geometry vs value & MC noise",
        color="#f0f6fc",
        fontsize=13,
        fontweight="bold",
    )

    sc = axes[0].scatter(b, price, c=wall, cmap="viridis", s=120, edgecolors="#30363d", linewidths=0.6, zorder=3)
    axes[0].plot(b, price, color="#58a6ff", lw=1.2, alpha=0.65, zorder=2)
    cbar = fig.colorbar(sc, ax=axes[0], shrink=0.82, pad=0.02)
    cbar.set_label("Wall time (s)", color="#c9d1d9")
    cbar.ax.yaxis.set_tick_params(color="#8b949e")
    axes[0].set_xlabel("Barrier $B$")
    axes[0].set_ylabel("MC price")
    axes[0].set_title("Price path near knock-out region")
    axes[0].grid(True)

    axes[1].fill_between(b, 0, rel, color="#a371f7", alpha=0.25)
    axes[1].plot(b, rel, "o-", color="#d2a8ff", lw=2, ms=7)
    axes[1].set_xlabel("Barrier $B$")
    axes[1].set_ylabel(r"Std error / MC price (relative noise)")
    axes[1].set_title("Relative statistical noise vs geometry")
    axes[1].grid(True)

    fig.savefig(out, dpi=200, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
