"""
Publication-style scaling dashboard: strong + weak speedup & efficiency in one figure.

Reads CSVs produced by `scripts/run_scaling.py` (columns: mode, ranks, elapsed_s, speedup, efficiency, ...).
"""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Multi-panel HPC scaling dashboard.")
    p.add_argument("--strong-csv", type=str, default="output/strong_scaling.csv")
    p.add_argument("--weak-csv", type=str, default="output/weak_scaling.csv")
    p.add_argument("--out", type=str, default="visuals/hpc_scaling_dashboard.png")
    return p.parse_args()


def read_mode(path: pathlib.Path, mode: str):
    ranks, elapsed, speedup, eff, tpaths, ppr = [], [], [], [], [], []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if row.get("mode", "").lower() != mode:
                continue
            ranks.append(int(row["ranks"]))
            elapsed.append(float(row["elapsed_s"]))
            speedup.append(float(row["speedup"]))
            eff.append(float(row["efficiency"]))
            tpaths.append(int(float(row["total_paths"])))
            ppr.append(int(float(row["paths_per_rank"])))
    order = np.argsort(ranks)
    return (
        np.array(ranks, dtype=int)[order],
        np.array(elapsed)[order],
        np.array(speedup)[order],
        np.array(eff)[order],
        np.array(tpaths)[order],
        np.array(ppr)[order],
    )


def apply_style():
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
            "grid.alpha": 0.6,
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
        }
    )


def main():
    args = parse_args()
    apply_style()

    strong_path = pathlib.Path(args.strong_csv)
    weak_path = pathlib.Path(args.weak_csv)
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rs, es, ss, ef, tp, pr = read_mode(strong_path, "strong")
    rw, ew, sw, efw, tpw, prw = read_mode(weak_path, "weak")

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5), constrained_layout=True)
    fig.suptitle("Monte Carlo option pricing — MPI scaling study", color="#f0f6fc", fontsize=15, fontweight="bold")

    # --- Strong speedup ---
    ax = axes[0, 0]
    ax.fill_between(rs, 0, ss, color="#388bfd", alpha=0.25)
    ax.plot(rs, ss, "o-", color="#58a6ff", lw=2.2, ms=8, label="Measured speedup")
    ax.plot(rs, rs, "--", color="#f0883e", lw=1.8, label="Ideal (linear)")
    ax.set_xlabel("MPI ranks $p$")
    ax.set_ylabel("Speedup $S(p)$")
    ax.set_title("Strong scaling — fixed total paths")
    ax.legend(loc="upper left", framealpha=0.2)
    ax.grid(True)
    ax.set_xticks(rs)

    # --- Strong efficiency ---
    ax = axes[0, 1]
    colors = plt.cm.plasma(np.linspace(0.25, 0.85, len(rs)))
    ax.bar(rs.astype(str), 100.0 * ef, color=colors, edgecolor="#30363d", linewidth=0.8)
    ax.axhline(100.0, color="#3fb950", linestyle="--", lw=1.5, label="100% (perfect)")
    ax.set_xlabel("MPI ranks $p$")
    ax.set_ylabel("Parallel efficiency $E(p)$ (%)")
    ax.set_title("Strong scaling — efficiency")
    ax.legend(loc="upper right", framealpha=0.2)
    ax.grid(True, axis="y")
    ax.set_ylim(0, 105)

    # --- Weak runtime growth ---
    ax = axes[1, 0]
    ax.plot(rw, ew, "s-", color="#a371f7", lw=2.2, ms=8, label="Wall time")
    ax.fill_between(rw, ew, alpha=0.2, color="#a371f7")
    ax.set_xlabel("MPI ranks $p$")
    ax.set_ylabel("Elapsed (s)")
    ax.set_title("Weak scaling — fixed paths / rank (runtime drift)")
    ax2 = ax.twinx()
    ax2.plot(rw, tpw, "D--", color="#79c0ff", lw=1.5, ms=6, alpha=0.9, label="Total paths")
    ax2.set_ylabel("Total paths", color="#79c0ff")
    ax2.tick_params(axis="y", colors="#79c0ff")
    ax.grid(True)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", framealpha=0.2)

    # --- Weak efficiency + annotation ---
    ax = axes[1, 1]
    ax.plot(rw, 100.0 * efw, "o-", color="#ff7b72", lw=2.2, ms=8)
    ax.fill_between(rw, 100.0 * efw, alpha=0.2, color="#ff7b72")
    ax.set_xlabel("MPI ranks $p$")
    ax.set_ylabel("Parallel efficiency $E(p)$ (%)")
    ax.set_title("Weak scaling — efficiency collapse (overhead)")
    ax.grid(True)
    ax.set_xticks(rw)
    note = (
        f"Strong: total paths = {tp[0]:,}\n"
        f"Weak: paths/rank = {prw[0]:,}\n"
        "Interpretation: reductions + launcher overhead\n"
        "dominate when batches shrink or ranks grow."
    )
    ax.text(
        0.04,
        0.06,
        note,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="bottom",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#21262d", edgecolor="#30363d", alpha=0.95),
        family="monospace",
        color="#c9d1d9",
    )

    fig.savefig(out_path, dpi=200, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
