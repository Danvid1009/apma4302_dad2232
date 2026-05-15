"""
Triptych: MC convergence — error vs N (log–log with O(N^{-1/2}) guide), price vs N, stderr vs N.
"""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="HPC + statistics convergence triptych.")
    p.add_argument("--infile", type=str, default="output/convergence.csv")
    p.add_argument("--out", type=str, default="visuals/hpc_convergence_triptych.png")
    return p.parse_args()


def apply_style():
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
            "grid.alpha": 0.55,
            "font.size": 10,
        }
    )


def main():
    args = parse_args()
    apply_style()
    in_path = pathlib.Path(args.infile)
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    paths, mae, rmse, mc_p, bs_p, stderr = [], [], [], [], [], []
    with in_path.open(newline="") as f:
        for row in csv.DictReader(f):
            paths.append(float(row["paths"]))
            mae.append(float(row["mean_abs_error"]))
            rmse.append(float(row["rmse_error"]))
            mc_p.append(float(row["mc_price"]))
            bs_p.append(float(row["bs_price"]))
            stderr.append(float(row["std_error"]))

    n = np.array(paths)
    mae = np.array(mae)
    rmse = np.array(rmse)
    mc_p = np.array(mc_p)
    bs_p = np.array(bs_p)
    stderr = np.array(stderr)

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), constrained_layout=True)
    fig.suptitle(
        "European MC convergence — statistical error vs path budget $N$",
        color="#f0f6fc",
        fontsize=13,
        fontweight="bold",
    )

    # log-log error with slope -1/2 reference
    ax = axes[0]
    ref = mae[0] * (n[0] / n) ** 0.5
    ax.loglog(n, mae, "o-", color="#58a6ff", lw=2, ms=7, label=r"Mean $|error|$")
    ax.loglog(n, rmse, "s--", color="#a371f7", lw=1.8, ms=6, label="RMSE (reps)")
    ax.loglog(n, ref, ":", color="#3fb950", lw=2, label=r"$O(N^{-1/2})$ guide")
    ax.set_xlabel(r"Paths $N$")
    ax.set_ylabel("Error")
    ax.set_title("Log–log slope check")
    ax.legend(loc="upper right", framealpha=0.15)
    ax.grid(True, which="both")

    ax = axes[1]
    bs0 = float(bs_p[0])
    ax.axhline(bs0, color="#3fb950", lw=2, label="Black–Scholes")
    ax.fill_between(n, bs0 * 0.998, bs0 * 1.002, color="#3fb950", alpha=0.12)
    ax.plot(n, mc_p, "o-", color="#f0883e", lw=2, ms=7, label="MC mean price")
    ax.set_xscale("log")
    ax.set_xlabel(r"Paths $N$")
    ax.set_ylabel("Price")
    ax.set_title("MC vs analytic anchor")
    ax.legend(loc="lower right", framealpha=0.15)
    ax.grid(True, which="both")

    ax = axes[2]
    ax.loglog(n, stderr, "D-", color="#ff7b72", lw=2, ms=7)
    ax.set_xlabel(r"Paths $N$")
    ax.set_ylabel(r"MC std error $\hat\sigma/\sqrt{N}$")
    ax.set_title("Posterior uncertainty vs budget")
    ax.grid(True, which="both")

    fig.savefig(out_path, dpi=200, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
