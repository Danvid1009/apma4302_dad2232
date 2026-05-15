"""Plot roll-forward sigma and MC price from run_phase2_roll_forward.py output."""

import argparse
import csv
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np

_scripts = pathlib.Path(__file__).resolve().parent
if str(_scripts) not in sys.path:
    sys.path.insert(0, str(_scripts))
from plot_archive import archive_prior_png


def parse_args():
    p = argparse.ArgumentParser(description="Plot phase2 roll-forward CSV.")
    p.add_argument("--infile", type=str, default="output/phase2_roll_forward.csv")
    p.add_argument("--out", type=str, default="visuals/phase2_roll_forward_panel.png")
    return p.parse_args()


def main():
    args = parse_args()
    path = pathlib.Path(args.infile)
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    xs, sig, mc, bs, s0 = [], [], [], [], []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            xs.append(row["asof_date"])
            sig.append(float(row["sigma_roll"]))
            mc.append(float(row["mc_price"]))
            s0.append(float(row["s0"]))
            try:
                bs.append(float(row["bs_price"]))
            except (KeyError, ValueError):
                bs.append(float("nan"))

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

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True, constrained_layout=True)
    fig.suptitle("Phase 2 depth — roll-forward calibrated $\\sigma$ and repricing", color="#f0f6fc", fontsize=13, fontweight="bold")

    idx = np.arange(len(sig), dtype=float)

    axes[0].plot(idx, sig, color="#f0883e", marker="o", lw=2.0, markersize=8)
    axes[0].set_ylabel("Rolling $\\sigma$", color="#f0883e")
    axes[0].tick_params(axis="y", colors="#f0883e")
    axes[0].grid(True)
    axes[0].set_title("Trailing realized vol (annualized)")
    for x, v in zip(idx, sig):
        axes[0].annotate(f"{v:.3f}", (x, v), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8, color="#8b949e")

    # Spot on left scale; option prices on twin (otherwise MC is visually crushed by S0).
    ax_spot = axes[1]
    ax_opts = ax_spot.twinx()
    ax_spot.plot(idx, s0, color="#3fb950", marker="D", lw=2.2, markersize=7, label="$S_0$ (spot)")
    ax_spot.set_ylabel("Spot $S_0$ (USD)", color="#3fb950")
    ax_spot.tick_params(axis="y", colors="#3fb950")

    ax_opts.plot(idx, mc, color="#58a6ff", marker="s", lw=2.0, markersize=7, label="MC call")
    bs_arr = np.asarray(bs, dtype=float)
    if np.all(np.isfinite(bs_arr)):
        ax_opts.plot(idx, bs_arr, color="#d2a8ff", marker="^", lw=1.6, markersize=6, ls="--", label="BS call")
    ax_opts.set_ylabel("Option price (USD)", color="#58a6ff")
    ax_opts.tick_params(axis="y", colors="#58a6ff")

    h1, l1 = ax_spot.get_legend_handles_labels()
    h2, l2 = ax_opts.get_legend_handles_labels()
    ax_spot.legend(h1 + h2, l1 + l2, loc="upper left", framealpha=0.25, fontsize=9)
    ax_spot.grid(True)
    ax_spot.set_xticks(idx)
    ax_spot.set_xticklabels(xs, rotation=30, ha="right")
    ax_spot.set_xlabel("Snapshot date")
    ax_spot.set_title("Spot vs repriced ATM call (twin axes — read MC/BS on the right scale)")

    archive_prior_png(out)
    fig.savefig(out, dpi=200, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
