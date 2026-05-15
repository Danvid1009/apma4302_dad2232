#!/usr/bin/env python3
"""Plot KSP iteration counts along pseudo-time from heat / advection-diffusion solver logs."""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Plot diffusion solver history CSV.")
    p.add_argument(
        "--infile",
        type=str,
        required=True,
        help="CSV from run_heat_paraview.py or run_advection_diffusion_paraview.py "
        "(columns: time,iterations,residual_norm).",
    )
    p.add_argument("--out", type=str, default="")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    in_path = pathlib.Path(args.infile)
    stem = in_path.stem
    if stem.endswith("_solver_history"):
        label = stem[: -len("_solver_history")]
    else:
        label = stem
    out_path = (
        pathlib.Path(args.out)
        if str(args.out).strip()
        else pathlib.Path("visuals") / f"{label}_solver_iters.png"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)

    times, iters, res = [], [], []
    with in_path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            times.append(float(row["time"]))
            iters.append(int(float(row["iterations"])))
            res.append(float(row["residual_norm"]))

    t = np.array(times)
    k = np.array(iters, dtype=float)
    r = np.array(res)

    fig, ax1 = plt.subplots(figsize=(8.0, 4.5))
    ax1.plot(t, k, color="steelblue", linewidth=1.4, label="KSP iterations")
    ax1.set_xlabel("time")
    ax1.set_ylabel("iterations", color="steelblue")
    ax1.tick_params(axis="y", labelcolor="steelblue")
    ax1.grid(alpha=0.3)

    ax2 = ax1.twinx()
    ax2.semilogy(t, np.maximum(r, 1e-30), color="darkorange", linewidth=1.0, alpha=0.85, label="|residual|")
    ax2.set_ylabel("residual norm (log)", color="darkorange")
    ax2.tick_params(axis="y", labelcolor="darkorange")

    plt.title(f"Solver history: {label}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
