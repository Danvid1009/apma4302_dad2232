#!/usr/bin/env python3
"""Create phase 3 convergence and scaling figures."""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot convergence and scaling figures.")
    parser.add_argument("--conv", default="output/poisson_convergence.csv")
    parser.add_argument("--scaling", default="output/strong_scaling.csv")
    parser.add_argument("--scaling-repeats", default="")
    parser.add_argument("--outdir", default="visuals")
    parser.add_argument("--summary-out", default="output/strong_scaling_summary.csv")
    args = parser.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    conv_path = pathlib.Path(args.conv)
    if conv_path.exists():
        conv = np.genfromtxt(conv_path, delimiter=",", names=True)
        h = 1.0 / (conv["n"] - 1.0)
        plt.figure(figsize=(7.5, 5.0))
        plt.loglog(h, conv["error_l2"], "o-", linewidth=2.0, label="L2 error")
        if len(h) > 1:
            coeff = np.polyfit(np.log(h), np.log(conv["error_l2"]), 1)
            slope = coeff[0]
            ref = np.exp(coeff[1]) * h ** slope
            plt.loglog(h, ref, "--", label=f"fit slope={slope:.2f}")
        plt.xlabel("Grid spacing h")
        plt.ylabel("Error")
        plt.title("Poisson Spatial Convergence")
        plt.grid(alpha=0.3, which="both")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "spatial_convergence_loglog.png", dpi=220)
        plt.close()
        print(f"Wrote {outdir / 'spatial_convergence_loglog.png'}")

    scaling_path = pathlib.Path(args.scaling)
    if scaling_path.exists():
        sc = np.genfromtxt(scaling_path, delimiter=",", names=True)
        if getattr(sc, "shape", ()) == ():
            sc = np.array([sc], dtype=sc.dtype)
        order = np.argsort(sc["ranks"])
        ranks = sc["ranks"][order]
        times = sc["solve_time_s"][order]
        t1 = times[0]
        speedup = t1 / times
        efficiency = speedup / ranks

        plt.figure(figsize=(7.5, 5.0))
        plt.plot(ranks, speedup, "o-", linewidth=2.0, label="Measured")
        if len(ranks) > 1:
            plt.plot(ranks, ranks, "--", label="Ideal")
        plt.xlabel("MPI ranks")
        plt.ylabel("Speedup")
        plt.title("Strong Scaling Speedup")
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "strong_scaling_speedup.png", dpi=220)
        plt.close()
        print(f"Wrote {outdir / 'strong_scaling_speedup.png'}")

        plt.figure(figsize=(7.5, 5.0))
        plt.plot(ranks, efficiency, "o-", linewidth=2.0)
        plt.xlabel("MPI ranks")
        plt.ylabel("Parallel efficiency")
        plt.title("Strong Scaling Efficiency")
        plt.ylim(0.0, 1.1)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(outdir / "strong_scaling_efficiency.png", dpi=220)
        plt.close()
        print(f"Wrote {outdir / 'strong_scaling_efficiency.png'}")

    repeats_path = pathlib.Path(args.scaling_repeats) if args.scaling_repeats else pathlib.Path()
    if args.scaling_repeats and repeats_path.exists():
        grouped: dict[int, list[float]] = {}
        with repeats_path.open("r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                rank = int(row["ranks"])
                t = float(row["solve_time_s"])
                grouped.setdefault(rank, []).append(t)

        ranks = np.array(sorted(grouped.keys()), dtype=int)
        means = np.array([np.mean(grouped[r]) for r in ranks], dtype=float)
        stds = np.array([np.std(grouped[r], ddof=1) if len(grouped[r]) > 1 else 0.0 for r in ranks], dtype=float)
        counts = np.array([len(grouped[r]) for r in ranks], dtype=int)
        base_idx = int(np.where(ranks == 1)[0][0]) if 1 in ranks else 0
        base_mean = means[base_idx]
        speedup = base_mean / means
        efficiency = speedup / ranks
        speedup_std = speedup * np.sqrt((stds[base_idx] / max(base_mean, 1e-16)) ** 2 + (stds / np.maximum(means, 1e-16)) ** 2)

        # Save summary CSV for report tables.
        summary_path = pathlib.Path(args.summary_out)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with summary_path.open("w", encoding="utf-8") as fh:
            fh.write("ranks,count,mean_time_s,std_time_s,speedup,efficiency\n")
            for i, r in enumerate(ranks):
                fh.write(
                    f"{r},{counts[i]},{means[i]:.12e},{stds[i]:.12e},{speedup[i]:.8f},{efficiency[i]:.8f}\n"
                )
        print(f"Wrote {summary_path}")

        plt.figure(figsize=(7.5, 5.0))
        plt.errorbar(ranks, speedup, yerr=speedup_std, fmt="o-", linewidth=2.0, capsize=4, label="Mean +/- 1 std")
        plt.plot(ranks, ranks, "--", label="Ideal")
        plt.xlabel("MPI ranks")
        plt.ylabel("Speedup")
        plt.title("Strong Scaling Speedup (Repeated Runs)")
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "strong_scaling_speedup_errorbars.png", dpi=220)
        plt.close()
        print(f"Wrote {outdir / 'strong_scaling_speedup_errorbars.png'}")

        plt.figure(figsize=(7.5, 5.0))
        plt.errorbar(ranks, means, yerr=stds, fmt="o-", linewidth=2.0, capsize=4)
        plt.xlabel("MPI ranks")
        plt.ylabel("Solve time (s)")
        plt.title("Strong Scaling Runtime (Mean +/- 1 std)")
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(outdir / "strong_scaling_runtime_errorbars.png", dpi=220)
        plt.close()
        print(f"Wrote {outdir / 'strong_scaling_runtime_errorbars.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
