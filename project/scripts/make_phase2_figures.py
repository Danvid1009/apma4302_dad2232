#!/usr/bin/env python3
"""Generate core Phase 2 PNG figures from solver history CSV files."""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def load_history(path: pathlib.Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return data["time"], data["iterations"], data["residual_norm"]


def main() -> int:
    parser = argparse.ArgumentParser(description="Build phase 2 residual/history figures.")
    parser.add_argument("--heat-history", default="output/paraview/vtp/heat_solver_history.csv")
    parser.add_argument(
        "--ad-history",
        default="output/paraview/vtp/advection_diffusion_solver_history.csv",
    )
    parser.add_argument("--outdir", default="visuals")
    args = parser.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ht, hit, hr = load_history(pathlib.Path(args.heat_history))
    at, ait, ar = load_history(pathlib.Path(args.ad_history))
    hr_safe = np.maximum(hr, 1e-16)
    ar_safe = np.maximum(ar, 1e-16)

    plt.figure(figsize=(8.5, 5.0))
    plt.semilogy(ht, hr_safe, label="Heat residual norm", linewidth=2.0)
    plt.semilogy(at, ar_safe, label="Advection-diffusion residual norm", linewidth=2.0)
    plt.xlabel("Time")
    plt.ylabel("Residual norm")
    plt.title("PETSc KSP Residual History")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "poisson_residual_history.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8.5, 5.0))
    plt.plot(ht, hit, label="Heat KSP iterations", linewidth=2.0)
    plt.plot(at, ait, label="Advection-diffusion KSP iterations", linewidth=2.0)
    plt.xlabel("Time")
    plt.ylabel("Iterations per step")
    plt.title("Solver Iterations per Time Step")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "solver_compare_runtime.png", dpi=200)
    plt.close()

    print(f"Wrote {outdir / 'poisson_residual_history.png'}")
    print(f"Wrote {outdir / 'solver_compare_runtime.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
