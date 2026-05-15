#!/usr/bin/env python3
"""Bar charts from ``run_poisson_pc_sweep.py`` CSV (solve time + KSP iterations per PC)."""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Plot Poisson PC sweep CSV.")
    p.add_argument("--infile", type=str, default="output/poisson_pc_sweep.csv")
    p.add_argument("--out-prefix", type=str, default="visuals/poisson_pc_sweep")
    return p.parse_args()


def _truthy(val: str) -> bool:
    return str(val).strip().lower() in ("1", "true", "yes")


def main() -> int:
    args = parse_args()
    in_path = pathlib.Path(args.infile)
    out_prefix = pathlib.Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(in_path.open(newline="")))
    ok_rows = [r for r in rows if _truthy(r.get("ok", "true"))]
    if not ok_rows:
        print(f"No ok rows in {in_path}; not writing figures.")
        return 1

    labels = [r["pc_type"] for r in ok_rows]
    times = np.array([float(r["solve_time_s"]) for r in ok_rows])
    iters = np.array([int(float(r["iterations"])) for r in ok_rows])

    x = np.arange(len(labels))

    plt.figure(figsize=(7.5, 4.5))
    plt.bar(x, times, color="steelblue")
    plt.xticks(x, labels)
    plt.ylabel("Wall time (s)")
    plt.title("Poisson solve: PC type vs time")
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    p_time = out_prefix.with_name(out_prefix.name + "_time.png")
    plt.savefig(p_time, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"wrote {p_time}")

    plt.figure(figsize=(7.5, 4.5))
    plt.bar(x, iters, color="darkorange")
    plt.xticks(x, labels)
    plt.ylabel("KSP iterations")
    plt.title("Poisson solve: PC type vs iterations")
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    p_it = out_prefix.with_name(out_prefix.name + "_iters.png")
    plt.savefig(p_it, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"wrote {p_it}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
