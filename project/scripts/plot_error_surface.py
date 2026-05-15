"""Heatmap of abs_err vs (sigma, n_paths) from run_error_surface.py CSV."""

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Plot error surface CSV as heatmap.")
    p.add_argument("--infile", type=str, default="output/error_surface_european.csv")
    p.add_argument("--out", type=str, default="output/error_surface_abs_err.png")
    return p.parse_args()


def main():
    args = parse_args()
    in_path = pathlib.Path(args.infile)
    rows = list(csv.DictReader(in_path.open("r", newline="")))
    sigmas = sorted({float(r["sigma"]) for r in rows})
    ns = sorted({int(float(r["n_paths"])) for r in rows})

    mat = np.full((len(sigmas), len(ns)), np.nan, dtype=float)
    lookup_s = {s: i for i, s in enumerate(sigmas)}
    lookup_n = {n: j for j, n in enumerate(ns)}
    for r in rows:
        i = lookup_s[float(r["sigma"])]
        j = lookup_n[int(float(r["n_paths"]))]
        mat[i, j] = float(r["abs_err"])

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure()
    im = plt.imshow(
        mat,
        aspect="auto",
        origin="lower",
        extent=[0, len(ns) - 1, 0, len(sigmas) - 1],
    )
    plt.colorbar(im, label="|MC - BS|")
    plt.xlabel("n_paths index (see xticks)")
    plt.ylabel("sigma index (see yticks)")
    plt.xticks(range(len(ns)), [str(n) for n in ns], rotation=45, ha="right")
    plt.yticks(range(len(sigmas)), [f"{s:g}" for s in sigmas])
    plt.title("European MC absolute error vs (N, sigma)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
