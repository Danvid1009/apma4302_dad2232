import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Plot convergence error vs paths on log-log scale.")
    parser.add_argument("--infile", type=str, default="output/convergence.csv")
    parser.add_argument("--out", type=str, default="output/convergence.png")
    return parser.parse_args()


def main():
    args = parse_args()
    in_path = pathlib.Path(args.infile)
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    paths = []
    errors = []
    with in_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            paths.append(float(row["paths"]))
            # Prefer RMSE if available (more stable), fallback to abs_error.
            err = row.get("rmse_error") or row.get("abs_error")
            errors.append(float(err))

    x = np.array(paths, dtype=float)
    y = np.array(errors, dtype=float)
    y_safe = np.maximum(y, 1e-16)

    slope, intercept = np.polyfit(np.log(x), np.log(y_safe), 1)
    y_fit = np.exp(intercept) * x ** slope

    plt.figure()
    plt.loglog(x, y_safe, marker="o", label="Measured abs error")
    plt.loglog(x, y_fit, linestyle="--", label=f"Fit slope={slope:.3f}")
    plt.loglog(x, y_safe[0] * (x / x[0]) ** (-0.5), linestyle=":", label="Reference slope=-0.5")
    plt.xlabel("Number of paths (N)")
    plt.ylabel("Absolute error |V_MC - V_BS|")
    plt.title("European Call Monte Carlo Convergence")
    plt.legend()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"fitted_loglog_slope={slope:.6f}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
