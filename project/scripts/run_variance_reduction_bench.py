"""
Compare plain vs antithetic Monte Carlo on the same nominal path count.

Emphasizes the HPC-relevant quantity: statistical error (stderr, CI width) per
unit of replicated work at fixed N, before layering parallel speedup.
"""

import argparse
import csv
import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.black_scholes import black_scholes_call
from montecarlo.engine import petsc_monte_carlo_price
from montecarlo.payoffs import european_call_payoff


def parse_args():
    parser = argparse.ArgumentParser(description="Bench plain vs antithetic European MC variance.")
    parser.add_argument("--paths", type=int, default=200000)
    parser.add_argument("--reps", type=int, default=20)
    parser.add_argument("--steps", type=int, default=252)
    parser.add_argument("--s0", type=float, default=100.0)
    parser.add_argument("--strike", type=float, default=100.0)
    parser.add_argument("--r", type=float, default=0.05)
    parser.add_argument("--sigma", type=float, default=0.2)
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--out", type=str, default="output/variance_reduction_bench.csv")
    return parser.parse_args()


def main():
    args = parse_args()
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    payoff_builder = lambda strike: lambda paths: european_call_payoff(paths, strike)
    bs = black_scholes_call(args.s0, args.strike, args.r, args.sigma, args.t)

    rows = []
    for label, use_anti in [("plain", False), ("antithetic", True)]:
        stderrs = []
        ci_widths = []
        wall_times = []
        for rep in range(args.reps):
            seed = args.seed + rep * 7919
            t0 = time.perf_counter()
            mc = petsc_monte_carlo_price(
                args.s0,
                args.strike,
                args.r,
                args.sigma,
                args.t,
                args.steps,
                args.paths,
                payoff_builder=payoff_builder,
                seed=seed,
                antithetic=use_anti,
            )
            wall_times.append(time.perf_counter() - t0)
            if mc["rank"] != 0:
                return
            stderrs.append(mc["std_error"])
            ci_widths.append(mc["ci_high"] - mc["ci_low"])

        rows.append(
            {
                "method": label,
                "mean_stderr": float(np.mean(stderrs)),
                "mean_ci_width": float(np.mean(ci_widths)),
                "mean_wall_s": float(np.mean(wall_times)),
                "paths_effective": mc["n_paths"],
                "ranks": mc["size"],
                "backend": mc["backend"],
                "bs_price": bs,
            }
        )
        print(
            f"{label:12s} mean_stderr={np.mean(stderrs):.6e} mean_ci_width={np.mean(ci_widths):.6e} "
            f"mean_wall_s={np.mean(wall_times):.4f} paths={mc['n_paths']}"
        )

    plain = rows[0]
    anti = rows[1]
    if plain["mean_stderr"] > 0:
        ratio = anti["mean_stderr"] / plain["mean_stderr"]
        print(f"stderr_ratio_antithetic_over_plain={ratio:.4f}")
    if plain["mean_ci_width"] > 0:
        print(f"ci_width_ratio_antithetic_over_plain={anti['mean_ci_width'] / plain['mean_ci_width']:.4f}")

    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "method",
                "mean_stderr",
                "mean_ci_width",
                "mean_wall_s",
                "paths_effective",
                "ranks",
                "backend",
                "bs_price",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
