"""
Sweep barrier level B for an up-and-out barrier call: price, stderr, CI width, wall time.

HPC angle: payoff variance and effective "hardness" of MC change with knock-out geometry;
pair with scaling experiments in the registry.
"""

import argparse
import csv
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.engine import petsc_monte_carlo_price
from montecarlo.payoffs import barrier_up_and_out_call_payoff


def parse_args():
    p = argparse.ArgumentParser(description="Sweep barrier level for MC barrier option.")
    p.add_argument("--s0", type=float, default=100.0)
    p.add_argument("--strike", type=float, default=100.0)
    p.add_argument("--r", type=float, default=0.05)
    p.add_argument("--sigma", type=float, default=0.2)
    p.add_argument("--t", type=float, default=1.0)
    p.add_argument("--steps", type=int, default=252)
    p.add_argument("--paths", type=int, default=150000)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument(
        "--barriers",
        type=float,
        nargs="+",
        default=[108.0, 115.0, 122.0, 130.0, 140.0, 155.0, 175.0],
        help="Barrier levels B to test (up-and-out).",
    )
    p.add_argument("--antithetic", action="store_true")
    p.add_argument("--out", type=str, default="output/barrier_sweep.csv")
    return p.parse_args()


def main():
    args = parse_args()
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    payoff_builder = lambda strike, barrier: lambda paths: barrier_up_and_out_call_payoff(paths, strike, barrier)

    rows = []
    rank0 = True
    for barrier in args.barriers:
        t0 = time.perf_counter()
        result = petsc_monte_carlo_price(
            args.s0,
            args.strike,
            args.r,
            args.sigma,
            args.t,
            args.steps,
            args.paths,
            payoff_builder=payoff_builder,
            seed=args.seed,
            barrier=float(barrier),
            antithetic=args.antithetic,
        )
        wall = time.perf_counter() - t0
        if result["rank"] != 0:
            rank0 = False
            continue
        ci_w = result["ci_high"] - result["ci_low"]
        rows.append(
            {
                "barrier": float(barrier),
                "price": result["price"],
                "std_error": result["std_error"],
                "ci_width": ci_w,
                "n_paths": result["n_paths"],
                "wall_s": wall,
                "backend": result["backend"],
                "ranks": result["size"],
            }
        )
        print(
            f"B={barrier:.3f} price={result['price']:.6f} stderr={result['std_error']:.6e} "
            f"CI_width={ci_w:.6e} wall_s={wall:.4f}"
        )

    if not rank0 or not rows:
        return

    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "barrier",
                "price",
                "std_error",
                "ci_width",
                "n_paths",
                "wall_s",
                "backend",
                "ranks",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
