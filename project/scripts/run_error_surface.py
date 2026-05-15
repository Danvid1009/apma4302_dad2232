"""
European call: grid in (N, sigma). Record MC price, Black–Scholes, absolute error, stderr, wall time.

Export a flat CSV suitable for ParaView Table To Points / 2D plotting (methods + HPC workload preview).
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

from montecarlo.black_scholes import black_scholes_call
from montecarlo.engine import petsc_monte_carlo_price
from montecarlo.payoffs import european_call_payoff


def parse_args():
    p = argparse.ArgumentParser(description="European MC error surface over (N, sigma).")
    p.add_argument("--s0", type=float, default=100.0)
    p.add_argument("--strike", type=float, default=100.0)
    p.add_argument("--r", type=float, default=0.05)
    p.add_argument("--t", type=float, default=1.0)
    p.add_argument("--steps", type=int, default=252)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument(
        "--paths-list",
        type=int,
        nargs="+",
        default=[25000, 100000, 400000],
        help="Path counts N for the surface.",
    )
    p.add_argument(
        "--sigma-list",
        type=float,
        nargs="+",
        default=[0.12, 0.2, 0.35],
        help="Volatilities sigma for the surface.",
    )
    p.add_argument("--antithetic", action="store_true")
    p.add_argument("--out", type=str, default="output/error_surface_european.csv")
    return p.parse_args()


def main():
    args = parse_args()
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    payoff_builder = lambda strike: lambda paths: european_call_payoff(paths, strike)

    rows = []
    rank0 = True
    for sigma in args.sigma_list:
        bs = black_scholes_call(args.s0, args.strike, args.r, float(sigma), args.t)
        for n_paths in args.paths_list:
            t0 = time.perf_counter()
            result = petsc_monte_carlo_price(
                args.s0,
                args.strike,
                args.r,
                float(sigma),
                args.t,
                args.steps,
                int(n_paths),
                payoff_builder=payoff_builder,
                seed=args.seed,
                antithetic=args.antithetic,
            )
            wall = time.perf_counter() - t0
            if result["rank"] != 0:
                rank0 = False
                continue
            err = abs(result["price"] - bs)
            rows.append(
                {
                    "n_paths": result["n_paths"],
                    "sigma": float(sigma),
                    "mc_price": result["price"],
                    "bs_price": bs,
                    "abs_err": err,
                    "std_error": result["std_error"],
                    "ci_width": result["ci_high"] - result["ci_low"],
                    "wall_s": wall,
                    "backend": result["backend"],
                    "ranks": result["size"],
                }
            )
            print(
                f"N={result['n_paths']} sigma={sigma:.3f} mc={result['price']:.6f} bs={bs:.6f} "
                f"abs_err={err:.6e} wall_s={wall:.4f}"
            )

    if not rank0 or not rows:
        return

    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "n_paths",
                "sigma",
                "mc_price",
                "bs_price",
                "abs_err",
                "std_error",
                "ci_width",
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
