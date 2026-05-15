import argparse
import csv
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.black_scholes import black_scholes_call
from montecarlo.engine import petsc_monte_carlo_price
from montecarlo.payoffs import european_call_payoff


def parse_args():
    parser = argparse.ArgumentParser(description="Run European convergence sweep and fit error slope.")
    parser.add_argument("--paths-list", type=int, nargs="+", default=[20000, 50000, 100000, 200000, 400000])
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--steps", type=int, default=252)
    parser.add_argument("--s0", type=float, default=100.0)
    parser.add_argument("--strike", type=float, default=100.0)
    parser.add_argument("--r", type=float, default=0.05)
    parser.add_argument("--sigma", type=float, default=0.2)
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--antithetic", action="store_true", help="Antithetic path pairs per replicate.")
    parser.add_argument("--out", type=str, default="output/convergence.csv")
    return parser.parse_args()


def main():
    args = parse_args()
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    bs = black_scholes_call(args.s0, args.strike, args.r, args.sigma, args.t)
    payoff_builder = lambda strike: lambda paths: european_call_payoff(paths, strike)

    rows = []
    rank0 = True
    size = 1
    backend = "serial"

    for n_paths in args.paths_list:
        rep_errors = []
        rep_prices = []
        rep_stderr = []
        result = None

        for rep in range(args.reps):
            result = petsc_monte_carlo_price(
                args.s0,
                args.strike,
                args.r,
                args.sigma,
                args.t,
                args.steps,
                n_paths,
                payoff_builder=payoff_builder,
                seed=args.seed + rep * 9973,
                antithetic=args.antithetic,
            )

            if result["rank"] != 0:
                rank0 = False
                continue

            err = result["price"] - bs
            rep_errors.append(err)
            rep_prices.append(result["price"])
            rep_stderr.append(result["std_error"])

        if not rank0 or not rep_errors:
            continue

        abs_error = float(abs(np.mean(rep_errors)))
        mean_abs_error = float(np.mean(np.abs(rep_errors)))
        rmse_error = float(np.sqrt(np.mean(np.square(rep_errors))))
        mean_price = float(np.mean(rep_prices))
        mean_stderr = float(np.mean(rep_stderr))
        rows.append(
            {
                "paths": n_paths,
                "mc_price": mean_price,
                "bs_price": bs,
                "abs_error": abs_error,
                "mean_abs_error": mean_abs_error,
                "rmse_error": rmse_error,
                "std_error": mean_stderr,
                "ci_low": result["ci_low"],
                "ci_high": result["ci_high"],
                "backend": result["backend"],
                "ranks": result["size"],
                "antithetic": bool(args.antithetic),
            }
        )
        size = result["size"]
        backend = result["backend"]
        print(
            f"paths={n_paths} mc_mean={mean_price:.8f} bs={bs:.8f} "
            f"mean_abs_error={mean_abs_error:.8e} rmse_error={rmse_error:.8e}"
        )

    if not rank0:
        return

    # Fit log(error) = a + b log(N), expected b ~ -0.5.
    n = np.array([row["paths"] for row in rows], dtype=float)
    e_abs = np.array([max(row["mean_abs_error"], 1e-16) for row in rows], dtype=float)
    e_rmse = np.array([max(row["rmse_error"], 1e-16) for row in rows], dtype=float)
    slope_abs, intercept_abs = np.polyfit(np.log(n), np.log(e_abs), 1)
    slope_rmse, intercept_rmse = np.polyfit(np.log(n), np.log(e_rmse), 1)

    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "paths",
                "mc_price",
                "bs_price",
                "abs_error",
                "mean_abs_error",
                "rmse_error",
                "std_error",
                "ci_low",
                "ci_high",
                "backend",
                "ranks",
                "antithetic",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"fitted_loglog_slope_mean_abs={slope_abs:.6f}")
    print(f"fitted_loglog_intercept_mean_abs={intercept_abs:.6f}")
    print(f"fitted_loglog_slope_rmse={slope_rmse:.6f}")
    print(f"fitted_loglog_intercept_rmse={intercept_rmse:.6f}")
    print(f"expected_slope_near=-0.5")
    print(f"backend={backend} ranks={size}")
    print(f"reps={args.reps}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
