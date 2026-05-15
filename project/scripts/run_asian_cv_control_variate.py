#!/usr/bin/env python3
"""Compare plain Asian MC vs Asian with European control variate (same paths, MPI-aware)."""

import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.engine import petsc_asian_call_cv_european_control, petsc_monte_carlo_price
from montecarlo.payoffs import asian_call_payoff


def parse_args():
    p = argparse.ArgumentParser(description="Asian call: naive MC vs European control variate.")
    p.add_argument("--s0", type=float, default=100.0)
    p.add_argument("--strike", type=float, default=100.0)
    p.add_argument("--r", type=float, default=0.05)
    p.add_argument("--sigma", type=float, default=0.2)
    p.add_argument("--t", type=float, default=1.0)
    p.add_argument("--steps", type=int, default=252)
    p.add_argument("--paths", type=int, default=200000)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--antithetic", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    payoff_builder = lambda k: lambda paths: asian_call_payoff(paths, k)
    naive = petsc_monte_carlo_price(
        args.s0,
        args.strike,
        args.r,
        args.sigma,
        args.t,
        args.steps,
        args.paths,
        payoff_builder=payoff_builder,
        seed=args.seed,
        antithetic=args.antithetic,
    )
    cv = petsc_asian_call_cv_european_control(
        args.s0,
        args.strike,
        args.r,
        args.sigma,
        args.t,
        args.steps,
        args.paths,
        seed=args.seed,
        antithetic=args.antithetic,
    )
    if naive["rank"] != 0:
        return
    print("=== Asian call: naive MC vs European control variate ===")
    print(f"naive_price={naive['price']:.8f} naive_stderr={naive['std_error']:.8e}")
    print(
        f"cv_price={cv['price_cv']:.8f} cv_stderr={cv['std_error_cv']:.8e} "
        f"c={cv['cv_coef_c']:.6f} bs_european={cv['european_bs']:.8f}"
    )
    print(f"stderr_ratio_cv_over_naive={cv['std_error_cv'] / max(naive['std_error'], 1e-20):.4f}")


if __name__ == "__main__":
    main()
