import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.engine import petsc_monte_carlo_price
from montecarlo.payoffs import asian_call_payoff, barrier_up_and_out_call_payoff, european_call_payoff


def parse_args():
    parser = argparse.ArgumentParser(description="Run PETSc-based parallel Monte Carlo option pricing.")
    parser.add_argument("--option", choices=["european", "asian", "barrier"], default="european")
    parser.add_argument("--s0", type=float, default=100.0)
    parser.add_argument("--strike", type=float, default=100.0)
    parser.add_argument("--r", type=float, default=0.05)
    parser.add_argument("--sigma", type=float, default=0.2)
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=252)
    parser.add_argument("--paths", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--barrier", type=float, default=130.0)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.option == "european":
        payoff_builder = lambda strike: lambda paths: european_call_payoff(paths, strike)
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
        )
    elif args.option == "asian":
        payoff_builder = lambda strike: lambda paths: asian_call_payoff(paths, strike)
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
        )
    else:
        payoff_builder = lambda strike, barrier: lambda paths: barrier_up_and_out_call_payoff(
            paths, strike, barrier
        )
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
            barrier=args.barrier,
        )

    if result["rank"] == 0:
        print(f"option={args.option}")
        print(f"price={result['price']:.8f}")
        print(f"std_error={result['std_error']:.8f}")
        print(f"95% CI=[{result['ci_low']:.8f}, {result['ci_high']:.8f}]")
        print(f"paths={result['n_paths']} ranks={result['size']} backend={result['backend']}")


if __name__ == "__main__":
    main()
