import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.black_scholes import black_scholes_call
from montecarlo.engine import petsc_monte_carlo_price
from montecarlo.payoffs import european_call_payoff


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate PETSc parallel Monte Carlo European call vs Black-Scholes."
    )
    parser.add_argument("--s0", type=float, default=100.0)
    parser.add_argument("--strike", type=float, default=100.0)
    parser.add_argument("--r", type=float, default=0.05)
    parser.add_argument("--sigma", type=float, default=0.2)
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=252)
    parser.add_argument("--paths", type=int, default=200000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument(
        "--antithetic",
        action="store_true",
        help="Use antithetic path pairs (total paths forced even).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    payoff_builder = lambda strike: lambda paths: european_call_payoff(paths, strike)
    mc = petsc_monte_carlo_price(
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

    if mc["rank"] != 0:
        return

    bs = black_scholes_call(args.s0, args.strike, args.r, args.sigma, args.t)
    abs_err = abs(mc["price"] - bs)
    rel_err = abs_err / abs(bs) if bs != 0.0 else float("nan")

    print(f"mc_price={mc['price']:.8f}")
    print(f"bs_price={bs:.8f}")
    print(f"abs_error={abs_err:.8e}")
    print(f"rel_error={rel_err:.8e}")
    print(f"std_error={mc['std_error']:.8e}")
    print(f"95% CI=[{mc['ci_low']:.8f}, {mc['ci_high']:.8f}]")
    print(f"backend={mc['backend']} ranks={mc['size']}")


if __name__ == "__main__":
    main()
