"""
Phase 2 (realism): estimate annualized volatility from historical OHLCV CSV, then run Monte Carlo
with baseline vs calibrated sigma side-by-side.

Data you need
-------------
- A CSV with **Date** (recommended) and **Adj Close** or **Close** (see `docs/DATA_SOURCES.md`).
- Quick path: `pip install yfinance pandas` then
  `python scripts/fetch_yahoo_ohlcv.py --ticker SPY --start 2019-01-01 --end 2024-12-31 --out data/raw/SPY.csv`
- Offline demo: `data/raw/example_synthetic_ohlcv.csv` in this repo.
"""

import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.black_scholes import black_scholes_call
from montecarlo.engine import petsc_monte_carlo_price
from montecarlo.historical_data import load_close_prices_from_csv
from montecarlo.payoffs import asian_call_payoff, barrier_up_and_out_call_payoff, european_call_payoff
from montecarlo.vol_calibration import annualized_volatility_from_closes, last_window


def parse_args():
    p = argparse.ArgumentParser(description="Phase 2: calibrate sigma from CSV, compare MC pricing.")
    p.add_argument("--csv", type=str, required=True, help="OHLCV CSV (Yahoo-style).")
    p.add_argument(
        "--window",
        type=int,
        default=0,
        help="Trailing calendar rows for vol (0 = use full series).",
    )
    p.add_argument("--trading-days", type=float, default=252.0)
    p.add_argument("--sigma-baseline", type=float, default=0.2, help="Synthetic benchmark sigma.")
    p.add_argument(
        "--s0",
        type=float,
        default=None,
        help="Spot for pricing; default = last close in CSV.",
    )
    p.add_argument("--strike", type=float, default=100.0)
    p.add_argument("--r", type=float, default=0.05)
    p.add_argument("--t", type=float, default=1.0)
    p.add_argument("--steps", type=int, default=252)
    p.add_argument("--paths", type=int, default=150000)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--option", choices=["european", "asian", "barrier"], default="european")
    p.add_argument("--barrier", type=float, default=130.0)
    p.add_argument("--antithetic", action="store_true")
    p.add_argument("--out-json", type=str, default="", help="Optional summary JSON path.")
    return p.parse_args()


def _mc_once(
    *,
    s0: float,
    strike: float,
    r: float,
    sigma: float,
    t: float,
    steps: int,
    n_paths: int,
    seed: int,
    option: str,
    barrier: float,
    antithetic: bool,
):
    if option == "european":
        payoff_builder = lambda k: lambda paths: european_call_payoff(paths, k)
        return petsc_monte_carlo_price(
            s0, strike, r, sigma, t, steps, n_paths, payoff_builder=payoff_builder, seed=seed, antithetic=antithetic
        )
    if option == "asian":
        payoff_builder = lambda k: lambda paths: asian_call_payoff(paths, k)
        return petsc_monte_carlo_price(
            s0, strike, r, sigma, t, steps, n_paths, payoff_builder=payoff_builder, seed=seed, antithetic=antithetic
        )
    payoff_builder = lambda k, b: lambda paths: barrier_up_and_out_call_payoff(paths, k, b)
    return petsc_monte_carlo_price(
        s0,
        strike,
        r,
        sigma,
        t,
        steps,
        n_paths,
        payoff_builder=payoff_builder,
        seed=seed,
        barrier=barrier,
        antithetic=antithetic,
    )


def main():
    args = parse_args()
    closes_all, meta = load_close_prices_from_csv(args.csv)
    series = last_window(closes_all, args.window)
    sigma_cal = annualized_volatility_from_closes(series, trading_days_per_year=args.trading_days)
    sigma_base = float(args.sigma_baseline)
    s0 = float(args.s0) if args.s0 is not None else float(series[-1])

    results = {}
    for label, sig in [("baseline", sigma_base), ("calibrated", sigma_cal)]:
        r_mc = _mc_once(
            s0=s0,
            strike=args.strike,
            r=args.r,
            sigma=sig,
            t=args.t,
            steps=args.steps,
            n_paths=args.paths,
            seed=args.seed,
            option=args.option,
            barrier=args.barrier,
            antithetic=args.antithetic,
        )
        if r_mc["rank"] == 0:
            row = {
                "sigma": sig,
                "price": r_mc["price"],
                "std_error": r_mc["std_error"],
                "ci_low": r_mc["ci_low"],
                "ci_high": r_mc["ci_high"],
                "n_paths": r_mc["n_paths"],
            }
            if args.option == "european":
                row["bs_price"] = black_scholes_call(s0, args.strike, args.r, sig, args.t)
            results[label] = row

    if r_mc["rank"] != 0:
        return

    summary = {
        "csv": str(pathlib.Path(args.csv).resolve()),
        "csv_meta": meta,
        "window": args.window,
        "trading_days_per_year": args.trading_days,
        "s0_used": s0,
        "sigma_baseline": sigma_base,
        "sigma_calibrated": sigma_cal,
        "option": args.option,
        "results": results,
    }

    print("=== Phase 2: historical vol + MC ===")
    print(f"csv={summary['csv']}")
    print(f"price_column={meta['price_column']} n_closes={meta['n_rows']} window={args.window or 'full'}")
    print(f"s0_used={s0:.6f} (from {'--s0' if args.s0 is not None else 'last close'})")
    print(f"sigma_baseline={sigma_base:.6f} sigma_calibrated={sigma_cal:.6f}")
    for label in ("baseline", "calibrated"):
        r = results[label]
        line = (
            f"{label:12s} sigma={r['sigma']:.6f} MC={r['price']:.8f} stderr={r['std_error']:.6e} "
            f"95%CI=[{r['ci_low']:.8f},{r['ci_high']:.8f}]"
        )
        if "bs_price" in r:
            line += f" BS={r['bs_price']:.8f}"
        print(line)

    if args.out_json:
        out = pathlib.Path(args.out_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(summary, indent=2))
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
