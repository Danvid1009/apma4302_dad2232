#!/usr/bin/env python3
"""
Roll-forward Phase 2: along a historical price series, re-estimate trailing annualized
vol every ``--step`` trading days and re-price a European call (constant T, r).

Connects Phase 1 (same MC + MPI) to Phase 2 (time-varying calibrated inputs).
"""

import argparse
import csv
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.black_scholes import black_scholes_call
from montecarlo.engine import petsc_monte_carlo_price
from montecarlo.historical_data import load_adj_close_series_from_csv
from montecarlo.payoffs import european_call_payoff
from montecarlo.vol_calibration import annualized_volatility_from_closes


def parse_args():
    p = argparse.ArgumentParser(description="Roll-forward vol + European repricing along history.")
    p.add_argument("--csv", type=str, default="data/raw/SPY.csv")
    p.add_argument("--vol-window", type=int, default=252)
    p.add_argument("--step", type=int, default=42, help="Trading days between snapshots (~2 months at 21d/mo).")
    p.add_argument("--t", type=float, default=0.25)
    p.add_argument("--r", type=float, default=0.045)
    p.add_argument("--strike-pct", type=float, default=1.0, help="Strike = pct * spot (1.0 = ATM).")
    p.add_argument("--paths", type=int, default=40000)
    p.add_argument("--steps", type=int, default=252)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--out", type=str, default="output/phase2_roll_forward.csv")
    return p.parse_args()


def load_dates_closes(path: pathlib.Path):
    return load_adj_close_series_from_csv(path)


def main():
    args = parse_args()
    path = pathlib.Path(args.csv)
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    dates, closes, _meta = load_dates_closes(path)
    n = closes.size
    w = int(args.vol_window)
    payoff_builder = lambda k: lambda paths: european_call_payoff(paths, k)

    rows = []
    for i in range(w, n, int(args.step)):
        seg = closes[i - w : i]
        sig = annualized_volatility_from_closes(seg)
        s0 = float(closes[i])
        strike = float(args.strike_pct * s0)
        res = petsc_monte_carlo_price(
            s0,
            strike,
            args.r,
            sig,
            args.t,
            args.steps,
            args.paths,
            payoff_builder=payoff_builder,
            seed=args.seed + i,
        )
        if res["rank"] != 0:
            return
        bs = black_scholes_call(s0, strike, args.r, sig, args.t)
        rows.append(
            {
                "asof_date": dates[i] if i < len(dates) else str(i),
                "s0": s0,
                "sigma_roll": sig,
                "strike": strike,
                "mc_price": res["price"],
                "bs_price": bs,
                "std_error": res["std_error"],
            }
        )
        print(f"{dates[i] if i < len(dates) else i} s0={s0:.2f} sig={sig:.4f} mc={res['price']:.4f} bs={bs:.4f}")

    if not rows:
        return

    with out.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["asof_date", "s0", "sigma_roll", "strike", "mc_price", "bs_price", "std_error"],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
