#!/usr/bin/env python3
"""
Compare European MC prices to quoted mids (CSV template: data/raw/quotes_template.csv).

Each row: spot, strike, T_years, r, mid, sigma_mc — runs one MC and prints absolute / relative error.
"""

import argparse
import csv
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.engine import petsc_monte_carlo_price
from montecarlo.payoffs import european_call_payoff


def parse_args():
    p = argparse.ArgumentParser(description="MC vs market mid quotes (European).")
    p.add_argument("--quotes", type=str, default="data/raw/quotes_template.csv")
    p.add_argument("--steps", type=int, default=252)
    p.add_argument("--paths", type=int, default=250000)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--antithetic", action="store_true")
    p.add_argument("--out", type=str, default="output/mc_vs_quotes.csv")
    return p.parse_args()


def main():
    args = parse_args()
    qpath = pathlib.Path(args.quotes)
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payoff_builder = lambda k: lambda paths: european_call_payoff(paths, k)

    rows_out = []
    with qpath.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not row.get("spot"):
                continue
            try:
                s0 = float(row["spot"])
            except (TypeError, ValueError, KeyError):
                continue
            strike = float(row["strike"])
            t = float(row["T_years"])
            r = float(row["r"])
            mid = float(row["mid"])
            sigma = float(row["sigma_mc"])
            res = petsc_monte_carlo_price(
                s0,
                strike,
                r,
                sigma,
                t,
                args.steps,
                args.paths,
                payoff_builder=payoff_builder,
                seed=args.seed,
                antithetic=args.antithetic,
            )
            if res["rank"] != 0:
                return
            mc = res["price"]
            err = mc - mid
            rel = err / abs(mid) if mid != 0 else float("nan")
            rows_out.append(
                {
                    "spot": s0,
                    "strike": strike,
                    "T_years": t,
                    "r": r,
                    "mid": mid,
                    "sigma_mc": sigma,
                    "mc_price": mc,
                    "abs_err": abs(err),
                    "rel_err": rel,
                    "std_error": res["std_error"],
                }
            )
            print(
                f"K={strike} T={t} mid={mid:.4f} mc={mc:.4f} abs_err={abs(err):.4f} "
                f"rel_err={rel:.4f} stderr={res['std_error']:.4e}"
            )

    if not rows_out:
        print("No quote rows processed.", file=sys.stderr)
        return 1

    with out_path.open("w", newline="") as g:
        w = csv.DictWriter(
            g,
            fieldnames=[
                "spot",
                "strike",
                "T_years",
                "r",
                "mid",
                "sigma_mc",
                "mc_price",
                "abs_err",
                "rel_err",
                "std_error",
            ],
        )
        w.writeheader()
        w.writerows(rows_out)
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
