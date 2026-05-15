#!/usr/bin/env python3
"""
Download daily OHLCV from Yahoo (via yfinance) for volatility calibration experiments.

Install once (optional; listed in requirements-phase2.txt):
  pip install -r requirements-phase2.txt

Example:
  python scripts/fetch_yahoo_ohlcv.py --ticker SPY --start 2020-01-01 --end 2024-12-31 --out data/raw/SPY.csv
"""

import argparse
import pathlib
import sys


def parse_args():
    p = argparse.ArgumentParser(description="Fetch Yahoo OHLCV via yfinance.")
    p.add_argument("--ticker", type=str, default="SPY")
    p.add_argument("--start", type=str, required=True)
    p.add_argument("--end", type=str, required=True)
    p.add_argument("--out", type=str, required=True)
    return p.parse_args()


def main():
    args = parse_args()
    try:
        import pandas as pd
        import yfinance as yf
    except ImportError:
        print("Missing dependency: pip install -r requirements-phase2.txt", file=sys.stderr)
        return 1

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    df = yf.download(args.ticker, start=args.start, end=args.end, progress=False, auto_adjust=False)
    if df.empty:
        print("No rows returned; check ticker and date range.", file=sys.stderr)
        return 2

    if isinstance(df.columns, pd.MultiIndex):
        df = df.copy()
        df.columns = df.columns.get_level_values(0)
    df = df.reset_index()
    if "Date" in df.columns:
        df["Date"] = pd.to_datetime(df["Date"]).dt.strftime("%Y-%m-%d")
    df.to_csv(out, index=False)
    print(f"wrote {out} rows={len(df)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
