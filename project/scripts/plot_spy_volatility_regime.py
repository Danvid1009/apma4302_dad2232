"""
Phase 2 visual: SPY adjusted close + rolling annualized volatility (twin axes, dark theme).
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

_scripts = pathlib.Path(__file__).resolve().parent
if str(_scripts) not in sys.path:
    sys.path.insert(0, str(_scripts))
from plot_archive import archive_prior_png


def parse_args():
    p = argparse.ArgumentParser(description="Plot SPY price + rolling annualized vol.")
    p.add_argument("--csv", type=str, default="data/raw/SPY.csv")
    p.add_argument("--window", type=int, default=252, help="Rolling window in trading days.")
    p.add_argument("--out", type=str, default="visuals/phase2_spy_price_and_vol.png")
    return p.parse_args()


def load_date_close(path: pathlib.Path):
    dates, closes = [], []
    with path.open(newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            d = row.get("Date") or row.get("date")
            ac = row.get("Adj Close") or row.get("AdjClose") or row.get("Close")
            if not d or not ac:
                continue
            try:
                closes.append(float(ac))
                dates.append(datetime.strptime(str(d)[:10], "%Y-%m-%d"))
            except (ValueError, TypeError):
                continue
    return np.array(dates), np.array(closes, dtype=np.float64)


def rolling_ann_vol(closes: np.ndarray, win: int, trading_days: float = 252.0) -> np.ndarray:
    n = closes.size
    out = np.full(n, np.nan)
    lr = np.diff(np.log(closes))
    for i in range(win, n):
        seg = lr[i - win : i]
        if seg.size < 2:
            continue
        out[i] = float(np.std(seg, ddof=1) * np.sqrt(trading_days))
    return out


def main():
    args = parse_args()
    path = pathlib.Path(args.csv)
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    dates, closes = load_date_close(path)
    vol = rolling_ann_vol(closes, args.window)

    plt.rcParams.update(
        {
            "figure.facecolor": "#0d1117",
            "axes.facecolor": "#161b22",
            "axes.edgecolor": "#30363d",
            "axes.labelcolor": "#c9d1d9",
            "text.color": "#c9d1d9",
            "xtick.color": "#8b949e",
            "ytick.color": "#8b949e",
            "grid.color": "#30363d",
            "grid.alpha": 0.5,
        }
    )

    fig, ax1 = plt.subplots(figsize=(13, 5.5), constrained_layout=True)
    # Subtle context band (COVID crash window) — optional visual anchor for the vol spike.
    covid0 = datetime(2020, 2, 15)
    covid1 = datetime(2020, 5, 15)
    ax1.axvspan(covid0, covid1, color="#f0883e", alpha=0.08, zorder=0, linewidth=0)

    ax1.fill_between(dates, closes, alpha=0.10, color="#58a6ff", zorder=1)
    ax1.plot(dates, closes, color="#79c0ff", lw=1.6, label="SPY Adj Close", zorder=2)
    ax1.set_ylabel("Price (USD)", color="#79c0ff", fontweight="medium")
    ax1.tick_params(axis="y", colors="#79c0ff")
    ax1.spines["left"].set_color("#79c0ff")
    ax1.spines["left"].set_alpha(0.65)
    ax1.set_title("SPY — price path and rolling realized volatility (Phase 2 inputs)", color="#f0f6fc", fontsize=13, pad=12)

    ax2 = ax1.twinx()
    ax2.plot(dates, vol, color="#ffa657", lw=1.85, label=f"{args.window}d ann. vol", zorder=2)
    ax2.set_ylabel("Annualized $\\sigma$ (rolling)", color="#ffa657", fontweight="medium")
    ax2.tick_params(axis="y", colors="#ffa657")
    ax2.spines["right"].set_color("#ffa657")
    ax2.spines["right"].set_alpha(0.65)

    ax1.xaxis.set_major_locator(mdates.YearLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax1.grid(True)
    plt.setp(ax1.get_xticklabels(), rotation=22, ha="right")

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper left", framealpha=0.22, fontsize=9)

    ax1.annotate(
        f"{args.window}-day vol is undefined until the window fills. Orange tint: Feb–May 2020.",
        xy=(0, 0),
        xycoords="axes fraction",
        xytext=(0, -0.14),
        textcoords="axes fraction",
        fontsize=8,
        color="#8b949e",
        ha="left",
        va="top",
        clip_on=False,
    )

    archive_prior_png(out)
    fig.savefig(out, dpi=200, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
