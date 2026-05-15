"""
Load close prices from a CSV (e.g. Yahoo / yfinance export) for volatility calibration.

Header-based columns. Prefers **Adj Close**, then **Close** (case-insensitive).
Sorts by **Date** when that column exists and parses.
"""

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np


def _norm_key(s: str) -> str:
    return re.sub(r"\s+", "", str(s).strip().lower())


def _parse_date(s: str) -> Optional[datetime]:
    s = str(s).strip()
    if not s:
        return None
    for fmt in ("%Y-%m-%d", "%Y/%m/%d", "%m/%d/%Y"):
        try:
            return datetime.strptime(s[:10], fmt)
        except ValueError:
            continue
    try:
        return datetime.fromisoformat(s[:10])
    except ValueError:
        return None


def _pick_header(fieldnames: list[str], *candidates: str) -> Optional[str]:
    cmap = {_norm_key(h): h for h in fieldnames}
    for c in candidates:
        if _norm_key(c) in cmap:
            return cmap[_norm_key(c)]
    return None


def load_close_prices_from_csv(
    path: str | Path,
    *,
    prefer_adj_close: bool = True,
) -> tuple[np.ndarray, dict]:
    """
    Returns
    -------
    closes : ndarray (n,), oldest -> newest
    meta : dict with path, n_rows, price_column, sorted_by_date
    """
    path = Path(path)
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        fieldnames = list(reader.fieldnames)
        records = list(reader)

    price_h = None
    if prefer_adj_close:
        price_h = _pick_header(fieldnames, "Adj Close", "AdjClose", "Adjusted Close", "adjusted_close")
    if price_h is None:
        price_h = _pick_header(fieldnames, "Close")
    if price_h is None:
        raise ValueError(f"No Adj Close or Close column in {path}; found {fieldnames}")

    date_h = _pick_header(fieldnames, "Date", "Datetime", "timestamp")

    parsed: list[tuple[Optional[datetime], float]] = []
    for rec in records:
        if price_h not in rec or rec[price_h] is None or str(rec[price_h]).strip() == "":
            continue
        try:
            px = float(rec[price_h])
        except (TypeError, ValueError):
            continue
        if not np.isfinite(px) or px <= 0:
            continue
        dt: Optional[datetime] = None
        if date_h is not None and date_h in rec and rec[date_h] not in (None, ""):
            dt = _parse_date(str(rec[date_h]))
        parsed.append((dt, px))

    if not parsed:
        raise ValueError(f"No valid price rows in {path}")

    sortable = date_h is not None and all(t is not None for t, _ in parsed)
    if sortable:
        parsed.sort(key=lambda x: x[0])
    closes = np.array([p for _, p in parsed], dtype=np.float64)

    meta = {
        "path": str(path.resolve()),
        "n_rows": int(closes.size),
        "price_column": price_h,
        "sorted_by_date": sortable,
    }
    return closes, meta


def load_adj_close_series_from_csv(
    path: str | Path,
    *,
    prefer_adj_close: bool = True,
) -> tuple[list[str], np.ndarray, dict]:
    """
    Like ``load_close_prices_from_csv`` but also returns ISO date strings per row (oldest -> newest).
    If dates are missing or unsortable, returns empty strings for those rows.
    """
    path = Path(path)
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        fieldnames = list(reader.fieldnames)
        records = list(reader)

    price_h = None
    if prefer_adj_close:
        price_h = _pick_header(fieldnames, "Adj Close", "AdjClose", "Adjusted Close", "adjusted_close")
    if price_h is None:
        price_h = _pick_header(fieldnames, "Close")
    if price_h is None:
        raise ValueError(f"No Adj Close or Close column in {path}; found {fieldnames}")

    date_h = _pick_header(fieldnames, "Date", "Datetime", "timestamp")

    parsed: list[tuple[Optional[datetime], float]] = []
    for rec in records:
        if price_h not in rec or rec[price_h] is None or str(rec[price_h]).strip() == "":
            continue
        try:
            px = float(rec[price_h])
        except (TypeError, ValueError):
            continue
        if not np.isfinite(px) or px <= 0:
            continue
        dt: Optional[datetime] = None
        if date_h is not None and date_h in rec and rec[date_h] not in (None, ""):
            dt = _parse_date(str(rec[date_h]))
        parsed.append((dt, px))

    if not parsed:
        raise ValueError(f"No valid price rows in {path}")

    sortable = date_h is not None and all(t is not None for t, _ in parsed)
    if sortable:
        parsed.sort(key=lambda x: x[0])

    date_strs: list[str] = []
    for t, _ in parsed:
        date_strs.append(t.strftime("%Y-%m-%d") if t is not None else "")
    closes = np.array([p for _, p in parsed], dtype=np.float64)

    meta = {
        "path": str(path.resolve()),
        "n_rows": int(closes.size),
        "price_column": price_h,
        "sorted_by_date": sortable,
    }
    return date_strs, closes, meta


__all__ = ["load_close_prices_from_csv", "load_adj_close_series_from_csv"]
