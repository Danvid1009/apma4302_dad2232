"""Historical volatility from log returns (numpy only)."""

from __future__ import annotations

import numpy as np


def annualized_volatility_from_closes(
    closes: np.ndarray,
    *,
    trading_days_per_year: float = 252.0,
) -> float:
    """
    Sample standard deviation of log returns, annualized by sqrt(trading_days_per_year).

    Parameters
    ----------
    closes
        Strictly positive prices, oldest to newest (e.g. daily).
    """
    x = np.asarray(closes, dtype=np.float64).ravel()
    x = x[np.isfinite(x) & (x > 0.0)]
    if x.size < 2:
        raise ValueError("Need at least two positive closes for volatility.")
    lr = np.diff(np.log(x))
    if lr.size < 1:
        raise ValueError("Need at least two closes.")
    sigma_daily = float(np.std(lr, ddof=1))
    return float(sigma_daily * np.sqrt(float(trading_days_per_year)))


def last_window(closes: np.ndarray, window: int) -> np.ndarray:
    """Return the last `window` closes; if window <= 0, return full series."""
    if window is None or window <= 0:
        return closes
    w = int(window)
    if w >= closes.size:
        return closes
    return closes[-w:]


__all__ = ["annualized_volatility_from_closes", "last_window"]
