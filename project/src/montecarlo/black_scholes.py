import numpy as np
from scipy.stats import norm


def black_scholes_call(s0: float, strike: float, r: float, sigma: float, t: float) -> float:
    if t <= 0.0:
        return max(s0 - strike, 0.0)
    if sigma <= 0.0:
        forward = s0 * np.exp(r * t)
        return np.exp(-r * t) * max(forward - strike, 0.0)

    d1 = (np.log(s0 / strike) + (r + 0.5 * sigma * sigma) * t) / (sigma * np.sqrt(t))
    d2 = d1 - sigma * np.sqrt(t)
    return s0 * norm.cdf(d1) - strike * np.exp(-r * t) * norm.cdf(d2)


def black_scholes_call_vega(s0: float, strike: float, r: float, sigma: float, t: float) -> float:
    """∂C/∂σ (Black–Scholes European call). Used e.g. for SNES / Newton Jacobians."""
    if t <= 0.0 or sigma <= 0.0:
        return 0.0
    d1 = (np.log(s0 / strike) + (r + 0.5 * sigma * sigma) * t) / (sigma * np.sqrt(t))
    return float(s0 * norm.pdf(d1) * np.sqrt(t))
