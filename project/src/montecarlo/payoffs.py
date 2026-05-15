import numpy as np


def european_call_payoff(paths: np.ndarray, strike: float) -> np.ndarray:
    st = paths[:, -1]
    return np.maximum(st - strike, 0.0)


def asian_call_payoff(paths: np.ndarray, strike: float) -> np.ndarray:
    # Exclude initial point in average to avoid bias toward S0.
    mean_price = np.mean(paths[:, 1:], axis=1)
    return np.maximum(mean_price - strike, 0.0)


def barrier_up_and_out_call_payoff(paths: np.ndarray, strike: float, barrier: float) -> np.ndarray:
    st = paths[:, -1]
    knocked_out = np.any(paths[:, 1:] >= barrier, axis=1)
    vanilla = np.maximum(st - strike, 0.0)
    return np.where(knocked_out, 0.0, vanilla)
