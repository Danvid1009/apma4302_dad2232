import numpy as np


def simulate_gbm_paths(
    s0: float,
    r: float,
    sigma: float,
    t: float,
    n_steps: int,
    n_paths: int,
    rng: np.random.Generator,
) -> np.ndarray:
    dt = t / n_steps
    drift = (r - 0.5 * sigma * sigma) * dt
    vol = sigma * np.sqrt(dt)

    z = rng.standard_normal((n_paths, n_steps))
    log_returns = drift + vol * z
    log_paths = np.cumsum(log_returns, axis=1)
    paths = s0 * np.exp(log_paths)

    # Include initial point S_0 at column 0.
    s0_col = np.full((n_paths, 1), s0, dtype=float)
    return np.hstack((s0_col, paths))


def simulate_gbm_paths_antithetic(
    s0: float,
    r: float,
    sigma: float,
    t: float,
    n_steps: int,
    n_pairs: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Simulate 2 * n_pairs GBM paths using n_pairs independent normal draws per step.

    For each draw Z, paths use +Z and -Z increments, producing an antithetic pair
    that shares RNG work and typically reduces payoff variance.
    """
    dt = t / n_steps
    drift = (r - 0.5 * sigma * sigma) * dt
    vol = sigma * np.sqrt(dt)

    z = rng.standard_normal((n_pairs, n_steps))
    log_a = np.cumsum(drift + vol * z, axis=1)
    log_b = np.cumsum(drift - vol * z, axis=1)
    paths_a = s0 * np.exp(log_a)
    paths_b = s0 * np.exp(log_b)
    s0_col = np.full((n_pairs, 1), s0, dtype=float)
    left = np.hstack((s0_col, paths_a))
    right = np.hstack((s0_col, paths_b))
    return np.vstack((left, right))
