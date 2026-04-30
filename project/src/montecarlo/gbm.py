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
