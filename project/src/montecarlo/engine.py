from typing import Callable, Iterator, Optional

import numpy as np

from .gbm import simulate_gbm_paths, simulate_gbm_paths_antithetic

try:
    from petsc4py import PETSc
except ImportError:
    PETSc = None

try:
    from mpi4py import MPI
except ImportError:
    MPI = None

import os
from contextlib import contextmanager


def _mc_petsc_log_events() -> bool:
    """Opt-in PETSc Log.Event regions around pathgen / payoff / reductions."""
    v = os.environ.get("MC_PETSC_LOG_EVENTS", "")
    return v.lower() in ("1", "true", "yes", "on")


@contextmanager
def _petsc_log_event(name: str) -> Iterator[None]:
    if PETSc is None or not _mc_petsc_log_events():
        yield
        return
    try:
        ev = PETSc.Log.Event(name)
    except Exception:
        yield
        return
    try:
        ev.begin()
        yield
    finally:
        try:
            ev.end()
        except Exception:
            pass


PayoffFn = Callable[[np.ndarray], np.ndarray]


def _petsc_comm():
    if PETSc is None:
        return None, None, 0, 1
    comm = PETSc.COMM_WORLD
    mpi_comm = comm.tompi4py() if MPI is not None else None
    return comm, mpi_comm, comm.getRank(), comm.getSize()


def _split_paths(n_paths: int, size: int, rank: int) -> int:
    base = n_paths // size
    rem = n_paths % size
    return base + (1 if rank < rem else 0)


def petsc_monte_carlo_price(
    s0: float,
    strike: float,
    r: float,
    sigma: float,
    t: float,
    n_steps: int,
    n_paths: int,
    payoff_builder: Callable[[float], PayoffFn],
    seed: int = 12345,
    barrier: Optional[float] = None,
    antithetic: bool = False,
):
    petsc_comm, mpi_comm, rank, size = _petsc_comm()
    if antithetic:
        n_paths = max((n_paths // 2) * 2, 2)
        n_pairs_total = n_paths // 2
        local_pairs = _split_paths(n_pairs_total, size, rank)
    else:
        local_pairs = 0
    local_seed = seed + 10007 * rank
    rng = np.random.default_rng(local_seed)

    with _petsc_log_event("MC_paths"):
        if antithetic:
            paths = simulate_gbm_paths_antithetic(s0, r, sigma, t, n_steps, local_pairs, rng)
        else:
            local_n = _split_paths(n_paths, size, rank)
            paths = simulate_gbm_paths(s0, r, sigma, t, n_steps, local_n, rng)
    with _petsc_log_event("MC_payoff"):
        if barrier is None:
            payoff = payoff_builder(strike)(paths)
        else:
            payoff = payoff_builder(strike, barrier)(paths)
        discounted = np.exp(-r * t) * payoff

    with _petsc_log_event("MC_reduce"):
        local_sum = float(np.sum(discounted))
        local_sumsq = float(np.sum(discounted * discounted))
        local_count = int(discounted.size)

        if mpi_comm is None:
            total_sum, total_sumsq, total_count = local_sum, local_sumsq, local_count
        else:
            total_sum = mpi_comm.allreduce(local_sum, op=MPI.SUM)
            total_sumsq = mpi_comm.allreduce(local_sumsq, op=MPI.SUM)
            total_count = mpi_comm.allreduce(local_count, op=MPI.SUM)

    mean = total_sum / total_count
    second_moment = total_sumsq / total_count
    variance = max(second_moment - mean * mean, 0.0)
    std_error = np.sqrt(variance / total_count)
    ci_half = 1.96 * std_error

    return {
        "price": mean,
        "std_error": std_error,
        "ci_low": mean - ci_half,
        "ci_high": mean + ci_half,
        "n_paths": total_count,
        "rank": rank,
        "size": size,
        "backend": "petsc" if petsc_comm is not None else "serial",
    }


def petsc_asian_call_cv_european_control(
    s0: float,
    strike: float,
    r: float,
    sigma: float,
    t: float,
    n_steps: int,
    n_paths: int,
    seed: int = 12345,
    antithetic: bool = False,
):
    """
    Arithmetic Asian call with a European-call control variate.

    Uses the same GBM paths to form discounted European payoff X and Asian payoff Y,
    then the CV estimator mean(Z) with Z_i = Y_i + c (X_i - E[X]) and E[X] = Black–Scholes
    European price (same strike/maturity inputs).
    """
    from .black_scholes import black_scholes_call

    petsc_comm, mpi_comm, rank, size = _petsc_comm()
    if antithetic:
        n_paths = max((n_paths // 2) * 2, 2)
        n_pairs_total = n_paths // 2
        local_pairs = _split_paths(n_pairs_total, size, rank)
    else:
        local_pairs = 0
    local_seed = seed + 10007 * rank
    rng = np.random.default_rng(local_seed)

    with _petsc_log_event("MC_cv_paths"):
        if antithetic:
            paths = simulate_gbm_paths_antithetic(s0, r, sigma, t, n_steps, local_pairs, rng)
        else:
            local_n = _split_paths(n_paths, size, rank)
            paths = simulate_gbm_paths(s0, r, sigma, t, n_steps, local_n, rng)

    with _petsc_log_event("MC_cv_local"):
        st = paths[:, -1]
        as_mean = np.mean(paths[:, 1:], axis=1)
        disc = np.exp(-r * t)
        X = disc * np.maximum(st - strike, 0.0)
        Y = disc * np.maximum(as_mean - strike, 0.0)

        sx = float(np.sum(X))
        sy = float(np.sum(Y))
        sxx = float(np.sum(X * X))
        syy = float(np.sum(Y * Y))
        sxy = float(np.sum(X * Y))
        nloc = int(X.size)

    with _petsc_log_event("MC_cv_reduce_xy"):
        if mpi_comm is None:
            total_x = sx
            total_y = sy
            total_xx = sxx
            total_yy = syy
            total_xy = sxy
            total_n = nloc
        else:
            total_x = mpi_comm.allreduce(sx, op=MPI.SUM)
            total_y = mpi_comm.allreduce(sy, op=MPI.SUM)
            total_xx = mpi_comm.allreduce(sxx, op=MPI.SUM)
            total_yy = mpi_comm.allreduce(syy, op=MPI.SUM)
            total_xy = mpi_comm.allreduce(sxy, op=MPI.SUM)
            total_n = mpi_comm.allreduce(nloc, op=MPI.SUM)

    mu_x = total_x / total_n
    mu_y = total_y / total_n
    var_x = max(total_xx / total_n - mu_x * mu_x, 0.0)
    cov_xy = total_xy / total_n - mu_x * mu_y
    ex_bs = black_scholes_call(s0, strike, r, sigma, t)
    if var_x > 1e-24:
        c = -cov_xy / var_x
    else:
        c = 0.0

    with _petsc_log_event("MC_cv_Z"):
        Z = Y + c * (X - ex_bs)
        sz = float(np.sum(Z))
        szz = float(np.sum(Z * Z))

    with _petsc_log_event("MC_cv_reduce_Z"):
        if mpi_comm is None:
            tz, tzz = sz, szz
        else:
            tz = mpi_comm.allreduce(sz, op=MPI.SUM)
            tzz = mpi_comm.allreduce(szz, op=MPI.SUM)

    mean_z = tz / total_n
    second_moment_z = tzz / total_n
    var_z = max(second_moment_z - mean_z * mean_z, 0.0)
    std_cv = np.sqrt(var_z / total_n)
    ci_half = 1.96 * std_cv

    var_y = max(total_yy / total_n - mu_y * mu_y, 0.0)
    std_naive = np.sqrt(var_y / total_n)

    return {
        "price_cv": mean_z,
        "price_naive": mu_y,
        "std_error_cv": std_cv,
        "std_error_naive": std_naive,
        "ci_low": mean_z - ci_half,
        "ci_high": mean_z + ci_half,
        "cv_coef_c": c,
        "european_bs": ex_bs,
        "n_paths": total_n,
        "rank": rank,
        "size": size,
        "backend": "petsc" if petsc_comm is not None else "serial",
    }


# Backward-compatible alias for earlier scripts/imports.
mpi_monte_carlo_price = petsc_monte_carlo_price
