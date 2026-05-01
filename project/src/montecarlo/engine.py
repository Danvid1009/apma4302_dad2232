from typing import Callable, Optional

import numpy as np

from .gbm import simulate_gbm_paths

try:
    from petsc4py import PETSc
except ImportError:
    PETSc = None

try:
    from mpi4py import MPI
except ImportError:
    MPI = None


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
):
    petsc_comm, mpi_comm, rank, size = _petsc_comm()
    local_n = _split_paths(n_paths, size, rank)
    local_seed = seed + 10007 * rank
    rng = np.random.default_rng(local_seed)

    paths = simulate_gbm_paths(s0, r, sigma, t, n_steps, local_n, rng)
    if barrier is None:
        payoff = payoff_builder(strike)(paths)
    else:
        payoff = payoff_builder(strike, barrier)(paths)
    discounted = np.exp(-r * t) * payoff

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


# Backward-compatible alias for earlier scripts/imports.
mpi_monte_carlo_price = petsc_monte_carlo_price
