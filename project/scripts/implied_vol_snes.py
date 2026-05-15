#!/usr/bin/env python3
"""
Implied volatility via PETSc SNES (scalar root of BS(σ) − market = 0).

Course tie-in: Newton / SNES weeks (nonlinear equations). Optional **analytic
1×1 Jacobian** (vega) vs PETSc finite-difference Jacobian.

Under ``mpirun -n p``, only rank 0 builds SNES on ``COMM_SELF``; implied σ is
broadcast so all ranks exit consistently (embed-friendly).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from petsc4py import PETSc

from montecarlo.black_scholes import black_scholes_call, black_scholes_call_vega


def parse_args():
    p = argparse.ArgumentParser(description="Implied vol with PETSc SNES (European call).")
    p.add_argument("--s0", type=float, default=100.0)
    p.add_argument("--strike", type=float, default=100.0)
    p.add_argument("--r", type=float, default=0.05)
    p.add_argument("--t", type=float, default=1.0)
    p.add_argument("--market-price", type=float, required=True, help="Target call price (e.g. mid).")
    p.add_argument("--sigma0", type=float, default=0.25, help="Initial guess for σ.")
    p.add_argument(
        "--analytic-jacobian",
        action="store_true",
        help="Use hand-coded 1×1 Jacobian (vega); else SNES finite differences.",
    )
    p.add_argument("--out", type=str, default="", help="Optional JSON path with solved σ and residual norm.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    comm = PETSc.COMM_WORLD
    rank = comm.getRank()
    mp = comm.tompi4py()

    sigma_star = 0.0
    iters = 0
    res_norm = 0.0

    if rank == 0:
        s0, K, r, T = args.s0, args.strike, args.r, args.t
        target = args.market_price

        def form_function(snes, X: PETSc.Vec, F: PETSc.Vec) -> None:
            sig = float(X.getArray(readonly=True)[0])
            F.zeroEntries()
            F.setValue(0, black_scholes_call(s0, K, r, sig, T) - target)
            F.assemble()

        def form_jacobian(snes, X: PETSc.Vec, J: PETSc.Mat, P: PETSc.Mat) -> None:
            sig = float(X.getArray(readonly=True)[0])
            v = black_scholes_call_vega(s0, K, r, sig, T)
            P.zeroEntries()
            P.setValue(0, 0, v)
            P.assemble()
            if J != P:
                J.zeroEntries()
                J.assemble()

        sigma = PETSc.Vec().create(comm=PETSc.COMM_SELF)
        sigma.setSizes(1)
        sigma.setFromOptions()
        sigma.set(args.sigma0)
        F = PETSc.Vec().create(comm=PETSc.COMM_SELF)
        F.setSizes(1)
        F.setFromOptions()

        snes = PETSc.SNES().create(comm=PETSc.COMM_SELF)
        snes.setOptionsPrefix("iv_")
        snes.setFunction(form_function, F)
        if args.analytic_jacobian:
            J = PETSc.Mat().createAIJ([1, 1], nnz=1, comm=PETSc.COMM_SELF)
            J.setUp()
            snes.setJacobian(form_jacobian, J, J)
        else:
            snes.setUseFD(True)
        snes.setFromOptions()
        snes.solve(None, sigma)
        sigma_star = float(sigma.getArray(readonly=True)[0])
        iters = snes.getIterationNumber()
        res_norm = float(F.norm())
        model = black_scholes_call(s0, K, r, sigma_star, T)
        print(f"implied_sigma={sigma_star:.12f} BS(model)={model:.10f} market={target:.10f}")
        print(f"snes_iters={iters} ||F||_2={res_norm:.3e} analytic_jac={args.analytic_jacobian}")

        if args.out:
            out = Path(args.out)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(
                json.dumps(
                    {
                        "implied_sigma": sigma_star,
                        "market_price": target,
                        "bs_at_sigma": model,
                        "snes_iters": iters,
                        "residual_norm": res_norm,
                        "analytic_jacobian": bool(args.analytic_jacobian),
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            print(f"wrote {out}")

    buf = np.empty(1, dtype=np.float64)
    if rank == 0:
        buf[0] = sigma_star
    mp.Bcast(buf, root=0)
    sigma_star = float(buf[0])

    comm.barrier()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
