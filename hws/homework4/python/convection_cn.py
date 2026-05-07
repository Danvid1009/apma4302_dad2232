# Extra credit (HW4 handout item 3): Crank–Nicolson in time, Q1 FE in space, plain Firedrake.
# One coupled SNES per step for (T, omega, psi) with CN on the temperature equation and
# instantaneous streamfunction–vorticity elliptic blocks at the new time level.
#
#   (T^{n+1}-T^n)/dt + v^{n+1/2}·∇T^{n+1/2} = (1/Ra)∇²T^{n+1/2},  T^{n+1/2}=(T^n+T^{n+1})/2,
#   v^{n+1/2} from ψ^{n+1/2}=(ψ^n+ψ^{n+1})/2;  -∇²ω^{n+1} = ∂T^{n+1}/∂x,  -∇²ψ^{n+1} = ω^{n+1}.
#
# Run inside the same Firedrake environment as other HW4 scripts (no firedrake-ts).
import argparse
import csv
import os
import sys

import petsc4py

petsc4py.init(sys.argv)
from firedrake import *  # noqa: E402


def solve_initial_omega_psi(V, T_expr):
    ME = MixedFunctionSpace([V, V])
    w = Function(ME)
    o, p = split(w)
    ot, pt = TestFunctions(ME)
    F = inner(grad(o), grad(ot)) * dx - inner(T_expr.dx(0), ot) * dx
    F += inner(grad(p), grad(pt)) * dx - inner(o, pt) * dx
    bcs = [
        DirichletBC(ME.sub(0), Constant(0.0), "on_boundary"),
        DirichletBC(ME.sub(1), Constant(0.0), "on_boundary"),
    ]
    params = {
        "snes_type": "ksponly",
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }
    prob = NonlinearVariationalProblem(F, w, bcs=bcs)
    NonlinearVariationalSolver(prob, solver_parameters=params).solve()
    return w.subfunctions[0], w.subfunctions[1]


def nusselt_number(mesh, T, ds_top_id=3, ds_bot_id=1):
    n = FacetNormal(mesh)
    top_flux = assemble(inner(grad(T), n) * ds(ds_top_id))
    bot_flux = assemble(inner(grad(T), n) * ds(ds_bot_id))
    if abs(float(bot_flux)) < 1e-30:
        return float("nan")
    return -float(top_flux) / float(bot_flux)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ra", type=float, default=1.0e4)
    parser.add_argument("--n", type=int, default=32)
    parser.add_argument("--t-max", type=float, default=500.0)
    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--A", type=float, default=0.1)
    parser.add_argument("--output-dir", type=str, default="result/bonus_cn")
    parser.add_argument("--ds-top", type=int, default=3)
    parser.add_argument("--ds-bot", type=int, default=1)
    parser.add_argument("--nu-every", type=int, default=10)
    args = parser.parse_args()

    mesh = UnitSquareMesh(args.n, args.n, quadrilateral=True)
    V = FunctionSpace(mesh, "CG", 1)
    ME = MixedFunctionSpace([V, V, V])
    u = Function(ME)
    u_n = Function(ME)

    x, y = SpatialCoordinate(mesh)
    T_ic = (1.0 - y) + Constant(args.A) * cos(pi * x)
    u.subfunctions[0].project(T_ic)
    o0, p0 = solve_initial_omega_psi(V, u.subfunctions[0])
    u.subfunctions[1].assign(o0)
    u.subfunctions[2].assign(p0)
    u_n.assign(u)

    T, o, p = split(u)
    q, chi, phi = TestFunctions(ME)
    tn = u_n.subfunctions[0]
    pn = u_n.subfunctions[2]
    T_mid = 0.5 * (tn + T)
    psi_mid = 0.5 * (pn + p)
    v_mid = as_vector((psi_mid.dx(1), -psi_mid.dx(0)))
    Ra = Constant(args.ra)
    inv_dt = Constant(1.0 / args.dt)

    F_T = inv_dt * inner(T - tn, q) * dx
    F_T += inner(inner(v_mid, grad(T_mid)), q) * dx + (1.0 / Ra) * inner(grad(T_mid), grad(q)) * dx
    F_o = inner(grad(o), grad(chi)) * dx - inner(T.dx(0), chi) * dx
    F_p = inner(grad(p), grad(phi)) * dx - inner(o, phi) * dx
    F = F_T + F_o + F_p

    bcs = [
        DirichletBC(ME.sub(0), Constant(1.0), "bottom"),
        DirichletBC(ME.sub(0), Constant(0.0), "top"),
        DirichletBC(ME.sub(1), Constant(0.0), "on_boundary"),
        DirichletBC(ME.sub(2), Constant(0.0), "on_boundary"),
    ]

    snes_params = {
        "snes_type": "newtonls",
        "snes_monitor": None,
        "snes_rtol": 1e-8,
        "snes_atol": 1e-10,
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }
    problem = NonlinearVariationalProblem(F, u, bcs=bcs)
    solver = NonlinearVariationalSolver(problem, solver_parameters=snes_params)

    os.makedirs(args.output_dir, exist_ok=True)
    csv_path = os.path.join(args.output_dir, "nu_history_cn.csv")
    nu_rows = []
    t = 0.0
    Th = u.subfunctions[0]
    nu0 = nusselt_number(mesh, Th, ds_top_id=args.ds_top, ds_bot_id=args.ds_bot)
    nu_rows.append((t, nu0))
    print(f"t={t:.6g}  Nu={nu0:.6f}")

    n_steps = int(round(args.t_max / args.dt))
    for k in range(1, n_steps + 1):
        u.assign(u_n)
        solver.solve()
        u_n.assign(u)
        t = k * args.dt
        if k % args.nu_every == 0 or k == n_steps:
            nu = nusselt_number(mesh, u.subfunctions[0], ds_top_id=args.ds_top, ds_bot_id=args.ds_bot)
            nu_rows.append((t, nu))
            print(f"t={t:.6g}  Nu={nu:.6f}")

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time", "Nu"])
        w.writerows(nu_rows)
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
