# HW4 Q4: streamfunction–vorticity thermal convection as a DAE with firedrake-ts.
# Matches handout Eqs. (4)–(6): T_t + v·∇T = (1/Ra)∇²T,  -∇²ω = ∂T/∂x,  -∇²ψ = ω,
# v = (∂ψ/∂y, -∂ψ/∂x). BCs: T=1 on y=0, T=0 on y=1, natural (no-flux) sides ⇒ ∂T/∂x=0;
# ω=0, ψ=0 on ∂Ω. IC: T = (1-y) + A cos(πx); ω,ψ from coupled elliptic solve at t=0.
# Time: BDF-2 (PETSc TS); each step: monolithic LU with MUMPS (ksp preonly / pc lu).
# Nu = -(∫ ∂T/∂n on top)/(∫ ∂T/∂n on bottom) via FacetNormal + ds markers (--ds-top/--ds-bot).
# Assignment batching: hw4/scripts/run_q4_assignment.sh (parts a–c) or run_q4.sh for one-off runs.
import argparse
import csv
import os
import sys

import petsc4py

petsc4py.init(sys.argv)
from firedrake import *  # noqa: E402

import firedrake_ts  # noqa: E402


def solve_initial_omega_psi(V, T_expr):
    """Given T in V, solve -∇²ω = ∂T/∂x, -∇²ψ = ω with homogeneous Dirichlet BCs."""
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
    sol = NonlinearVariationalSolver(prob, solver_parameters=params)
    sol.solve()
    return w.subfunctions[0], w.subfunctions[1]


def nusselt_number(mesh, T, ds_top_id=3, ds_bot_id=1):
    """Nu = -(∫ ∂T/∂y on top)/(∫ ∂T/∂y on bottom) using facet normals and ds markers."""
    n = FacetNormal(mesh)
    top_flux = assemble(inner(grad(T), n) * ds(ds_top_id))
    bot_flux = assemble(inner(grad(T), n) * ds(ds_bot_id))
    if abs(float(bot_flux)) < 1e-30:
        return float("nan")
    return -float(top_flux) / float(bot_flux)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ra", type=float, default=1.0e2, help="Rayleigh number")
    parser.add_argument("--n", type=int, default=64, help="Mesh cells per direction (NxN Q1)")
    parser.add_argument("--t-max", type=float, default=1.0e5, help="End time for TS")
    parser.add_argument("--dt", type=float, default=0.1, help="Fixed time step")
    parser.add_argument("--A", type=float, default=0.1, help="Initial T perturbation amplitude")
    parser.add_argument("--vtk-every", type=float, default=0.0, help="Write VTK every this many time units (0=off)")
    parser.add_argument("--nu-every", type=int, default=10, help="Print / log Nu every this many TS steps")
    parser.add_argument("--output-dir", type=str, default="result/q4")
    parser.add_argument(
        "--ds-top",
        type=int,
        default=3,
        help="Facet marker for y=1 (UnitSquareMesh default; use with Dirichlet BCs and Nu)",
    )
    parser.add_argument(
        "--ds-bot",
        type=int,
        default=1,
        help="Facet marker for y=0 (UnitSquareMesh default)",
    )
    args = parser.parse_args()

    mesh = UnitSquareMesh(args.n, args.n, quadrilateral=True)
    V = FunctionSpace(mesh, "CG", 1)
    ME = MixedFunctionSpace([V, V, V])
    u = Function(ME)
    u.subfunctions[0].rename("temperature")
    u.subfunctions[1].rename("vorticity")
    u.subfunctions[2].rename("streamfunction")
    u_dot = Function(ME)

    x, y = SpatialCoordinate(mesh)
    T_ic = (1.0 - y) + Constant(args.A) * cos(pi * x)
    u.subfunctions[0].project(T_ic)
    o0, p0 = solve_initial_omega_psi(V, u.subfunctions[0])
    u.subfunctions[1].assign(o0)
    u.subfunctions[2].assign(p0)

    u_dot.assign(0.0)

    T, o, p = split(u)
    Td, od, pd = split(u_dot)
    q, chi, phi = TestFunctions(ME)

    Ra = Constant(args.ra)
    v = as_vector((p.dx(1), -p.dx(0)))

    F_T = inner(Td, q) * dx + inner(inner(v, grad(T)), q) * dx + (1.0 / Ra) * inner(grad(T), grad(q)) * dx
    F_o = inner(grad(o), grad(chi)) * dx - inner(T.dx(0), chi) * dx
    F_p = inner(grad(p), grad(phi)) * dx - inner(o, phi) * dx
    F = F_T + F_o + F_p

    # Use facet IDs, not "bottom"/"top" strings — those are invalid on non-extruded
    # UnitSquareMesh in current Firedrake (raises ValueError in boundary_nodes).
    bcs = [
        DirichletBC(ME.sub(0), Constant(1.0), args.ds_bot),
        DirichletBC(ME.sub(0), Constant(0.0), args.ds_top),
        DirichletBC(ME.sub(1), Constant(0.0), "on_boundary"),
        DirichletBC(ME.sub(2), Constant(0.0), "on_boundary"),
    ]

    t_init = 0.0
    t_max = args.t_max
    params = {
        "ts_type": "bdf",
        "ts_bdf_order": 2,
        "ts_dt": args.dt,
        "ts_monitor": None,
        "ts_rtol": 1e-6,
        "ts_atol": 1e-10,
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "ts_max_time": t_max,
        "ts_adapt_dt_min": 1.0e-9,
        "ts_exact_final_time": "matchstep",
    }

    os.makedirs(args.output_dir, exist_ok=True)
    csv_path = os.path.join(args.output_dir, "nu_history.csv")
    vtk_every = args.vtk_every
    outfile = None
    if vtk_every > 0:
        outfile = VTKFile(os.path.join(args.output_dir, "convection.pvd"))

    nu_rows = []
    Th0 = u.subfunctions[0]
    nu0 = nusselt_number(mesh, Th0, ds_top_id=args.ds_top, ds_bot_id=args.ds_bot)
    nu_rows.append((0.0, nu0))
    print(f"t=0  Nu={nu0:.6f}")

    step_count = [0]
    last_vtk_t = [-1.0e300]
    Vv = VectorFunctionSpace(mesh, "CG", 1)
    vel = Function(Vv, name="velocity")

    def monitor(ts, step, t, xvec):
        step_count[0] += 1
        Th = u.subfunctions[0]
        log_nu = (step_count[0] % args.nu_every == 0) or (t >= t_max - 1e-9)
        if log_nu and t > 1e-12:
            nu = nusselt_number(mesh, Th, ds_top_id=args.ds_top, ds_bot_id=args.ds_bot)
            nu_rows.append((float(t), nu))
            print(f"t={float(t):.6g}  Nu={nu:.6f}")
        if outfile is not None and (t - last_vtk_t[0] >= vtk_every - 1e-12):
            last_vtk_t[0] = float(t)
            psi_h = u.subfunctions[2]
            vel.interpolate(as_vector((psi_h.dx(1), -psi_h.dx(0))))
            outfile.write(u.subfunctions[0], u.subfunctions[1], u.subfunctions[2], vel, time=t)

    if outfile is not None:
        psi_h = u.subfunctions[2]
        vel.interpolate(as_vector((psi_h.dx(1), -psi_h.dx(0))))
        outfile.write(u.subfunctions[0], u.subfunctions[1], u.subfunctions[2], vel, time=0.0)
        last_vtk_t[0] = 0.0

    problem = firedrake_ts.DAEProblem(F, u, u_dot, (t_init, t_max), bcs=bcs)
    solver = firedrake_ts.DAESolver(problem, solver_parameters=params, monitor_callback=monitor)
    solver.solve()

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time", "Nu"])
        w.writerows(nu_rows)
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
