# HW4 Q3: biharmonic with RHS f = dT/dx, T = (1-y) + A*cos(pi*x), A=0.1
# Output VTK for ParaView (temperature, vorticity, streamfunction, velocity).
import argparse
import os
import sys

import petsc4py

petsc4py.init(sys.argv)
from firedrake import *  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-base", type=int, default=2, help="Coarse UnitSquareMesh(n,n) before refinement")
    parser.add_argument("--levels", type=int, default=8, help="MeshHierarchy refinement levels")
    parser.add_argument("--A", type=float, default=0.1, help="Amplitude in T = (1-y) + A*cos(pi*x)")
    args = parser.parse_args()

    base_mesh = UnitSquareMesh(args.n_base, args.n_base, quadrilateral=True)
    hierarchy = MeshHierarchy(base_mesh, args.levels)
    mesh = hierarchy[-1]
    V = FunctionSpace(mesh, "Lagrange", 1)
    ME = MixedFunctionSpace([V, V], name=["vorticity", "streamfunction"])
    Vv = VectorFunctionSpace(mesh, "Lagrange", 1)

    x, y = SpatialCoordinate(mesh)
    T = (1.0 - y) + Constant(args.A) * cos(pi * x)
    f_rhs = T.dx(0)

    omega_t, psi_t = TestFunctions(ME)
    u = Function(ME)
    u.subfunctions[0].rename("vorticity")
    u.subfunctions[1].rename("streamfunction")
    omega, psi = split(u)

    Fomega = inner(grad(omega_t), grad(omega)) * dx - omega_t * f_rhs * dx
    Fpsi = inner(grad(psi_t), grad(psi)) * dx - psi_t * omega * dx
    F = Fomega + Fpsi

    u.subfunctions[0].interpolate(0.0)
    u.subfunctions[1].interpolate(0.0)
    bcs = [
        DirichletBC(ME.sub(0), 0.0, "on_boundary"),
        DirichletBC(ME.sub(1), 0.0, "on_boundary"),
    ]

    params = {
        "snes_type": "ksponly",
        "snes_monitor": None,
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }
    problem = NonlinearVariationalProblem(F, u, bcs=bcs)
    solver = NonlinearVariationalSolver(problem, solver_parameters=params)
    solver.solve()

    T_vis = Function(V, name="temperature")
    T_vis.project(T)
    v = Function(Vv, name="velocity")
    v.interpolate(curl(psi))

    outdir = os.environ.get("HW4_RESULT_DIR", "result/q3")
    os.makedirs(outdir, exist_ok=True)
    outfile = VTKFile(os.path.join(outdir, "biharm_T_rhs.pvd"))
    outfile.write(T_vis, u.subfunctions[0], u.subfunctions[1], v)
    print(f"Wrote {os.path.join(outdir, 'biharm_T_rhs.pvd')}")


if __name__ == "__main__":
    main()
