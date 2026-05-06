# HW4 Q3: temperature-driven biharmonic — RHS f = dT/dx, no manufactured closed form.
# VTK: biharm_temp.pvd (temperature, vorticity, streamfunction, velocity).
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
    args = parser.parse_args()

    base_mesh = UnitSquareMesh(args.n_base, args.n_base, quadrilateral=True)
    hierarchy = MeshHierarchy(base_mesh, args.levels)
    mesh = hierarchy[-1]
    V = FunctionSpace(mesh, "Lagrange", 1)
    ME = MixedFunctionSpace([V, V], name=["vorticity", "streamfunction"])
    Vv = VectorFunctionSpace(mesh, "Lagrange", 1)

    x, z = SpatialCoordinate(mesh)
    A = Constant(0.1)
    T = Function(V, name="temperature")
    T.interpolate((1.0 - z) + A * cos(pi * x))

    f = Function(V, name="dT_dx")
    f.interpolate(T.dx(0))

    omega_t, psi_t = TestFunctions(ME)
    u = Function(ME)
    u.subfunctions[0].rename("vorticity")
    u.subfunctions[1].rename("streamfunction")
    omega, psi = split(u)

    Fomega = inner(grad(omega_t), grad(omega)) * dx - omega_t * f * dx
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

    v = Function(Vv, name="velocity")
    v.interpolate(curl(psi))

    outdir = os.environ.get("HW4_RESULT_DIR", "result/biharm_temp")
    os.makedirs(outdir, exist_ok=True)
    outfile = VTKFile(os.path.join(outdir, "biharm_temp.pvd"))
    outfile.write(T, u.subfunctions[0], u.subfunctions[1], v)
    print(f"Wrote {os.path.join(outdir, 'biharm_temp.pvd')}")


if __name__ == "__main__":
    main()
