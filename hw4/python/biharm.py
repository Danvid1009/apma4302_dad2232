# Solve the biharmonic equation as a coupled system of two Poisson equations.
# Reproduces the manufactured solution from the FD code in C, using Firedrake.
# Run: python biharm.py --preset direct|split_direct|split_mg
import argparse
import os
import sys

import petsc4py

petsc4py.init(sys.argv)
from firedrake import *  # noqa: E402


def solver_parameters_for_preset(preset: str) -> dict:
    """Match PETSc option files in hw4/c/ (monolithic LU vs fieldsplit direct vs fieldsplit MG)."""
    common = {
        "snes_type": "ksponly",
        "snes_monitor": None,
        "ksp_monitor": None,
        "ksp_rtol": 1.0e-6,
        "ksp_atol": 1.0e-10,
    }
    if preset == "direct":
        p = {
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        }
    elif preset == "split_direct":
        p = {
            "ksp_type": "fgmres",
            "pc_type": "fieldsplit",
            "pc_fieldsplit_type": "multiplicative",
            "fieldsplit_0_ksp_type": "preonly",
            "fieldsplit_0_pc_type": "lu",
            "fieldsplit_0_pc_factor_mat_solver_type": "mumps",
            "fieldsplit_1_ksp_type": "preonly",
            "fieldsplit_1_pc_type": "lu",
            "fieldsplit_1_pc_factor_mat_solver_type": "mumps",
        }
    elif preset == "split_mg":
        p = {
            "ksp_type": "fgmres",
            "pc_type": "fieldsplit",
            "pc_fieldsplit_type": "multiplicative",
            "fieldsplit_0_ksp_type": "preonly",
            "fieldsplit_0_pc_type": "mg",
            "fieldsplit_0_pc_mg_levels": 8,
            "fieldsplit_0_pc_mg_galerkin": None,
            "fieldsplit_1_ksp_type": "preonly",
            "fieldsplit_1_pc_type": "mg",
            "fieldsplit_1_pc_mg_levels": 8,
            "fieldsplit_1_pc_mg_galerkin": None,
        }
    else:
        raise ValueError(preset)
    p.update(common)
    return p


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--preset",
        choices=("direct", "split_direct", "split_mg"),
        default="direct",
        help="Solver stack aligned with c/options_file_*",
    )
    args = parser.parse_args()

    N = 2
    levels = 8
    Nfine = N * 2**levels
    base_mesh = UnitSquareMesh(N, N, quadrilateral=True)
    hierarchy = MeshHierarchy(base_mesh, levels)
    mesh = hierarchy[-1]
    V = FunctionSpace(mesh, "Lagrange", 1)
    ME = MixedFunctionSpace([V, V], name=["vorticity", "streamfunction"])
    Vv = VectorFunctionSpace(mesh, "Lagrange", 1)

    omega_t, psi_t = TestFunctions(ME)
    u = Function(ME)
    u.subfunctions[0].rename("vorticity")
    u.subfunctions[1].rename("streamfunction")
    omega, psi = split(u)

    f = Function(V, name="rhs")
    x, z = SpatialCoordinate(mesh)
    cx = x**3 * (1.0 - x) ** 3
    cz = z**3 * (1.0 - z) ** 3
    ddcx = 6.0 * x * (1.0 - x) * (1.0 - 5.0 * x + 5.0 * x * x)
    ddcz = 6.0 * z * (1.0 - z) * (1.0 - 5.0 * z + 5.0 * z * z)
    d4cx = -72.0 * (1.0 - 5.0 * x + 5.0 * x * x)
    d4cz = -72.0 * (1.0 - 5.0 * z + 5.0 * z * z)
    f.interpolate(d4cx * cz + 2.0 * ddcx * ddcz + cx * d4cz)

    u_true = Function(ME, name="u_true")
    u_true.subfunctions[0].interpolate(-ddcx * cz - cx * ddcz)
    u_true.subfunctions[1].interpolate(cx * cz)

    Fomega = inner(grad(omega_t), grad(omega)) * dx - omega_t * f * dx
    Fpsi = inner(grad(psi_t), grad(psi)) * dx - psi_t * omega * dx
    F = Fomega + Fpsi

    u.subfunctions[0].interpolate(0.0)
    u.subfunctions[1].interpolate(0.0)
    bcs = [
        DirichletBC(ME.sub(0), 0.0, "on_boundary"),
        DirichletBC(ME.sub(1), 0.0, "on_boundary"),
    ]

    params = solver_parameters_for_preset(args.preset)
    problem = NonlinearVariationalProblem(F, u, bcs=bcs)
    solver = NonlinearVariationalSolver(problem, solver_parameters=params)
    solver.solve()

    v = Function(Vv, name="Velocity")
    v.interpolate(curl(psi))
    outdir = os.environ.get("HW4_RESULT_DIR", "result")
    os.makedirs(outdir, exist_ok=True)
    outfile = VTKFile(os.path.join(outdir, "biharm.pvd"))
    outfile.write(f, u.subfunctions[0], u.subfunctions[1], v)

    abs_error_o = assemble(
        sqrt(inner(u.subfunctions[0] - u_true.subfunctions[0], u.subfunctions[0] - u_true.subfunctions[0])) * dx
    )
    abs_error_p = assemble(
        sqrt(inner(u.subfunctions[1] - u_true.subfunctions[1], u.subfunctions[1] - u_true.subfunctions[1])) * dx
    )
    rel_error_o = abs_error_o / assemble(sqrt(inner(u_true.subfunctions[0], u_true.subfunctions[0])) * dx)
    rel_error_p = abs_error_p / assemble(sqrt(inner(u_true.subfunctions[1], u_true.subfunctions[1])) * dx)
    print(f"\nworking on ~{Nfine}x{Nfine} effective resolution, preset={args.preset}, ksp={params['ksp_type']}")
    print(f"L2 error (vorticity) -- abs: {float(abs_error_o):.3e}, rel: {float(rel_error_o):.3e}")
    print(f"L2 error (streamfun) -- abs: {float(abs_error_p):.3e}, rel: {float(rel_error_p):.3e}")


if __name__ == "__main__":
    main()
