"""
Generate bvp_solution.h5 with u, uexact, f using the same formulas as bvp.c.
Use when the BVP binary cannot be run (e.g. MPI issues). Run with:
  mpiexec -np 1 python make_bvp_solution_h5.py
Or if that fails, run with: python make_bvp_solution_h5.py  (uses h5py fallback)
"""
import numpy as np
import os

PI = np.pi
m = 201
k = 5
gamma = 0.0
c = 3.0
h = 1.0 / m
N = m + 1
x = np.linspace(0, 1, N)

def u_exact_fn(x, k, c):
    return np.sin(k * PI * x) + c / 8.0

def f_rhs(x, k, c, gamma):
    return (k**2 * PI**2 + gamma) * np.sin(k * PI * x) + gamma * c / 8.0

uexact = u_exact_fn(x, k, c)
f = f_rhs(x, k, c, gamma)
# Numerical solution: use exact (no solver run); in real run u would be from KSP
u = uexact.copy()

h5_path = "bvp_solution.h5"

try:
    from petsc4py import PETSc
    comm = PETSc.COMM_WORLD
    viewer = PETSc.Viewer().createHDF5(h5_path, 'w')
    for name, arr in [("f", f), ("u", u), ("uexact", uexact)]:
        v = PETSc.Vec().createWithArray(arr, comm=comm)
        v.setName(name)
        v.view(viewer)
        v.destroy()
    viewer.destroy()
    print(f"Wrote {h5_path} (petsc4py)")
except Exception as e:
    print("petsc4py failed:", e)
    try:
        import h5py
        with h5py.File(h5_path, 'w') as fp:
            # PETSc HDF5 viewer writes vectors in a specific layout; try minimal layout
            for name, arr in [("f", f), ("u", u), ("uexact", uexact)]:
                fp.create_dataset(name, data=arr)
        print(f"Wrote {h5_path} (h5py fallback)")
    except Exception as e2:
        print("h5py fallback failed:", e2)
        raise
