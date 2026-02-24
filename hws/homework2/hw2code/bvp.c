/*
  APMA 4302 - Homework 2, Problem 3(b)
  Modified from tri.c (p4pdes, Chapter 2).

  BVP: -u''(x) + gamma*u(x) = f(x) on [0,1], Dirichlet BCs.
  Manufactured solution: u(x) = sin(k*pi*x) + c/8  (from 3(a)).

  Run: mpiexec -np P ./bvp -options_file options_file
  Requirements:
  - PETSc options handling for command-line/options_file
  - Parallel assembly of matrix A and RHS f; BCs via MatZeroRowsColumns (preserve symmetry)
  - Compute and print relative error
  - HDF5 output: u, f, uexact (PETSc built with --download-hdf5)
*/

#include <petsc.h>
#include <petscviewerhdf5.h>
#include <math.h>

#define PI 3.14159265358979323846

/* Exact solution u(x) = sin(k*pi*x) + c/8 */
static inline PetscReal u_exact_fn(PetscReal x, PetscInt k, PetscReal c)
{
  return (PetscReal)sin((double)(k * PI * x)) + c / 8.0;
}

/* RHS f(x) = (k^2*pi^2 + gamma)*sin(k*pi*x) + gamma*c/8  so that -u''+gamma*u = f */
static inline PetscReal f_rhs(PetscReal x, PetscInt k, PetscReal c, PetscReal gamma)
{
  return (k * k * PI * PI + gamma) * (PetscReal)sin((double)(k * PI * x)) + gamma * c / 8.0;
}

int main(int argc, char **argv)
{
  MPI_Comm       comm = PETSC_COMM_WORLD;
  PetscInt       m = 40, k = 1, i, N, rStart, rEnd, cols[3];
  PetscReal      gamma = 0.0, c = 1.0, h, val, norm_err, norm_exact, relerr;
  PetscScalar    v[3];
  Mat            A;
  Vec            f, u, uexact;
  KSP            ksp;
  PetscViewer    viewer;
  PetscInt       bc_rows[2] = {0, 0};  /* rows to zero for Dirichlet: 0 and m */

  /* -------------------------------------------------------------------------
   * Step 1: Initialize PETSc and use options handling so the code can be run
   *         with:  mpiexec -np P ./bvp -options_file options_file
   * ------------------------------------------------------------------------- */
  PetscCall(PetscInitialize(&argc, &argv, NULL, "BVP -u'' + gamma*u = f on [0,1], Dirichlet BCs.\n\n"));

  PetscOptionsBegin(comm, NULL, "bvp", NULL);
  PetscCall(PetscOptionsInt("-bvp_m", "Number of intervals (grid points = m+1)", NULL, m, &m, NULL));
  PetscCall(PetscOptionsReal("-bvp_gamma", "Coefficient gamma in -u''+gamma*u=f", NULL, gamma, &gamma, NULL));
  PetscCall(PetscOptionsInt("-bvp_k", "Wavenumber k in manufactured solution", NULL, k, &k, NULL));
  PetscCall(PetscOptionsReal("-bvp_c", "Constant c in manufactured solution", NULL, c, &c, NULL));
  PetscOptionsEnd();

  PetscCheck(m > 0, comm, PETSC_ERR_USER_INPUT, "bvp_m must be > 0");
  N         = m + 1;
  h         = 1.0 / (PetscReal)m;
  bc_rows[1] = m;

  /* Create parallel matrix A (N x N) and vectors f, u, uexact (length N) */
  PetscCall(MatCreate(comm, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));
  PetscCall(MatGetOwnershipRange(A, &rStart, &rEnd));

  PetscCall(VecCreate(comm, &f));
  PetscCall(VecSetSizes(f, PETSC_DECIDE, N));
  PetscCall(VecSetFromOptions(f));
  PetscCall(VecDuplicate(f, &u));
  PetscCall(VecDuplicate(f, &uexact));

  /* -------------------------------------------------------------------------
   * Step 2: Assemble the matrix A and right-hand side vector f in parallel.
   *         Each process owns rows [rStart, rEnd). Stencil for -u'' + gamma*u:
   *         (1/h^2)*(-1, 2, -1) on off-diagonals/diagonal plus gamma on diagonal.
   *         Boundary rows 0 and m are assembled here; we will zero them and
   *         set BCs via MatZeroRowsColumns to preserve symmetry.
   * ------------------------------------------------------------------------- */
  for (i = rStart; i < rEnd; i++) {
    PetscReal xi = (PetscReal)i * h;
    PetscCall(VecSetValue(f, i, (PetscScalar)f_rhs(xi, k, c, gamma), INSERT_VALUES));
    PetscCall(VecSetValue(uexact, i, (PetscScalar)u_exact_fn(xi, k, c), INSERT_VALUES));

    if (i == 0) {
      /* First row: A(0,0) = 2/h^2 + gamma, A(0,1) = -1/h^2 */
      cols[0] = 0; cols[1] = 1;
      val = 2.0 / (h * h) + gamma;
      PetscCall(MatSetValue(A, 0, 0, (PetscScalar)val, INSERT_VALUES));
      PetscCall(MatSetValue(A, 0, 1, (PetscScalar)(-1.0 / (h * h)), INSERT_VALUES));
    } else if (i == (PetscInt)m) {
      /* Last row: A(m,m-1) = -1/h^2, A(m,m) = 2/h^2 + gamma */
      cols[0] = m - 1; cols[1] = m;
      val = 2.0 / (h * h) + gamma;
      PetscCall(MatSetValue(A, m, m - 1, (PetscScalar)(-1.0 / (h * h)), INSERT_VALUES));
      PetscCall(MatSetValue(A, m, m, (PetscScalar)val, INSERT_VALUES));
    } else {
      /* Interior: stencil [-1/h^2, 2/h^2+gamma, -1/h^2] */
      cols[0] = i - 1; cols[1] = i; cols[2] = i + 1;
      v[0] = -1.0 / (h * h);
      v[1] = 2.0 / (h * h) + gamma;
      v[2] = -1.0 / (h * h);
      PetscCall(MatSetValues(A, 1, &i, 3, cols, v, INSERT_VALUES));
    }
  }
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscCall(VecAssemblyBegin(f));
  PetscCall(VecAssemblyEnd(f));
  PetscCall(VecAssemblyBegin(uexact));
  PetscCall(VecAssemblyEnd(uexact));

  /* -------------------------------------------------------------------------
   * Step 3: Enforce Dirichlet BCs and preserve symmetry by calling
   *         MatZeroRowsColumns after assembling A, uexact, and f.
   *         This zeros rows and columns 0 and m in A, sets diagonal to 1,
   *         and sets f[0] = uexact[0], f[m] = uexact[m] (so the equations
   *         become u[0]=uexact[0], u[m]=uexact[m]).
   * ------------------------------------------------------------------------- */
  PetscCall(MatZeroRowsColumns(A, 2, bc_rows, 1.0, uexact, f));

  /* Solve A*u = f */
  PetscCall(KSPCreate(comm, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSolve(ksp, f, u));

  /* -------------------------------------------------------------------------
   * Step 4: Compute the relative error in the solution and print it to the
   *         console:  ||u - uexact|| / ||uexact||  (2-norm).
   * ------------------------------------------------------------------------- */
  PetscCall(VecAXPY(u, -1.0, uexact));   /* u := u - uexact (error vector) */
  PetscCall(VecNorm(u, NORM_2, &norm_err));
  PetscCall(VecAXPY(u, 1.0, uexact));    /* restore u = solution for HDF5 */
  PetscCall(VecNorm(uexact, NORM_2, &norm_exact));
  relerr = (norm_exact > 1e-30) ? (norm_err / norm_exact) : norm_err;
  PetscCall(PetscPrintf(comm, "Relative error (2-norm): %.6e\n", (double)relerr));

  /* -------------------------------------------------------------------------
   * Step 5: Output solution, RHS, and exact solution to HDF5 (requires
   *         PETSc built with --download-hdf5). Named datasets: "u", "f", "uexact"
   *         for use with the provided plot_bvp.py script.
   * ------------------------------------------------------------------------- */
#if defined(PETSC_HAVE_HDF5)
  PetscCall(PetscViewerHDF5Open(comm, "bvp_solution.h5", FILE_MODE_WRITE, &viewer));
  PetscCall(PetscObjectSetName((PetscObject)uexact, "uexact"));
  PetscCall(PetscObjectSetName((PetscObject)f, "f"));
  PetscCall(PetscObjectSetName((PetscObject)u, "u"));
  PetscCall(VecView(f, viewer));
  PetscCall(VecView(u, viewer));
  PetscCall(VecView(uexact, viewer));
  PetscCall(PetscViewerDestroy(&viewer));
#else
  PetscCall(PetscPrintf(comm, "HDF5 not available; skipping bvp_solution.h5\n"));
#endif

  PetscCall(KSPDestroy(&ksp));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&f));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&uexact));
  PetscCall(PetscFinalize());
  return 0;
}
