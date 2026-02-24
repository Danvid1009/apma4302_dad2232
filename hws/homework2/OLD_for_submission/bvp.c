/*
  APMA 4302 - Homework 2, Problem 3(a)+(b)
  Modified from tri.c (p4pdes, Chapter 2).

  BVP: -u''(x) + gamma*u(x) = f(x) on [0,1], Dirichlet BCs.

  ADDITIONS FOR HW2 (beyond original tri.c):
  - 3(a): Manufactured solution u(x) and derived f(x); two forms (see -bvp_form below).
  - 3(b): PETSc options (-bvp_m, -bvp_gamma, -bvp_k, -bvp_c) so run works with
    -options_file options_file; parallel assembly (ownership [rStart,rEnd]);
    MatZeroRowsColumns for Dirichlet BCs (preserve symmetry); relative error
    printed; HDF5 output (u, f, uexact) for plot_bvp.py.

  3(a) Manufactured solutions (choose with -bvp_form):
    Form 0 (default): u(x) = sin(k*pi*x) + c*(1-1/2)^3 = sin(k*pi*x) + c/8.
      Then -u'' = k^2*pi^2*sin(k*pi*x), so
      f(x) = (k^2*pi^2 + gamma)*sin(k*pi*x) + gamma*c/8.
    Form 1: u(x) = sin(k*pi*x) + c*(1-x)^3.
      Then u'' = -k^2*pi^2*sin(k*pi*x) + 6c*(1-x), so
      f(x) = (k^2*pi^2 + gamma)*sin(k*pi*x) + gamma*c*(1-x)^3 - 6c*(1-x).

  Run: mpiexec -np P ./bvp -options_file options_file
*/

#include <petsc.h>
#include <petscviewerhdf5.h>
#include <math.h>

#define PI 3.14159265358979323846

/* [3(a)] Exact u(x) and RHS f(x). form=0: problem statement (c/8); form=1: alternative c*(1-x)^3. */
static inline PetscReal u_exact_fn(PetscReal x, PetscInt k, PetscReal c, PetscInt form)
{
  PetscReal one_minus_x = 1.0 - x;
  if (form == 0)
    return (PetscReal)sin((double)(k * PI * x)) + c / 8.0;
  return (PetscReal)sin((double)(k * PI * x)) + c * (one_minus_x * one_minus_x * one_minus_x);
}

static inline PetscReal f_rhs(PetscReal x, PetscInt k, PetscReal c, PetscReal gamma, PetscInt form)
{
  PetscReal sx = (PetscReal)sin((double)(k * PI * x));
  PetscReal one_minus_x = 1.0 - x;
  if (form == 0)
    return (k * k * PI * PI + gamma) * sx + gamma * c / 8.0;
  return (k * k * PI * PI + gamma) * sx + gamma * c * (one_minus_x * one_minus_x * one_minus_x) - 6.0 * c * one_minus_x;
}

int main(int argc, char **argv)
{
  PetscInt       m = 40, k = 1, form = 0, i, N, rStart, rEnd, cols[3];
  PetscReal      gamma = 0.0, c = 1.0, h, val, norm_err, norm_exact, relerr;
  PetscScalar    v[3];
  Mat            A;
  Vec            f, u, uexact;
  KSP            ksp;
  PetscInt       bc_rows[2] = {0, 0};  /* rows to zero for Dirichlet: 0 and m */

  /* [3(b)] Step 1: Options handling so run works with -options_file options_file. */
  PetscCall(PetscInitialize(&argc, &argv, NULL, "BVP -u''+gamma*u=f [3(a) manufactured solution].\n\n"));

  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "bvp", NULL);
  PetscCall(PetscOptionsInt("-bvp_m", "Number of intervals (grid points = m+1)", NULL, m, &m, NULL));
  PetscCall(PetscOptionsReal("-bvp_gamma", "Coefficient gamma in -u''+gamma*u=f", NULL, gamma, &gamma, NULL));
  PetscCall(PetscOptionsInt("-bvp_k", "Wavenumber k in manufactured solution", NULL, k, &k, NULL));
  PetscCall(PetscOptionsReal("-bvp_c", "Constant c in manufactured solution", NULL, c, &c, NULL));
  PetscCall(PetscOptionsInt("-bvp_form", "3(a) form: 0 = sin(k*pi*x)+c/8, 1 = sin(k*pi*x)+c*(1-x)^3", NULL, form, &form, NULL));  /* [3(a)] */
  PetscOptionsEnd();

  PetscCheck(m > 0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "bvp_m must be > 0");
  PetscCheck(form == 0 || form == 1, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "bvp_form must be 0 or 1");
  N         = m + 1;
  h         = 1.0 / (PetscReal)m;
  bc_rows[1] = m;

  /* Create parallel matrix A (N x N) and vectors f, u, uexact (length N) */
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));
  PetscCall(MatGetOwnershipRange(A, &rStart, &rEnd));

  PetscCall(VecCreate(PETSC_COMM_WORLD, &f));
  PetscCall(VecSetSizes(f, PETSC_DECIDE, N));
  PetscCall(VecSetFromOptions(f));
  PetscCall(VecDuplicate(f, &u));
  PetscCall(VecDuplicate(f, &uexact));

  /* [3(b)] Step 2: Parallel assembly. Each process owns rows [rStart, rEnd).
   *         Stencil (1/h^2)*(-1,2,-1) + gamma on diagonal; fill f and uexact from 3(a). */
  for (i = rStart; i < rEnd; i++) {
    PetscReal xi = (PetscReal)i * h;
    PetscCall(VecSetValue(f, i, (PetscScalar)f_rhs(xi, k, c, gamma, form), INSERT_VALUES));
    PetscCall(VecSetValue(uexact, i, (PetscScalar)u_exact_fn(xi, k, c, form), INSERT_VALUES));

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

  /* [3(b)] Step 3: Dirichlet BCs with symmetry preserved (MatZeroRowsColumns on rows 0 and m). */
  PetscCall(MatZeroRowsColumns(A, 2, bc_rows, 1.0, uexact, f));

  /* Solve A*u = f */
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSolve(ksp, f, u));

  /* [3(b)] Step 4: Relative error ||u - uexact||_2 / ||uexact||_2, printed to console. */
  PetscCall(VecAXPY(u, -1.0, uexact));   /* u := u - uexact (error vector) */
  PetscCall(VecNorm(u, NORM_2, &norm_err));
  PetscCall(VecAXPY(u, 1.0, uexact));    /* restore u = solution for HDF5 */
  PetscCall(VecNorm(uexact, NORM_2, &norm_exact));
  relerr = (norm_exact > 1e-30) ? (norm_err / norm_exact) : norm_err;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Relative error (2-norm): %.6e\n", (double)relerr));

  /* [3(b)] Step 5: HDF5 output (PETSc --download-hdf5). Named "u", "f", "uexact" for plot_bvp.py. */
#if defined(PETSC_HAVE_HDF5)
  PetscViewer    viewer;
  PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, "bvp_solution.h5", FILE_MODE_WRITE, &viewer));
  PetscCall(PetscObjectSetName((PetscObject)uexact, "uexact"));
  PetscCall(PetscObjectSetName((PetscObject)f, "f"));
  PetscCall(PetscObjectSetName((PetscObject)u, "u"));
  PetscCall(VecView(f, viewer));
  PetscCall(VecView(u, viewer));
  PetscCall(VecView(uexact, viewer));
  PetscCall(PetscViewerDestroy(&viewer));
#else
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "HDF5 not available; skipping bvp_solution.h5\n"));
#endif

  PetscCall(KSPDestroy(&ksp));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&f));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&uexact));
  PetscCall(PetscFinalize());
  return 0;
}
