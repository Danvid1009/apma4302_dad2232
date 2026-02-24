#include <petsc.h>
#include <petscviewerhdf5.h>
int main(int argc, char **argv) {
  Mat A;
  PetscInt m = 40, k = 1;
  PetscReal gamma = 0.0, c = 1.0;
  PetscCall(PetscInitialize(&argc, &argv, NULL, "BVP -u'' + gamma*u = f on [0,1], Dirichlet BCs.\n\n"));
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "bvp", NULL);
  PetscCall(PetscOptionsInt("-bvp_m", "m", NULL, m, &m, NULL));
  PetscCall(PetscOptionsReal("-bvp_gamma", "gamma", NULL, gamma, &gamma, NULL));
  PetscCall(PetscOptionsInt("-bvp_k", "k", NULL, k, &k, NULL));
  PetscCall(PetscOptionsReal("-bvp_c", "c", NULL, c, &c, NULL));
  PetscOptionsEnd();
  PetscCheck(m > 0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "bvp_m must be > 0");
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatDestroy(&A));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "PETSc Options+Mat OK\n"));
  PetscCall(PetscFinalize());
  return 0;
}
