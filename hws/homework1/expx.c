/*
  APMA 4302 - Homework 1, Problem 1
  Parallel N-term Taylor approximation of exp(x). Based on p4pdes/c/ch1/e.c

  Run: mpiexec -n nP ./expx -x <x> -N <N>
  Uses block partition of [0,N) across ranks; for x<0 uses exp(x)=1/exp(-x).

  Notes:
  - Computes the Taylor sum in parallel with O(N/nP) work per process.
  - Uses a "balanced" strategy for negative x: exp(x) = 1/exp(-x).
*/

#include <petsc.h>
#include <math.h>

/*
  term_at_index(x, k): compute the k-th Taylor term  x^k / k!
  Uses log-gamma so we don't overflow:  x^k/k! = exp(k*log(x) - lgamma(k+1)).
  Static = visible only in this file.
*/
static PetscReal term_at_index(PetscReal x, PetscInt k)
{
  /* x^0/0! = 1; for x=0 all other terms are 0 (avoid log(0)). */
  if (x == (PetscReal)0) return (k == 0) ? (PetscReal)1 : (PetscReal)0;
  /* x^k/k! = exp(k*log(x) - ln(k!)) and ln(k!) = lgamma(k+1). */
  return (PetscReal)exp((double)k * log((double)x) - lgamma((double)k + 1.0));
}

int main(int argc, char **argv) {
  /* MPI: rank (this process id), size (total number of processes). */
  PetscMPIInt  rank, size;
  /* N = number of Taylor terms; kStart,kEnd = this rank's term range [kStart, kEnd). */
  PetscInt     N = 20, k, kStart, kEnd;
  /* x = evaluation point; y = |x| for the series; localSum/globalSum = partial/full sum. */
  PetscReal    x = 1.0, y, localSum, globalSum, term;
  /* approx = our sum (or 1/sum for x<0); exact = exp(x); relerr = relative error. */
  PetscReal    approx, exact, relerr;

  /* Initialize PETSc and MPI; register help string for -help. */
  PetscCall(PetscInitialize(&argc, &argv, NULL,
      "Approximate exp(x) with N Taylor terms in parallel.\n\n"));
  PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCall(MPI_Comm_size(PETSC_COMM_WORLD, &size));

  /* Parse command-line options: -x <value>, -N <value>. */
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "expx", NULL);
  PetscCall(PetscOptionsReal("-x", "Point at which to evaluate exp(x)", NULL, x, &x, NULL));
  PetscCall(PetscOptionsInt("-N", "Number of Taylor terms", NULL, N, &N, NULL));
  PetscOptionsEnd();
  PetscCheck(N > 0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "-N must be > 0 (got %" PetscInt_FMT ")", N);

  /* For x<0 we sum exp(|x|) then use exp(x)=1/exp(-x) for better numerics. */
  y = (x < 0) ? -x : x;
  /* Block partition: rank r gets indices [kStart, kEnd) with kStart = floor(r*N/size), kEnd = floor((r+1)*N/size). */
  kStart = (PetscInt)(((PetscInt64)rank * (PetscInt64)N) / (PetscInt64)size);
  kEnd   = (PetscInt)(((PetscInt64)(rank + 1) * (PetscInt64)N) / (PetscInt64)size);

  /* Compute this rank's chunk of the Taylor sum: sum_{k=kStart}^{kEnd-1} y^k/k! */
  localSum = 0.0;
  term = term_at_index(y, kStart);   /* first term in this chunk */
  for (k = kStart; k < kEnd; k++) {
    localSum += term;
    term *= y / (PetscReal)(k + 1);  /* recurrence: term_{k+1} = term_k * y/(k+1) */
  }

  /* Sum all ranks' localSum into globalSum (everyone gets the same value). */
  PetscCall(MPI_Allreduce(&localSum, &globalSum, 1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD));

  /* globalSum = sum of y^k/k! over k=0..N-1. For x>=0 that's approx exp(x); for x<0, exp(x)=1/exp(-x)=1/globalSum. */
  approx = (x < 0) ? ((PetscReal)1 / globalSum) : globalSum;
  exact  = (PetscReal)exp((double)x);
  relerr = (exact != (PetscReal)0)
      ? (PetscReal)fabs((double)(approx - exact) / (double)exact)
      : (PetscReal)fabs((double)(approx - exact));

  /* Per-rank report (like e.c): each rank prints how many terms it computed; order may vary. */
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "rank %d did %d terms\n", rank, (int)(kEnd - kStart)));
  /* Collective: one line for the result (printed by rank 0). */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
      "exp(x) is about %.17g  (x=%g, N=%" PetscInt_FMT ", nP=%d, relerr=%.6e)\n",
      (double)approx, (double)x, N, (int)size, (double)relerr));

  PetscCall(PetscFinalize());
  return 0;
}
