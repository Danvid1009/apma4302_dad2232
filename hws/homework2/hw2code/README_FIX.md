# Why you're seeing both errors

You have **MPICH** linked (we switched earlier), but **PETSc is still built with Open MPI**. So:

1. **"PETSc was configured with Open MPI but now ... non-Open MPI"**  
   → You're compiling with MPICH's `mpicc` against a PETSc that was built for Open MPI. They must match.

2. **MPI_ERR_COMM when running**  
   → You're either running an old `bvp` binary (from when Open MPI was linked), or the run is using Open MPI's `mpiexec`; either way it hits the Open MPI bug on Mac.

**Fix:** Rebuild PETSc with MPICH (configure + make all), then build and run bvp with MPICH. Configure did **not** finish last time (log stops early); it has to run to completion once.

---

## Option 1: Run configure and wait (~20 min)

In a terminal, run this and **leave it running** until it prints a summary and returns to the prompt (do not Ctrl+C):

```bash
brew link mpich --force
cd /Users/dan/Desktop/Columbia/HPC_4302/petsc
export PETSC_DIR=$PWD PETSC_ARCH=apma4302-base-opt
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
  --with-debugging=0 --with-shared-libraries=1 \
  --COPTFLAGS="-O3 -march=native" --CXXOPTFLAGS="-O3 -march=native" \
  --FOPTFLAGS="-O3 -march=native"
```

When it finishes, run:

```bash
make all
```

Then in this folder:

```bash
cd "/Users/dan/Desktop/Columbia/HPC_4302/apma4302_dad2232/hws/homework2/hw2code"
export PETSC_DIR=/Users/dan/Desktop/Columbia/HPC_4302/petsc PETSC_ARCH=apma4302-base-opt
make bvp
mpiexec -np 1 ./bvp -options_file options_file
```

---

## Option 2: Use Open MPI (build works, run still fails on Mac)

If you only need to **compile** (e.g. for submission) and will run elsewhere:

```bash
brew unlink mpich
brew link open-mpi --force
cd "/Users/dan/Desktop/Columbia/HPC_4302/apma4302_dad2232/hws/homework2/hw2code"
export PETSC_DIR=/Users/dan/Desktop/Columbia/HPC_4302/petsc PETSC_ARCH=apma4302-base-opt
make bvp
```

`make bvp` will succeed. `mpiexec ./bvp ...` will still hit MPI_ERR_COMM on your Mac; run on a cluster or with PETSc built with MPICH to get a working run.
