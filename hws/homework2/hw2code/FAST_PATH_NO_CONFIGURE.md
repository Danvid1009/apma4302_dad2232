# Get BVP running without waiting for PETSc configure

## Option A: Use the course cluster (fastest if you have access)

If APMA 4302 gives you access to a Columbia HPC cluster (e.g. via ssh):

1. Copy your hw2 folder to the cluster (or clone your repo there).
2. Log in and load the environment (example; adjust to your cluster’s modules):
   ```bash
   module load petsc   # or: module load openmpi petsc
   ```
3. Go to the directory with `bvp.c` and the Makefile, then:
   ```bash
   export PETSC_DIR=$PETSC_DIR   # often set by module
   export PETSC_ARCH=$PETSC_ARCH # often set by module; may be empty
   make bvp
   mpiexec -np 1 ./bvp -options_file options_file
   ```
No local configure; build + run usually under 2 minutes.

---

## Option B: Conda PETSc (~2–5 min install, no configure)

Use a pre-built PETSc from conda-forge so you don’t run `./configure` at all.

1. **Install PETSc and MPI** (in a terminal, with network):
   ```bash
   conda install -c conda-forge petsc mpich
   ```
   This can take 2–5 minutes.

2. **Find PETSc location** (conda puts it in your env):
   ```bash
   echo $CONDA_PREFIX
   ls $CONDA_PREFIX/lib/petsc/conf/variables 2>/dev/null && echo "Found"
   ```
   If `variables` is not there, try:
   ```bash
   ls $CONDA_PREFIX/lib/petsc/conf/
   ```
   If you see `variables` under a subdir (e.g. under an arch name), set `PETSC_ARCH` to that folder name.

3. **Build bvp** from the hw2 folder:
   ```bash
   cd "/Users/dan/Desktop/Columbia/HPC_4302/apma4302_dad2232/hws/homework2/hw2 copy"
   export PETSC_DIR=$CONDA_PREFIX
   export PETSC_ARCH=   # leave empty if conf is directly under $CONDA_PREFIX/lib/petsc/conf
   make bvp
   ```
   If `make` says it can’t find `variables`, set `PETSC_ARCH` to the single directory name you see under `$CONDA_PREFIX/lib/petsc/conf/` (e.g. `export PETSC_ARCH=arch-macosx-arm64` or similar).

4. **Run:**
   ```bash
   mpiexec -np 1 ./bvp -options_file options_file
   ```
   Conda’s `mpiexec` should match the PETSc from conda.

---

## If Option B fails

- Confirm you’re in the conda env that has `petsc`: `conda list petsc`.
- If the Makefile can’t find `variables`, the conda package layout may differ; then the only local option without waiting is **Option A (cluster)**.

Summary: **Cluster = no configure, ~2 min. Conda = one install, no configure, then build/run.**
