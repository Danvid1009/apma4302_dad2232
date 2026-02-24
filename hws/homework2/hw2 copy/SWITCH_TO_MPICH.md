# Switch to MPICH (fix MPI on macOS)

Homebrew's Open MPI 5 can misbehave on macOS (firewall/prte, MPI_ERR_COMM). **MPICH** is another MPI implementation that works with PETSc and avoids those issues.

**Note:** Homebrew's `mpich` and `open-mpi` conflict (both provide `mpicc`/`mpiexec`). You'll use MPICH for PETSc and this course; you can switch back to Open MPI later with `brew unlink mpich && brew link open-mpi`.

---

## Step 1: Install MPICH and make it active

```bash
# Install MPICH (if not already installed)
brew install mpich

# Use MPICH instead of Open MPI (so mpicc/mpiexec come from MPICH)
brew unlink open-mpi
brew link mpich --force
```

Verify:
```bash
which mpicc mpiexec
# Should show /opt/homebrew/bin/mpicc and /opt/homebrew/bin/mpiexec (from MPICH)
mpiexec --version
# Should say "MPICH" (or "Intel MPI" etc.), not "Open MPI"
```

---

## Step 2: Reconfigure PETSc with MPICH

PETSc was built with Open MPI; it must be rebuilt with MPICH's compilers.

```bash
cd /Users/dan/Desktop/Columbia/HPC_4302/petsc

# Use one of your configs (e.g. base-opt). Configure with current mpicc (now MPICH).
export PETSC_DIR=/Users/dan/Desktop/Columbia/HPC_4302/petsc
export PETSC_ARCH=apma4302-base-debug

# Re-run configure so PETSc picks up MPICH
./configure \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpif90 \
  --with-debugging=0 \
  --with-shared-libraries=1 \
  --COPTFLAGS="-O3 -march=native" \
  --CXXOPTFLAGS="-O3 -march=native" \
  --FOPTFLAGS="-O3 -march=native"

# Rebuild PETSc
make all
```

If `mpif90` is not provided by MPICH or fails, you can try configuring without Fortran (PETSc may allow that for a C-only build; if not, keep FC).

---

## Step 3: Rebuild bvp

```bash
cd "/Users/dan/Desktop/Columbia/HPC_4302/apma4302_dad2232/hws/homework2/hw2 copy"
export PETSC_DIR=/Users/dan/Desktop/Columbia/HPC_4302/petsc
export PETSC_ARCH=apma4302-base-debug
make clean-bvp
make bvp
```

---

## Step 4: Run bvp

```bash
mpiexec -np 1 ./bvp -options_file options_file
```

Use `mpirun` if your MPICH install prefers it: `mpirun -np 1 ./bvp -options_file options_file`.

---

## Switching back to Open MPI later

```bash
brew unlink mpich
brew link open-mpi --force
```

Then you’d need to reconfigure and rebuild PETSc again if you want to use Open MPI with PETSc.

---

## If PETSc configure fails on Fortran (mpif90)

MPICH from Homebrew may provide `mpifort` instead of `mpif90`, or no Fortran. Try:

```bash
# See what's available
ls /opt/homebrew/bin/mpi*

# If you have mpifort but not mpif90, create a symlink or use:
./configure ... --with-fc=mpifort ...

# If you must disable Fortran (PETSc may still build for C-only use):
./configure ... --with-fc=0
```

---

## Alternative: Open MPI 4.1 from source

If you prefer to stay with Open MPI but use 4.1.x (avoids the macOS 5.x issues):

1. Download: https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.8.tar.gz  
2. Build and install to a prefix (e.g. `$HOME/opt/openmpi-4.1`):
   ```bash
   tar xzf openmpi-4.1.8.tar.gz && cd openmpi-4.1.8
   ./configure --prefix=$HOME/opt/openmpi-4.1 --with-pmix=internal
   make -j4 && make install
   ```
3. Use that MPI for PETSc:
   ```bash
   export PATH=$HOME/opt/openmpi-4.1/bin:$PATH
   export LD_LIBRARY_PATH=$HOME/opt/openmpi-4.1/lib:$LD_LIBRARY_PATH
   # Then reconfigure and build PETSc (Step 2 above)
   ```
