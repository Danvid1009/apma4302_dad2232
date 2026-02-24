# Fixing MPI_ERR_COMM / Open MPI 5 on macOS

Your `bvp` fails with **MPI error 5 MPI_ERR_COMM: invalid communicator** inside PETSc’s `MatCreate` (which calls `MPI_Comm_dup`). Below are the most likely fixes.

---

## 1. macOS Firewall (try this first)

Open MPI 5 uses **prte** (PRRTE). If the macOS firewall is blocking it, you get odd MPI errors (including invalid communicator).

**Check and fix:**
1. **System Settings → Network → Firewall → Options**
2. If **prte** (or **prterun**) is listed and set to “Block incoming connections”, change it to **Allow**.
3. If you’ve never seen a firewall popup for MPI, run once and when the “prte wants to accept incoming connections” dialog appears, click **Allow** (you may need an admin password).
4. Optionally add **mpiexec** and your **bvp** binary to the firewall allow list.

**Quick test:** Temporarily turn the firewall off and run:
```bash
mpiexec -np 1 ./bvp -options_file options_file
```
If it works with the firewall off, the fix is to allow prte (and related MPI tools) with the firewall on.

---

## 2. Use Open MPI 4.1 instead of 5

Open MPI 5 on macOS can require a working network/loopback setup; Open MPI 4.1 has a fallback that often works better on Macs.

**Option A – Homebrew (if a 4.x formula exists):**
```bash
# Check if open-mpi@4 exists
brew search open-mpi

# If you see open-mpi@4:
brew install open-mpi@4
brew unlink open-mpi
brew link open-mpi@4 --force
```
Then **rebuild PETSc** with this MPI (re-run your configure script so PETSc picks up the new `mpicc`/`mpiexec`), and rebuild `bvp`.

**Option B – Build Open MPI 4.1 from source:**
1. Download Open MPI 4.1.x from https://www.open-mpi.org/software/ompi/v4.1/
2. Build and install to a prefix (e.g. `$HOME/opt/openmpi-4.1`).
3. Set `PATH` and `LD_LIBRARY_PATH` (or equivalent) so that `mpicc` and `mpiexec` come from this install.
4. Reconfigure and rebuild PETSc with this MPI, then rebuild `bvp`.

---

## 3. Try a different PETSc config

You have three PETSc builds; all use the same Open MPI 5, but one might behave differently:

```bash
export PETSC_DIR=/Users/dan/Desktop/Columbia/HPC_4302/petsc

# Try debug build
export PETSC_ARCH=apma4302-base-debug
make clean-bvp && make bvp
mpiexec -np 1 ./bvp -options_file options_file

# Or pkgs-opt
export PETSC_ARCH=apma4302-pkgs-opt
make clean-bvp && make bvp
mpiexec -np 1 ./bvp -options_file options_file
```

---

## Summary

1. **First:** Allow **prte** (and MPI apps) in the macOS firewall, or test with the firewall off.
2. **If that doesn’t help:** Switch to Open MPI 4.1 (Homebrew or from source), rebuild PETSc and `bvp` against it.
3. **Optional:** Try `apma4302-base-debug` or `apma4302-pkgs-opt` in case one of them avoids the bad path.

Reference: [Open MPI 5 doesn't work after update to macOS 15.3.1](https://github.com/open-mpi/ompi/issues/13129) (resolved by allowing prte in the firewall).
