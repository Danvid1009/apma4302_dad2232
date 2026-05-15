#!/usr/bin/env python3
"""Quick PETSc environment checker for apma4302-pkgs-opt."""

from __future__ import annotations

import os
import sys


def main() -> int:
    expected_arch = "apma4302-pkgs-opt"
    expected_dir = "/Users/dan/Desktop/Columbia/HPC_4302/petsc"
    petsc_arch = os.environ.get("PETSC_ARCH", "")
    petsc_dir = os.environ.get("PETSC_DIR", "")
    pythonpath = os.environ.get("PYTHONPATH", "")
    path = os.environ.get("PATH", "")

    print(f"PETSC_DIR={petsc_dir or '<unset>'}")
    print(f"PETSC_ARCH={petsc_arch or '<unset>'}")
    print(f"python_executable={sys.executable}")

    missing = []
    if petsc_dir != expected_dir:
        missing.append(f"PETSC_DIR should be '{expected_dir}'")
    if petsc_arch != expected_arch:
        missing.append(f"PETSC_ARCH should be '{expected_arch}'")

    if "apma4302-pkgs-opt/lib" not in pythonpath:
        missing.append("PYTHONPATH missing PETSc lib path for apma4302-pkgs-opt")

    if "apma4302-pkgs-opt/bin" not in path:
        missing.append("PATH missing PETSc bin path for apma4302-pkgs-opt")

    try:
        from petsc4py import PETSc  # pylint: disable=import-outside-toplevel
    except Exception as exc:  # pragma: no cover
        print(f"petsc4py_import=FAIL ({exc})")
        missing.append("petsc4py failed to import")
        PETSc = None  # type: ignore
    else:
        print("petsc4py_import=OK")
        print(f"PETSc version={PETSc.Sys.getVersion()}")
        print(f"MPI size={PETSc.COMM_WORLD.getSize()}")

    if missing:
        print("\nStatus: NOT READY")
        for item in missing:
            print(f"- {item}")
        return 1

    print("\nStatus: READY")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
