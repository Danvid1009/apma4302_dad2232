## APMA 4302 HW1 - build expx with PETSc
##
## Usage:
##   make expx
##   mpiexec -n 4 ./expx -x -10 -N 50
##
## Or set env vars instead of config.mk:
##   export PETSC_DIR=/path/to/petsc  PETSC_ARCH=your-arch
##   make expx


##terminal commands: 
# cd "/Users/dan/Desktop/Columbia/HPC_4302/apma4302_dad2232/hws/homework1"
# make expx
#ls -l
# mpiexec -n 4 ./expx -x -10 -N 50

#to access the config.mk file (you should just use environmental variables)
include config.mk

# -----------------------------------------------------------------------------
# PETSc build system: variables = compiler, flags, libs; rules = how to build .o
# These files come from your PETSc install; they define CLINKER, PETSC_LIB,
# and rules to compile .c -> .o with the right flags and MPI.
# -----------------------------------------------------------------------------
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Extra C flags: strict C99, pedantic warnings (optional; matches p4pdes style).
CFLAGS += -pedantic -std=c99

# -----------------------------------------------------------------------------
# Target: expx (the executable)
# Prerequisite: expx.o (built automatically by PETSc rules from expx.c).
# Recipe: link expx.o with PETSc/MPI libs into ./expx, then remove the .o.
# The leading - on the CLINKER line means "ignore link errors" (PETSc quirk).
# -----------------------------------------------------------------------------
expx: expx.o
	-${CLINKER} -o expx expx.o ${PETSC_LIB}
	${RM} expx.o

# -----------------------------------------------------------------------------
# Phony targets: not real files; make always runs their recipes when asked.
# -----------------------------------------------------------------------------
.PHONY: clean distclean

# distclean = same as clean (no separate "dist" files here).
distclean: clean

# clean: remove build products (expx, expx.o) and editor backups (*~, *tmp).
# The @ suppresses echoing the command; :: allows multiple clean rules elsewhere.
# -----------------------------------------------------------------------------
clean::
	@rm -f *~ expx expx.o *tmp
