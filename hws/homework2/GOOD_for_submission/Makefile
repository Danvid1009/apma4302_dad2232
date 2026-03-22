# Build BVP solver. Requires PETSc (with HDF5 for Q3 HDF5 output).
# Set PETSC_DIR; set PETSC_ARCH only if you use a named arch.
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

bvp: bvp.o
	$(CLINKER) -o bvp bvp.o $(PETSC_LIB)
	$(RM) bvp.o

petsc_minimal: petsc_minimal.o
	$(CLINKER) -o petsc_minimal petsc_minimal.o $(PETSC_LIB)
	$(RM) petsc_minimal.o

clean-bvp:
	$(RM) bvp bvp.o petsc_minimal petsc_minimal.o
