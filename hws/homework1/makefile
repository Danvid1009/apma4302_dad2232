## APMA 4302 HW1 - build expx with PETSc
##
## Usage:
##   make expx
##   mpiexec -n 4 ./expx -x -10 -N 50


##terminal commands: 
# cd "/Users/dan/Desktop/Columbia/HPC_4302/apma4302_dad2232/hws/homework1"
# make expx
#ls -l
# mpiexec -n 4 ./expx -x -10 -N 50

#to access the config.mk file (you should just use environmental variables)
include config.mk

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

expx: expx.o
	-${CLINKER} -o expx expx.o ${PETSC_LIB}
	${RM} expx.o

.PHONY: clean distclean

distclean: clean

clean::
	@rm -f *~ expx expx.o *tmp
