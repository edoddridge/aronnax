# External dependencies
HYPRE_DIR = lib/hypre/src

LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -lm

aronnax_core: aronnax.f90 Makefile
	gfortran -g -Ofast $< -o $@ -cpp

aronnax_test: aronnax.f90 Makefile
	gfortran -g -O1 -fcheck=all $< -o $@ -cpp

aronnax_prof: aronnax.f90 Makefile
	gfortran -g -pg -Ofast $< -o $@ -cpp

aronnax_external_solver: aronnax.f90 Makefile
  mpifort -g $< -o $@ -cpp -DuseExtSolver $(LIBS)
