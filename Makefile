# External dependencies
HYPRE_DIR = lib/hypre/src

LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -lm

aronnax_core: aronnax.f90 Makefile
	mpifort -g -Ofast $< -o $@ -cpp

aronnax_test: aronnax.f90 Makefile
	mpifort -g -O1 -fcheck=all $< -o $@ -cpp

aronnax_prof: aronnax.f90 Makefile
	mpifort -g -pg -Ofast $< -o $@ -cpp

aronnax_external_solver_test: aronnax.f90 Makefile
	mpifort -g $< -o $@ -cpp -DuseExtSolver $(LIBS)

aronnax_external_solver: aronnax.f90 Makefile
	mpifort -g -Ofast $< -o $@ -cpp -DuseExtSolver $(LIBS)