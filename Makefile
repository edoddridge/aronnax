# External dependencies
HYPRE_DIR = lib/hypre/src

LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -lm

aronnax_core: aronnax.f90 Makefile
	mpif90 -g -Ofast -fno-stack-arrays $< -o $@ -cpp

aronnax_test: aronnax.f90 Makefile
	mpif90 -g -fprofile-arcs -ftest-coverage -O1 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -Wuninitialized $< -o $@ -cpp

aronnax_prof: aronnax.f90 Makefile
	mpif90 -g -pg -Ofast $< -o $@ -cpp

aronnax_external_solver_test: aronnax.f90 Makefile
	mpif90 -g -fprofile-arcs -ftest-coverage -O1 -ffpe-trap=invalid,zero,overflow,underflow -fcheck=all -Wuninitialized $< -o $@ -cpp -DuseExtSolver $(LIBS)

aronnax_external_solver: aronnax.f90 Makefile
	mpif90 -g -Ofast -fno-stack-arrays $< -o $@ -cpp -DuseExtSolver $(LIBS)
