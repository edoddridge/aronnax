# External dependencies
HYPRE_DIR = lib/hypre/src

LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -lm

src_dir = src/

TEST_OPTS = -g -fprofile-arcs -ftest-coverage -O1 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -Wuninitialized -Werror

CORE_OPTS = -g -Ofast -fno-stack-arrays

PROF_OPTS = -g -pg -Ofast

FILES = declarations advection_schemes adams_bashforth end_run enforce_thickness boundaries vorticity momentum io thickness bernoulli state_deriv time_stepping barotropic_mode model_main aronnax

TEST_objects = $(patsubst %, $(src_dir)%_TEST.o, $(FILES))
CORE_objects = $(patsubst %, $(src_dir)%_CORE.o, $(FILES))
PROF_objects = $(patsubst %, $(src_dir)%_PROF.o, $(FILES))

HYPRE_TEST_objects = $(patsubst %, $(src_dir)%_HYPRE_TEST.o, $(FILES))
HYPRE_CORE_objects = $(patsubst %, $(src_dir)%_HYPRE_CORE.o, $(FILES))
HYPRE_PROF_objects = $(patsubst %, $(src_dir)%_HYPRE_PROF.o, $(FILES))

# Profiling execuable
aronnax_prof: $(PROF_objects) Makefile
	mpif90 $(PROF_OPTS) $< -o $@ -cpp

%_PROF.o: %.f90
	mpif90 $(PROF_OPTS) -J $(src_dir) -c $< -o $@ -cpp

# Testing executables with compile-time flags for catching errors
aronnax_test: $(TEST_objects) Makefile
	mpif90 $(TEST_objects) $(TEST_OPTS) -o $@ -cpp

%_TEST.o: %.f90
	mpif90 $(TEST_OPTS) -J $(src_dir) -c $< -o $@ -cpp

aronnax_external_solver_test: $(HYPRE_TEST_objects) Makefile
	mpif90 $(HYPRE_TEST_objects) $(TEST_OPTS) -o $@ -cpp -DuseExtSolver $(LIBS)

%_HYPRE_TEST.o: %.f90
	mpif90 $(TEST_OPTS) -J $(src_dir) -c $< -o $@ -cpp -DuseExtSolver $(LIBS)

# Highly optimised executables for faster simulations
aronnax_core: $(CORE_objects) Makefile
	mpif90 $(CORE_objects) $(CORE_OPTS) -o $@ -cpp

%_CORE.o: %.f90
	mpif90 $(CORE_OPTS) -J $(src_dir) -c $< -o $@ -cpp

aronnax_external_solver: $(HYPRE_CORE_objects) Makefile
	mpif90 $(HYPRE_CORE_objects) $(CORE_OPTS) -o $@ -cpp -DuseExtSolver $(LIBS)

%_HYPRE_CORE.o: %.f90
	mpif90 $(CORE_OPTS) -J $(src_dir) -c $< -o $@ -cpp -DuseExtSolver $(LIBS)

# shortcuts for removing compiled files
clean:
	rm $(src_dir)*.o $(src_dir)*.gcno $(src_dir)*.gcda $(src_dir)*.mod
