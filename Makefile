aronnax_core: aronnax.f90 Makefile
	gfortran -g -Ofast $< -o $@

aronnax_test: aronnax.f90 Makefile
	gfortran -g -O1 -fcheck=all $< -o $@

aronnax_prof: aronnax.f90 Makefile
	gfortran -g -pg -Ofast $< -o $@
