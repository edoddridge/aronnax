MIM: aronnax.f90 Makefile
	gfortran -g -Ofast $< -o $@

MIM_test: aronnax.f90 Makefile
	gfortran -g -O1 -fcheck=all $< -o $@

MIM_prof: aronnax.f90 Makefile
	gfortran -g -pg -Ofast $< -o $@
