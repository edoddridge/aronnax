MIM: MIM.f90 Makefile
	gfortran -g -Ofast $< -o $@

MIM_external_solver: MIM.f90 Makefile
	mpif77 -g -Ofast $< -o $@ -include 'mpif.h'

MIM_test: MIM.f90 Makefile
	gfortran -g -O1 -fcheck=all $< -o $@

MIM_prof: MIM.f90 Makefile
	gfortran -g -pg -Ofast $< -o $@
