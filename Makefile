MIM: MIM.f90 Makefile
	gfortran -g -Ofast $< -o $@

MIM_test: MIM.f90 Makefile
	gfortran -g -O1 -fcheck=all $< -o $@

MIM_prof: MIM.f90 Makefile
	gfortran -g -pg -Ofast $< -o $@
