MIM: MIM.f90 Makefile
	gfortran -g -Ofast $< -o $@ -cpp

MIM_external_solver: MIM.f90 Makefile
	mpifort -g -Ofast $< -o $@ -cpp -DuseMPI

MIM_test: MIM.f90 Makefile
	gfortran -g -O1 -fcheck=all $< -o $@ -cpp

MIM_prof: MIM.f90 Makefile
	gfortran -g -pg -Ofast $< -o $@ -cpp
