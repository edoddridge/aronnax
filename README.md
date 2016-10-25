# MIM
Minimalist Isopycnal Model: a simplistic isopycnal model with n layers and variable bathymetry.

The model is written in fortran 90 and exports the data as unformatted fortran data files. Most parameters are specified at runtime in the 'parameters.in' file. The exceptions are the number of grid points in the x and y directions, and the number of layers which are specified when the code is compiled.

This model is largely based on the n + 1/2 layer model described [here](https://doddridge.me/code/).