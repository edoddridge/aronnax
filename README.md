[![Build Status](https://travis-ci.org/edoddridge/MIM.svg?branch=master)](https://travis-ci.org/edoddridge/MIM)

# MIM
Minimalist Isopycnal Model: a simplistic isopycnal model that can either be run as a reduced gravity model with n + 1/2 layers, or with n layers and variable bathymetry.

The model is written in fortran 90 and exports the data as unformatted fortran data files. Most parameters are specified at runtime in the 'parameters.in' file. The exceptions are the number of grid points in the x and y directions, and the number of layers which must be specified when the code is compiled.

The online documentation can be accessed [here](https://edoddridge.github.io/MIM/), or it can be found in the docs folder.

# Developing

The `test` directory contains a suite of automated tests written in
Python and runnable in
[pytest](http://doc.pytest.org/en/latest/contents.html).  The tests
depend on numpy, scipy, and the `MIMutils` package from the `utils`
directory.

As of this writing, that test suite checks bit-for-bit consistency of
results across changes to MIM.  Thus, it is only useful for checking
that refactorings that are not expected to affect the numerics indeed
do not affect the numerics.
