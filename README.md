[![Build Status](https://travis-ci.org/edoddridge/MIM.svg?branch=master)](https://travis-ci.org/edoddridge/MIM)

# MIM
Minimalist Isopycnal Model: a simplistic isopycnal model that can either be run as a reduced gravity model with n + 1/2 layers, or with n layers and variable bathymetry.

The model is written in Fortran 90 and exports the data as unformatted fortran data files. All parameters, including grid size, are specified at runtime in the 'parameters.in' file.

The online documentation can be accessed [here](https://edoddridge.github.io/MIM/), or it can be found in the docs folder.

# Developing

The `test` directory contains a suite of automated tests written in
Python and runnable in
[pytest](http://doc.pytest.org/en/latest/contents.html).  The tests
depend on numpy, scipy, and the `MIMutils` package from the `utils`
directory.

As of this writing, that test suite checks that current outputs are consistent with previously accepted results.  Thus, it is only useful for checking
that refactorings that are not expected to affect the numerics indeed
do not affect the numerics.

# Code style

Summary of Fortran code style enforced:

* Indentation of two spaces for all block constructs.
* Four-space indentation of continuation lines (with some exceptions to clarify structure of long algebraic expressions).
* Lines should have no trailing whitespace.
* The last line in the file should have one newline, and no additional blank lines.
* Lines should be < 80 characters wide, but I was not successful in breaking all over-long lines. The remainder are predominantly line comments that I did not wish to move off their line, not knowing what exactly they referred to.
* Spell language keywords in lower case.
* All subroutines should declare implicit none
* Spell block terminators end do, end if, and end subroutine `<name>`, not enddo or endif.
* Commas should be followed by whitespace, as in written natural languages. 
 * (Exception: array subscripts in dense algebraic expressions, as they should be compact enough that the eye groups the array element as one thing.)
* The equal sign should be offset with spaces when it means assignment, and not offset when it means value for keyword argument.
