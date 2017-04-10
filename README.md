[![Build Status](https://travis-ci.org/edoddridge/aronnax.svg?branch=master)](https://travis-ci.org/edoddridge/aronnax)

# Aronnax

An idealized isopycnal ocean circulation model that can either be run
as a reduced gravity model with n + 1/2 layers, or with n layers and
variable bathymetry.

Aronnax is
- Easy to install on a laptop or a compute node, including without
  administrative privileges.
- Easy to configure.  All parameters, including grid size, are
  specified at runtime in a simple configuration file.
- [Easy to use](https://edoddridge.github.io/aronnax/examples.html).
  Aronnax can be called as a simple command-line program
  that reads and writes the standard NetCDF data format, or can be
  controlled programmatically as a Python library, communicating data
  through Numpy arrays.
- Easy to learn and understand, with extensive [online
  documentation](https://edoddridge.github.io/aronnax/), including a
  complete description of [the
  physics](https://edoddridge.github.io/aronnax/about_aronnax.html#the-physics)
  and [the
  numerics](https://edoddridge.github.io/aronnax/about_aronnax.html#discretisation).
- [Verified](https://edoddridge.github.io/aronnax/verification.html).
  Aronnax successfully reproduces multiple published results from
  idealized models appearing in the literature.
- [Fast](https://edoddridge.github.io/aronnax/benchmarks.html).  The
  main integration loop is a multi-core Fortran program, wrapped in
  Python for convenient use.

# Installation

todo

# Developing

The `test` directory contains a suite of automated tests written in
Python and runnable in
[pytest](http://doc.pytest.org/en/latest/contents.html).  The tests
depend on numpy and scipy.

As of this writing, that test suite checks that current outputs are
consistent with previously accepted results.  Thus, it is only useful
for checking that refactorings that are not expected to affect the
numerics indeed do not affect the numerics.

# Code style

Summary of Fortran code style enforced:

* Indentation of two spaces for all block constructs.
* Four-space indentation of continuation lines (with some exceptions
  to clarify structure of long algebraic expressions).
* Lines should have no trailing whitespace.
* The last line in the file should have one newline, and no additional blank lines.
* Lines should be < 80 characters wide, but I was not successful in
  breaking all over-long lines. The remainder are predominantly line
  comments that I did not wish to move off their line, not knowing
  what exactly they referred to.
* Spell language keywords in lower case.
* All subroutines should declare implicit none
* All subroutines should declare all their arguments, in order of
  receipt, with intents.  Any local variables are declared in a
  separate block.  (Note: "no intent" is actually different from any
  of the declarable intents, including `intent (inout)`, but so far,
  Aronnax contains no subroutines that take any arguments that exhibit
  that i/o behavior.)
* Spell block terminators end do, end if, and end subroutine `<name>`, not enddo or endif.
* Commas should be followed by whitespace, as in written natural languages. 
  * (Exception: array subscripts in dense algebraic expressions, as
    they should be compact enough that the eye groups the array
    element as one thing.)
* The equal sign should be offset with spaces when it means
  assignment, and not offset when it means value for keyword argument.
