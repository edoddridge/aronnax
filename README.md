[![Build Status](https://travis-ci.org/edoddridge/aronnax.svg?branch=master)](https://travis-ci.org/edoddridge/aronnax)

# Aronnax

An idealized isopycnal ocean circulation model that can either be run
as a reduced gravity model with n + 1/2 layers, or with n layers and
variable bathymetry.

Aronnax is
- [Easy to install](https://github.com/edoddridge/aronnax#install)
  on a laptop or a compute node, including without
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

# Install

todo

# Get Involved

- TODO Contact email, and/or list?
- Please report any bugs you come across in the [Github issue
  tracker](https://github.com/edoddridge/aronnax/issues)

- We are happy to receive pull requests for new features or bug fixes;
  check out the [issues](https://github.com/edoddridge/aronnax/issues) for
  stuff that we know needs doing, and [HACKING.md](HACKING.md) for a
  developer's intro to the repository.
