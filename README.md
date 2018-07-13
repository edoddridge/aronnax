[![Build Status](https://travis-ci.org/edoddridge/aronnax.svg?branch=master)](https://travis-ci.org/edoddridge/aronnax)
[![codecov](https://codecov.io/gh/edoddridge/aronnax/branch/master/graph/badge.svg)](https://codecov.io/gh/edoddridge/aronnax)
[![Documentation Status](http://readthedocs.org/projects/aronnax/badge/?version=latest)](http://aronnax.readthedocs.io/en/latest/?badge=latest)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00592/status.svg)](https://doi.org/10.21105/joss.00592)



# Aronnax

An idealised isopycnal ocean circulation model that can either be run
as a reduced gravity model with n + 1/2 layers, or with n layers and
variable bathymetry.

Aronnax is
- [Easy to install](http://aronnax.readthedocs.io/en/latest/installation.html)
  on a laptop or a compute node, including without
  administrative privileges.
- Easy to configure.  All parameters, including grid size, are
  specified at runtime in a simple configuration file.
- [Easy to use](http://aronnax.readthedocs.io/en/latest/examples.html).
  Aronnax can be controlled programmatically as a Python library.
- Easy to learn and understand, with extensive [online
  documentation](http://aronnax.readthedocs.io/en/latest/index.html), including a
  complete description of [the
  physics](http://aronnax.readthedocs.io/en/latest/aronnax_model.html#the-physics)
  and [the
  numerics](http://aronnax.readthedocs.io/en/latest/aronnax_model.html#numerical-algorithm).
- [Verified](http://aronnax.readthedocs.io/en/latest/verification.html).
  Aronnax successfully [reproduces published results](http://aronnax.readthedocs.io/en/latest/examples.html#reproducing-published-results) from
  idealised models appearing in the literature.
- [Fast](http://aronnax.readthedocs.io/en/latest/benchmarks.html).  The
  main integration loop is a Fortran program, wrapped in
  Python for convenient use.


# Get Involved

- Please report any bugs you come across in the [Github issue
  tracker](https://github.com/edoddridge/aronnax/issues)

- Feature requests should also be made through the [Github issue
  tracker](https://github.com/edoddridge/aronnax/issues)

- We are happy to receive pull requests for new features or bug fixes;
  check out the [issues](https://github.com/edoddridge/aronnax/issues) for
  stuff that we know needs doing, and [HACKING.md](HACKING.md) for a
  developer's intro to the repository.
