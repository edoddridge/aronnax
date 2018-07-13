.. Aronnax documentation master file, created by
   sphinx-quickstart on Wed Apr  5 22:01:04 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Aronnax
=======

Aronnax is an idealised and easily configurable isopycnal ocean circulation model.
Aronnax can be run as either a reduced gravity model with n + 1/2 layers, or with n layers and variable bathymetry.

Aronnax is

- `Easy to install <installation.html>`_
  on a laptop or a compute node, including without
  administrative privileges.

- Easy to configure.  All parameters, including grid size, are
  specified at runtime in a simple configuration file or function call.

- `Easy to use <examples.html>`_.
  Aronnax can be run from the shell as a Python script.

- Easy to learn and understand, with extensive online
  documentation, including a complete description of `the
  physics <aronnax_model.html#the-physics>`_
  and `the numerics <aronnax_model.html#discretisation>`_.

- `Verified <verification.html>`_.
  Aronnax successfully reproduces published results from
  idealised models appearing in the literature.

- `Fast <benchmarks.html>`_.  The
  main integration loop is a Fortran program, wrapped in
  Python for convenient use.

.. toctree::
   :numbered: 3
   :maxdepth: 2
   :caption: Contents:

   aronnax_model
   installation
   input_generation
   running_aronnax
   output
   examples
   verification
   benchmarks
   publications
   development


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
