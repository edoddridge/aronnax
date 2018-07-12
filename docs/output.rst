Output
******

Aronnax produces output in two formats. The full fields are saved in Fortran unformatted files, that are compact and efficient, but not particularly user friendly. Layerwise statistics are saved into one csv file per variable.


Reading the data
===================

Aronnax includes a number of helper functions for dealing with the unformatted Fortran output.

The simplest of these is `interpret_raw_file` which loads a single output of a single variable into a numpy array.

.. autofunction:: aronnax.interpret_raw_file

Aronnax also ships with a function for lazily loading multiple timestamps of a single variable. This function uses `xarray <http://xarray.pydata.org/en/stable/>`_ and `dask <http://dask.pydata.org/en/latest/docs.html>`_ to create labeled n-dimensional arrays.

.. autofunction:: aronnax.open_mfdataarray
