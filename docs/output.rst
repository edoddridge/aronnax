Output
******

Aronnax produces output in two formats. The full fields are saved in Fortran unformatted files, that are compact and efficient, but not particularly user friendly. Layerwise statistics are saved into one csv file per variable.


Reading the data
===================

Aronnax includes a helper function for dealing with the unformatted output.

.. autofunction:: aronnax.interpret_raw_file