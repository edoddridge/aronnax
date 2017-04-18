Input generators
*******************

Aronnax includes a number of helper functions for generating input fields.
These can be called with numerical arguments (generating very simple input fields)
directory from the aronnax.conf file, using the syntax `:generator_name:arg1,arg2,...,argn`.

Grid
=====
.. autoclass:: aronnax.Grid
   :members:
   :private-members:


Forcings
==========


.. autofunction:: aronnax.wind_x

.. autofunction:: aronnax.wind_y


Initial conditions
===================

.. autofunction:: aronnax.interpret_initial_heights


Domain fields
===============

.. autofunction:: aronnax.f_plane_u

.. autofunction:: aronnax.f_plane_v

.. autofunction:: aronnax.beta_plane_u

.. autofunction:: aronnax.beta_plane_v

.. autofunction:: aronnax.rectangular_pool
