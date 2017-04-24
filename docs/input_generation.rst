Input generators
*******************



Grid
=====
.. autoclass:: aronnax.Grid
   :members:
   :private-members:



Inputs
==========
Aronnax includes a number of helper functions for generating input fields.
These can be called with numerical arguments (generating very simple input fields)
directly from the aronnax.conf file, using the syntax `:generator_name:arg1,arg2,...,argn`.

.. autofunction:: aronnax.tracer_point_variable_3d

.. autofunction:: aronnax.u_point_variable_3d

.. autofunction:: aronnax.v_point_variable_3d




Coriolis fields
++++++++++++++++

.. autofunction:: aronnax.f_plane_f_u

.. autofunction:: aronnax.f_plane_f_v

.. autofunction:: aronnax.beta_plane_f_u

.. autofunction:: aronnax.beta_plane_f_v


Domain shape
+++++++++++++

.. autofunction:: aronnax.rectangular_pool


Custom generator functions
+++++++++++++++++++++++++++
The use of custom input generator functions, described in more detail in <running_aronnax.html>, allows the model to be passed user defined functions of `X` and `Y`, created from a `numpy.meshgrid` call on the appropriate axes, evaluate them to create the input fields and save these fields to the 'input' folder in the format required by the Fortran core.