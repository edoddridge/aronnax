Input generation
*******************

There are a number of helpers designed to ease the creation of input fields.


Grid object
============
The grid object contains the axes for variables on tracer points, both velocity points, and vorticity points. This allows users to create their inputs in a resolution independent manner by specifying functions based on physical location, rather than grid point number. Changing resolution therefore becomes a trivial exercise that requires changing only the `nx`, `ny`, `dx`, and `dy` parameters.


.. autoclass:: aronnax.Grid
   :members:
   :private-members:



Inputs
==========

Custom generator functions
+++++++++++++++++++++++++++
The use of custom input generator functions, described in more detail below, allows the model to be passed user defined functions of `X` and `Y`, or layerwise constant values for the input fields.

Aronnax will evaluate the user defined functions or constants to create the input fields, and save these fields to the 'input' folder in the format required by the Fortran core.


The generator functions can be called either directly from the aronnax.conf file, using the syntax `:generator_name:arg1,arg2,...,argn`, or they can be included in the call to :meth:`aronnax.driver.simulate` using the syntax `field_name=[arg1, arg2, ..., argn]`. 

If they are used from the configuration file then the arguments must be numerical (generating very simple layerwise constant value input fields). If the generator functions are included in the call to :meth:`aronnax.driver.simulate`, then the arguments may be either numerical or functions, but must be passed as a list. That is, they must be wrapped in square brackets [], even if that list has only one element. Aronnax will check that the length of the list equals the number of layers the field is expected to have and will throw an error if they are not equal.


Generic generator functions 
----------------------------
.. autofunction:: aronnax.tracer_point_variable

.. autofunction:: aronnax.u_point_variable

.. autofunction:: aronnax.v_point_variable

.. autofunction:: aronnax.time_series_variable




Coriolis fields
----------------

.. autofunction:: aronnax.f_plane_f_u

.. autofunction:: aronnax.f_plane_f_v

.. autofunction:: aronnax.beta_plane_f_u

.. autofunction:: aronnax.beta_plane_f_v


Domain shape
-------------

.. autofunction:: aronnax.rectangular_pool


