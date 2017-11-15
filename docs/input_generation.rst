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

Pre-existing Fortran files
++++++++++++++++++++++++++

Fortran unformatted files (often called Fortran binary files) can be used as input. To do this, they should be placed in the 'input' folder, and the name of the file given either in the `aronnax.conf` file, or passed as a keyword argument in the call to `aronnax.driver.simulate`.

The fields need to be the correct shape and size. If they aren't the error messages may be difficult to comprehend or nonexistent depending on whether the field is too big, too small, or the wrong shape. The parameter `DumpWind` can be helpful for ensuring that the wind stress has been set correctly.


Generator functions
+++++++++++++++++++

The use of input generator functions allows the model to be passed user defined functions describing the inputs. Aronnax will evaluate the user defined functions or constants to create the input fields, and save these fields to the 'input' folder in the format required by the Fortran core.


The generator functions can be called in two ways:

 - either directly from the aronnax.conf file, using the syntax `:generator_name:arg1,arg2,...,argn`. In this case the arguments must be numerical (generating very simple layerwise constant value fields).
 - or they can be included in the call to :meth:`aronnax.driver.simulate` using the syntax `field_name=[arg1, arg2, ..., argn]`. In this case the arguments may be either numerical or functions and must be passed as a list. That is, they must be wrapped in square brackets [], even if that list has only one element. Aronnax will check that the length of the list equals the number of layers the field is expected to have and will throw an error if they are not equal.

Regardless of the method used, each argument is evaluated for a single layer. That is, `arg1` is used to create the field for the upper layer. 


For use in `aronnax.driver.simulate` call
-----------------------------------------

When using generator functions to create the inputs through the call to :meth:`aronnax.driver.simulate` it is possible to use Python functions to create the fields. This is an extremely powerful method for creating arbitrary forcings and initial conditions.

If a function is passed to the generator functions it must depend on:

 - `X` and `Y`, created from a `np.meshgrid` call with the appropriate axis, if it is for a 2D field; or
 - `nTimeSteps` and `dt` if it is for a time series.

All other values used within the function must be set in a namespace that the function has access to, or be hard-coded.


For use in `aronnax.conf` file 
------------------------------

These generic generator functions can be used in the `aronnax.conf` file to create simple inputs or initial conditions using the syntax `:generator_name:arg1,arg2,...,argn`. The arguments must be numeric.


For specific grid locations
###########################

These generator functions are passed one numeric argument per layer.

.. autofunction:: aronnax.tracer_point_variable

.. autofunction:: aronnax.u_point_variable

.. autofunction:: aronnax.v_point_variable

.. autofunction:: aronnax.time_series_variable




Coriolis fields
###############

The :math:`f`-plane generator functions are passed one numeric argument, the value for :math:`f`, and the :math:`\beta`-plane functions are passed two numeric arguments.

.. autofunction:: aronnax.f_plane_f_u

.. autofunction:: aronnax.f_plane_f_v

.. autofunction:: aronnax.beta_plane_f_u

.. autofunction:: aronnax.beta_plane_f_v


Domain shape
############

This generator function is passed without arguments in the `aronnax.conf` file.

.. autofunction:: aronnax.rectangular_pool


