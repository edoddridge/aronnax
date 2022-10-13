Running Aronnax
*******************


To run a simulation with Aronnax one needs to have a Python session active in the folder for the simulation. This folder should contain a file, `aronnax.conf` that contains the configuration choices for the simulation. Some, or all, of these choices may be overridden by arguments passed to the simulate function, but it is likely simpler to specify many of the choices in a configuration file.

.. autofunction:: aronnax.driver.simulate

As described above, it is possible to define functions that can be passed to `aronnax.driver.simulate` and used to create input or forcing fields. The test suite, found in the 'test' folder uses this functionality to create the zonal wind stress for the :math:`\beta`-plane gyre tests. The relevant code is shown below:


   .. code-block:: python

      def wind(_, Y):
          return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))

      with working_directory(p.join(self_path, "beta_plane_gyre_red_grav")):
          drv.simulate(zonal_wind_file=[wind],
                       nx=10, ny=10, exe="aronnax_test", dx=xlen/10, dy=ylen/10)


.. warning::
    Parameters cannot be set to 0 or 1 in the call to `drv.simulate` because the Python wrapper gets confused between numerical and logical variables, as described in https://github.com/edoddridge/aronnax/issues/132

Executables
===========

Aronnax includes multiple targets for the Makefile. These produce executables intended either for testing or production runs. The different options for `exe` are discussed below. For a comparison on the execution speed of these executables see :ref:`benchmarking`.

For production runs, and simulations that span multiple processors, use either `aronnax_core` (for :math:`n+1/2`-layer mode) or `aronnax_external_solver` (for :math:`n`-layer mode).


- `aronnax_core`

  - Uses the internal Fortran code and aggressive compiler optimisation. This is the best executable to use for reduced gravity (:math:`n+1/2` layer) simulations. It may also be used for :math:`n` layer simulations but is considerably slower than `aronnax_external_solver` since it does not use the external Hypre library for the pressure solve. `aronnax_core` cannot run on multiple processors in :math:`n`-layer mode.


- `aronnax_test`
  
  - Uses the internal Fortran code and no compiler optimisations. This is the best executable to use for assessing code coverage when running tests. The Fortran error messages might be helpful.

- `aronnax_prof`

  - Executable intended solely for profiling the Fortran core. Unlikely to be of general use.

- `aronnax_external_solver_test`

  - Uses the external Hypre library to solve the matrix inversion in the :math:`n` layer mode. Uses no optimisations for the internal Fortran code and is intended for assessing code coverage when running tests.

- `aronnax_external_solver`

  - Uses the external Hypre library and aggressive compiler optimisations. This is the fastest executable to use in :math:`n` layer mode. It can also be used in :math:`n+1/2` layer mode, but there is no advantage over `aronnax_core`.



Parameters
===========
Parameters can be passed to the model in two ways. Either they can be included in a file called `aronnax.conf` in the working directory, or they can be passed as keyword arguments to :meth:`aronnax.driver.simulate`. The main directory of the repository includes an example `aronnax.conf` file.

.. note::
    Aronnax takes a deliberately strict approach towards specifying parameter choices. The model contains very few default values, except in situations where a default of zero can sensibly be applied. As such, you will need to specify your parameter choices, either through the configuration file, or the function call. However, parameters that are not required, for example, bottom drag in n+1/2 layer mode, need not be set.

The example file is reproduced here with comments describing the parameters and their units. All possible parameters are included, but they are not all assigned values. After modifying this file for your simulation, any unassigned parameters should be deleted.

.. include:: ../aronnax.conf
   :literal:

.. warning::
    The configuration file shown above includes all of the possible input parameters and fields since it forms part of the documentation. IT WILL NOT WORK AS IS. To use it in an actual simulation the file will need to be modified either by giving values to the parameters that are currently unspecified, or deleting them from the file. If you wish to see a configuration file that corresponds to a successful simulation, look in any of the test, benchmark, or reproduction directories.

debug_level
-----------
This parameter determines whether the model produces additional outputs. It should be set to an integer value greater than or equal to zero. The different values have the following effects:

 - 0: no additional outputs. Output frequency controlled by `dump_freq` and `av_freq`
 - 1: output tendencies at frequency given by `dump_freq`
 - 2: output tendencies and convergence diagnostics from the linear solve at frequency given by `dump_freq` (not implemented)
 - 3: output convergence diagnostics and tendencies before and after applying some or all of sponges, barotropic correction, winds, and boundary conditions at frequency controlled by `dump_freq` (not implemented)
 - 4: dump all of the above fields every time step (mostly implemented)
 - 5: dump everything every time step including the two initial RK4 steps (not implemented)

niter0
------
This parameter allows a simulation to be restarted from the given timestep. It requires that the appropriate files are in the 'checkpoints' directory. All parameters, except for the number of grid points in the domain, may be altered when restarting a simulation. This is intended for breaking long simulations into shorter, more manageable chunks, and for running perturbation experiments. 

wet_mask_file
-----------
The wetmask defines which grid points within the computational domain contain fluid. The wetmask is defined on the tracer points, and a value of 1 defines fluid, while a value of 0 defines land. The domain is doubly periodic in `x` and `y` by default. To produce a closed domain the wetmask should be set to 0 along the edges of the domain. Placing a barrier at either the southern or the northern boundary will remove periodicity in the meridional direction. Similarly, a barrier along either the eastern or western boundary will prevent the model from being periodic in the zonal direction.

relative_wind
------------
If this is false, then the wind input file is given in N m\ :sup:`--2`. If true, then the wind input file is in m s\ :sup:`--1` and a quadratic drag law is used to accelerate the fluid with `Cd` as the quadratic drag coefficient. 

h_advec_scheme
------------
`h_advec_scheme` is an integer that selects the advection scheme used for thickness in the continuity equation. Currently two options are implemented:

- `h_advec_scheme` = 1 (default) uses a first-order centered stencil
- `h_advec_scheme` = 2 uses a first-order upwind stencil

ts_algorithm
------------
`ts_algorithm` is an integer that selects the timestepping algorithm used by the model. The default behaviour is to use a third-order Adams-Bashforth scheme (`ts_algorithm` = 3), with the initialisation performed by a second-order Runge-Kutta method (`ts_algorithm` = 12).

 - ts_algorithm = 1:  Forward Euler
 - ts_algorithm = 2: Second-order Adams-Bashfort
 - ts_algorithm = 3: Third-order Adams-Bashfort (default)
 - ts_algorithm = 4: Fourth-order Adams-Bashfort
 - ts_algorithm = 5: Fifth-order Adams-Bashfort
 - ts_algorithm = 12: Second-order Runge-Kutta
 - ts_algorithm = 13: Third-order Runge-Kutta (not implemented)
 - ts_algorithm = 14: Fourth-order Runge-Kutta (not implemented)


More complicated forcing fields
===============================

It is possible to use surface forcing fields that vary both in space and time.

If the pattern of the surface forcing is always the same, the `wind_mag_time_series_file` parameter can be used to make the magnitude of the forcing vary. The provided `wind_mag_time_series_file` needs an entry for each timestep of the simulation. Like all inputs, this can be specified as a Python function that will be evaluated when the model initialises and the output saved to disk so that the Fortran core can obtain the input.

Spatially varying forcing is slightly more complex. In this situation `zonal_wind_file` and `meridional_wind_file` will now point to three-dimensional arrays. These arrays can be generated by the Python wrapper if a function is passed in the call to :meth:`aronnax.driver.simulate`. This function must depend only on  `X, Y, wind_n_records`. The parameters that must be specified are:

 - `wind_n_records`: the number of time slices in the provided wind field. An integer.
 - `wind_period`: the time (in seconds) between subsequent wind forcing patterns. A float.
 - `wind_loop_fields`: whether the wind forcing should loop back through the provided fields. A Boolean.
 - `wind_interpolate`: whether the wind forcing should be interpolated between slices (provides a continuously varying field), or jump discontinuously between the slices. A Boolean.


Discussion
==========

The following points may be helpful when using the model but don't fit above.

Outcropping
-----------

Aronnax allows for layers to become very thin, which simulates outcropping of isopycnals at the surface, and grounding of isopycnals at the sea floor. Unfortunately outcropping adversely impacts layerwise volume conservation. Over a long simulation the change in volume of a layer may be substantial.

To allow for outcropping the following parameters should be set:

- `hmin` very small
- `wind_depth` set to >10 m
