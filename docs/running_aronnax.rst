Running Aronnax
*******************


To run a simulation with Aronnax one needs to have a Python session active in the folder for the simulation. This folder should contain a file, `aronnax.conf` that contains the configuration choices for the simulation. Some, or all, of these choices may be overridden by arguments passed to the simulate function, but it is likely simpler to specify many of the choices in a configuration file.


.. autofunction:: aronnax.driver.simulate

As described above, it is possible to define functions that can be passed to `aronnax.driver.simulate` and used to create input or forcing fields. The test suite, found in the 'test' folder uses this functionality to create the zonal wind stress for the :math:`\beta`-plane gyre tests. The relevant code is shown below:


   .. code-block:: python

      def wind(_, Y):
          return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))

      with working_directory(p.join(self_path, "beta_plane_gyre_red_grav")):
          drv.simulate(zonalWindFile=[wind],
                       nx=10, ny=10, exe="aronnax_test", dx=xlen/10, dy=ylen/10)


.. warning::
    Parameters cannot be set to 0 or 1 in the call to `drv.simulate` because the Python wrapper gets confused between numerical and logical variables, as described in https://github.com/edoddridge/aronnax/issues/132



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

 - 0: no additional outputs. Output frequency controlled by `DumpFreq` and `AvFreq`
 - 1: output tendencies at frequency given by `DumpFreq`
 - 2: output tendencies and convergence diagnostics from the linear solve at frequency given by `DumpFreq` (not implemented)
 - 3: output convergence diagnostics and tendencies before and after applying some or all of sponges, barotropic correction, winds, and boundary conditions at frequency controlled by `DumpFreq` (not implemented)
 - 4: dump all of the above fields every time step (mostly implemented)
 - 5: dump everything every time step including the two initial RK4 steps (not implemented)

niter0
------
This parameter allows a simulation to be restarted from the given timestep. It requires that the appropriate files are in the 'checkpoints' directory. All parameters, except for the number of grid points in the domain, may be altered when restarting a simulation. This is intended for breaking long simulations into shorter, more manageable chunks, and for running perturbation experiments. 

wetMaskFile
-----------
The wetmask defines which grid points within the computational domain contain fluid. The wetmask is defined on the tracer points, and a value of 1 defines fluid, while a value of 0 defines land. The domain is doubly periodic in `x` and `y` by default. To produce a closed domain the wetmaks should be set to 0 along the edges of the domain. To close the domain it is sufficient to place a strip of land along either the northern or southern boundary and either the western or eastern boundary. You may find it conceptually easier to close both edges.

RelativeWind
------------
If this is false, then the wind input file is given in N m\ :sup:`--2`. If true, then the wind input file is in m s\ :sup:`--1` and a quadratic drag law is used to accelerate the fluid with `Cd` as the quadratic drag coefficient. 
