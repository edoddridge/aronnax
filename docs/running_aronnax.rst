Running Aronnax
*******************


To run a simulation with Aronnax one needs to have a Python session active in the folder for the simulation. This folder should contain a file, `aronnax.conf` that contains the configuration choices for the simulation. Some, or all, of these choices may be overridden by arguments passed to the simulate function, but it is likely simpler to specify many of the choices in a configuration file.


.. autofunction:: aronnax.driver.simulate

As described above, it is possible to define functions that can be passed to `aronnax.driver.simulate` and used to create input or forcing fields. The test suite, found in the 'test' folder uses this functionality to create the zonal wind stress for the :math:`\beta`-plane gyre tests. The relevant code is shown below:


   .. code-block:: python

      def wind(_, Y):
          return 0.05 * (1 - np.cos(2*np.pi * Y/np.max(grid.y)))

      with working_directory(p.join(self_path, "beta_plane_gyre_red_grav")):
          drv.simulate(zonalWindFile=wind, valgrind=True,
                       nx=10, ny=10, exe="aronnax_test", dx=xlen/10, dy=ylen/10)
