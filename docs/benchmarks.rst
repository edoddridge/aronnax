Benchmarking
************************

Aronnax runs in two different modes. These two modes have
substantially different runtimes due to the equations being solved.

These benchmarks are run on a single processor.


Reduced gravity mode
========================
The reduced gravity mode is substantially faster than the n-layer mode.

.. figure:: beta_plane_bump_red_grav_scaling.png
   :alt: run time for reduced gravity mode

   Time taken to simulate 502 timesteps in the reduced gravity mode.


n-layer mode
==========================
Because it requires a linear equation solve at every timestep, modeling the ocean
floor is substantially more expensive.

.. figure:: beta_plane_bump_scaling.png
   :alt: run time for n-layer mode

   Time taken to simulate 502 timesteps in the n-layer mode.
