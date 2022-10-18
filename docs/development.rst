Development
===========


Since latest release
--------------------

- Add capacity to simulate an "upside down" reduced gravity model with a quiescent abyss and an active lower layer `GH241 <https://github.com/edoddridge/aronnax/pull/241>`_ (18 October 2022)
- Implement snake case for parameter and subroutine names `GH239 <https://github.com/edoddridge/aronnax/pull/239>`_ (14 October 2022)
- Update continuous integration to GitHub actions `GH234 <https://github.com/edoddridge/aronnax/pull/234>`_ (13 October 2022)
- open_mfdataarray bugfix. Was loading copies of last array, now loads all timestamps as expected `GH223 <https://github.com/edoddridge/aronnax/pull/223>`_ (20 March 2019)
- Halos. Improve parallel performance by making better use of tile halos to reduce number of exchange calls `GH220 <https://github.com/edoddridge/aronnax/pull/220>`_ (14 March 2019)
- Benchmarking. Add figure with runtime for 1.5 and 2 layer configurations on 1-64 processors `GH216 <https://github.com/edoddridge/aronnax/pull/216>`_ (13 March 2019)
- Docs improvements by @navidcy and @edoddridge `GH215 <https://github.com/edoddridge/aronnax/pull/215>`_ and GH217 <https://github.com/edoddridge/aronnax/pull/217>`_ (13 March 2019)
- Checkpointing bugfix. Subroutine was using global array indicies for local array `GH214 <https://github.com/edoddridge/aronnax/pull/214>`_ (12 March 2019)

Version 0.3.0 (12 March 2019)
--------------------------------

- MULTICORE Aronnax can now run in parallel across multiple processors using MPI `GH212 <https://github.com/edoddridge/aronnax/pull/212>`_ (12 March 2019)
- Remove unused subroutine for enforcing minimum layer thickness `GH211 <https://github.com/edoddridge/aronnax/pull/211>`_ (9 March 2019)
- Bug fix for bottom drag (didn't apply drag in one-layer simulations) `GH209 <https://github.com/edoddridge/aronnax/pull/209>`_ (9 March 2019)
- Bug fix for open_mfdatarray, now respects variable location on C-grid `GH208 <https://github.com/edoddridge/aronnax/pull/208>`_ (21 December 2018)
- Improved initial guess for external solver routine `GH206 <https://github.com/edoddridge/aronnax/pull/206>`_ (20 November 2018)
- Python functions now have axis order specified by Comodo conventions `GH204 <https://github.com/edoddridge/aronnax/pull/204>`_ (13 July 2018)
- Python wrapper can lazily load multiple timestamps into a 4D array (t,z,y,x) `GH204 <https://github.com/edoddridge/aronnax/pull/204>`_ (13 July 2018)
- Aronnax paper published in Journal of Open Source Software (15 June 2018)

Version 0.2.0 (15 June 2018)
--------------------------------

- Allow outcropping `GH196 <https://github.com/edoddridge/aronnax/pull/196>`_ (30 April 2018)
- Add Python 3 compatibility `GH189 <https://github.com/edoddridge/aronnax/pull/189>`_ (23 April 2018)
- Add Adams-Bashforth family of timestepping algorithms up to fifth-order `GH184 <https://github.com/edoddridge/aronnax/pull/184>`_ (15 March 2018)
- Make test suite plot all differences when a test fails, as suggested in `GH136 <https://github.com/edoddridge/aronnax/issues/136>`_ and implemented in `GH181 <https://github.com/edoddridge/aronnax/pull/181>`_ (12 March 2018)
- Refactor Fortran code into modules in src/ directory `GH181 <https://github.com/edoddridge/aronnax/pull/181>`_ (12 March 2018)
- Add option for first-order upwind advection scheme `GH179 <https://github.com/edoddridge/aronnax/pull/179>`_ (8 March 2018)
- Add vertical diffusion of mass between layers `GH177 <https://github.com/edoddridge/aronnax/pull/177>`_ (2 March 2018)


Version 0.1.0 (11 December 2017)
--------------------------------

Initial release