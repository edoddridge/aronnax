.. role:: bash(code)
   :language: bash

Installation
************************


While we aspire for the installation process for Aronnax to be as simple as :code:`pip install aronnax`, it is not yet that easy.

Dependencies
============
Aronnax depends on several external libraries, namely, make, gfortran >= 4.7.4, and mpi. This means you will need a working installation of each of these. In particular, Aronnax will need the :bash:`mpif90` command to work in order for it to compile its Fortran core.

Aronnax also depends on numpy and scipy.


Installation instructions
=========================

 #. Clone the repository to a local directory

    - :bash:`git clone --recursive https://github.com/edoddridge/aronnax.git`

 #. Compile Hypre

    - move to the directory 'lib/hypre/src'

      - :bash:`cd lib/hypre/src`
    
    - configure the Hypre installer

      - :bash:`./configure`

    - compile Hypre. This will take a few minutes
      
      - :bash:`make install`

 #. install Aronnax
   
    - :code:`pip install -e .`

Aronnax is now installed and ready to use. To verify that everything is working, you may wish to run the test suite. Do this by executing :code:`pytest` in the base directory of the repository. This requires that the pytest module is installed.

