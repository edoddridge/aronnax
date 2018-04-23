.. role:: bash(code)
   :language: bash

Installation
************************


While we aspire for the installation process for Aronnax to be as simple as :code:`pip install aronnax`, it is not yet that easy.

Dependencies
============

Aronnax is tested with Python 2.7 and 3.6, and depends on several external libraries. 

The Python components of Aronnax depend on 

| :bash:`numpy`
| :bash:`scipy`
|

We recommend installing these with your favourite package manager before installing Aronnax.

Additionally, the Fortran core requires

| :bash:`make`
| :bash:`gfortran >= 4.7.4`
| :bash:`mpi`
| 

In particular, Aronnax will need the :bash:`mpif90` command in order for it to compile its Fortran core. You will need to manually ensure that you have a working installation of each of these.

In addition to these dependencies, the automated tests also depend on :bash:`pytest` and :bash:`matplotlib`.

Installation instructions
=========================

 #. Clone the repository to a local directory

    - :bash:`git clone --recursive https://github.com/edoddridge/aronnax.git`

 #. Move into the base directory of the repository

    - :bash:`cd aronnax`

 #. Compile Hypre

    - move to the directory 'lib/hypre/src'

      - :bash:`cd lib/hypre/src`
    
    - configure the Hypre installer

      - :bash:`./configure`

    - compile Hypre. This will take a few minutes
      
      - :bash:`make install`

    - move back to root directory of the repository

      - :bash:`cd ../../../`

 #. install Aronnax
   
    - :code:`pip install -e ./`

Aronnax is now installed and ready to use. To verify that everything is working, you may wish to run the test suite. Do this by executing :code:`pytest` in the base directory of the repository. This requires that the :bash:`pytest` module is installed.


.. note:: 
    Installing in HPC environments: If your cluster requires programs to be compiled on the compute cores, then you will need to perform step 3 on the compute cores.
