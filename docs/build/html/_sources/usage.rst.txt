Usage
=====

Installation
------------

Installing genereporter is straight forward--just clone the repository and install the package using pip. 

First, clone the github repository.

.. code-block:: console

   $ git clone https://github.com/chickaroo/genereporter

Next, you might want to make a conda environment to install the package in.

.. code-block:: console

   $ conda create -n genereporter python=3.11
   $ conda activate genereporter

Then, navigate to the directory and install the package using pip.

.. code-block:: console

   $ cd genereporter
   $ pip install .

This will install the package and all of its dependencies using pip. The standard dependencies are
only geared towards the core functionality of the package. If you want to generate any of the GRN data
using `SCENIC <https://scenic.aertslab.org/>`_, for example, you will need to install additional dependencies according to the pySCENIC `documentation <https://pyscenic.readthedocs.io/en/latest/installation.html>`_.

You can now test if the package is installed correctly by running the following command:

.. code-block:: console

   $ python
   >>>import genereporter

If you don't get any errors, the package is installed correctly.


Obtaining Data
--------------

The genereporter repository does not house the VEO-IBD data for storage and safety reasons. 
The data can be obtained from the Helmholtz HPC cluster. 


Vignettes
----------

Check out the notebooks for examples on how to use the package:
:doc:`Cell_Example` and :doc:`GRN_Example`