.. CHEEREIO documentation master file, created by
   sphinx-quickstart on Fri Aug 27 11:24:47 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CHEEREIO's documentation!
====================================

Welcome to the **CHEEREIO** ReadTheDocs documentation! This site provides a manual for the installation, use, and modification of CHEEREIO for a variety of scientific applications.

The GEOS-Chem CHEmistry and Emissions REanalysis Interface with Observations (CHEEREIO) is a set of Python and shell scripts that support data assimilation and emissions inversions for arbitrary runs of the GEOS-Chem chemical transport model via an ensemble approach (i.e. without the model adjoint). CHEEREIO follows five design principles:

#. **Easy to customize**: Assimilate anything, in any GEOS-Chem configuration or simulation.
#. **Easy to maintain**: Science automatically aligned with latest model version.
#. **Easy to deploy**: One configuration file controls installation and settings
#. **Easy to link observations**: Object-oriented observation operator implementation allows the user to rapidly add new kinds of data with minimal programming required.
#. **Quick to run**: Wall runtime should be no more than 2x longer than vanilla GEOS-Chem (4D-Var limit) assuming resources are available.

This manual assumes that you are familiar with the basics of setting up and running the GEOS-Chem chemical transport model. For more information, please see the `GEOS-Chem Wiki <http://wiki.seas.harvard.edu/geos-chem/index.php/Main_Page>`__.

.. toctree::
   :maxdepth: 3   
   :caption: Basics of CHEEREIO

   About-CHEEREIO
   Installing-CHEEREIO
   Overview-Of-Capabilities

.. toctree::
   :maxdepth: 3   
   :caption: Using the ensemble

   Deploying-the-ensemble
   Guide-to-the-ensemble-directory
   Running-the-ensemble

.. toctree::
   :maxdepth: 3   
   :caption: The assimilation toolkit

   The-letkf-utils-module
   The-observation-operator-module

.. toctree::
   :maxdepth: 3   
   :caption: Postprocessing

   Combining-ensemble-runs
   Plotting-with-GCPy

