.. CHEEREIO documentation master file, created by
   sphinx-quickstart on Fri Aug 27 11:24:47 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CHEEREIO's documentation!
====================================

Welcome to the **CHEEREIO** website and documentation! This site provides a manual for the installation, use, and modification of CHEEREIO for a variety of scientific applications.

The GEOS-Chem CHEmistry and Emissions REanalysis Interface with Observations (CHEEREIO) is a set of Python and shell scripts that support data assimilation and emissions inversions for arbitrary runs of the GEOS-Chem chemical transport model via an ensemble approach (i.e. without the model adjoint). CHEEREIO follows five design principles:

#. **Easy to customize**: Assimilate anything, in any GEOS-Chem configuration or simulation.
#. **Easy to maintain**: Science automatically aligned with latest model version.
#. **Easy to deploy**: One configuration file controls installation and settings
#. **Easy to link observations**: Object-oriented observation operator implementation in Python allows the user to rapidly add new kinds of data with minimal programming required.
#. **Quick to run**: Wall runtime should be comparable to vanilla GEOS-Chem assuming resources are available.

This manual assumes that you are familiar with setting up and running the GEOS-Chem chemical transport model. For more information and background, please see the `GEOS-Chem Wiki <http://wiki.seas.harvard.edu/geos-chem/index.php/Main_Page>`__. Although this manual offers an intuitive explanation of CHEEREIO's data assimilation algorithm in the :ref:`How CHEEREIO works` section, you should probably consider the original algorithm paper `Hunt. et. al., [2007] <https://doi.org/10.1016/j.physd.2006.11.008>`__ as a prerequisite to using CHEEREIO. Don't worry, it's quite readable!

I am more than happy to collaborate with you, however you choose to use CHEEREIO: just send an email to me at pendergrass [AT] g [DOT] harvard [DOT] edu. Questions of any kind are also more than welcome. This can be on anything from interpreting results to getting the ensemble installed on your machine. Bugs should be reported by opening an `issue <https://github.com/drewpendergrass/CHEEREIO/issues>`__ on Github, as this will allow all users to see the solution.

.. toctree::
   :maxdepth: 4   
   :caption: Basics of CHEEREIO

   About-CHEEREIO
   Installing-CHEEREIO
   Overview-Of-Capabilities

.. toctree::
   :maxdepth: 4  
   :caption: Using the ensemble

   The-ensemble-configuration-file
   The-setup-ensemble-script
   Guide-to-the-ensemble-directory
   Running-the-ensemble

.. toctree::
   :maxdepth: 4 
   :caption: The assimilation toolkit

   Field-guide-to-scripts
   The-letkf-utils-module
   Adding-new-observations

.. toctree::
   :maxdepth: 4
   :caption: Postprocessing

   Postprocess-workflow
   Postprocessing-API-guidebook

