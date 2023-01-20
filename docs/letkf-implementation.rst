.. _LETKF modules:

The CHEEREIO LETKF implementation
==========

The 4D localized ensemble transform Kalman filter (4D-LETKF) algorithm lies at the heart of CHEEREIO. Everything else merely supports the assimilation operations conducted by this workflow.

Technical overview of the 4D-LETKF implementation
-------------

LETKF assimilation is implemented in CHEEREIO using a structure of nested Python objects, designed primarily to ensure that new observation operators can immediately plug into CHEEREIO and work automatically, without requiring users to have deep knowledge of the CHEEREIO code structure.

CHEEREIO works by creating a suite of objects called translators, which load data from gridded NetCDF files used by GEOS-Chem runs, form one-dimensional ensemble state vectors :math:`x_i^b`and background vectors of simulated observations :math:`y_i^b` that are acceptable to the CHEEREIO LETKF routine, and convert assimilated state vectors back into a format acceptable to GEOS-Chem. Translator objects are assembled in a nested structure, with low-level translators performing IO operations and basic calculations to form vectors, which are then passed to objects that operate at a higher level of abstraction. Abstract objects do the actual LETKF calculations without any knowledge of the GEOS-Chem simulation or even the user-defined rules on how to construct the state vector, enabled by the fully general nature of the LETKF. Because all the details of a specific simulation are handled by low-level translators, which are designed to easily expand to include new capabilities added by the community, users are able to modify only one small part of CHEEREIO and remain assured that workflow will continue to operate as expected.

Now for the details.

State vectors are handled by an array of GC_Translator objects (:ref:`GC Translator`), one per ensemble member, which wrap around the GEOS-Chem restart file and scaling factors. These objects translate between gridded NetCDF files used by GEOS-Chem and the one-dimensional vectors needed by the LETKF process. The Assimilator object (:ref:`Assimilator`) contains an array of GC_Translator objects and uses them to form the background perturbation matrix :math:`X_i^b`, a key input into the LETKF algorithm.

To facilitate the 4D-LETKF calculation, the HIST_Translator object (:ref:`HIST Translator`, one per ensemble member) wraps around GEOS-Chem species concentration output files along with other files necessary for calculating simulated observations, such as meteorological information. Because the 4D version of LETKF requires aligning observations with the model state in time, the HIST_Translator object accesses model data from throughout the assimilation window. HIST_Translator objects are stored in an array within a HIST_Ens object (:ref:`HIST Ensemble`, one per Assimilator object), which is responsible for aligning observations with simulated observations obtained by applying an observation operator to model output. 

The HIST_Ens object handles observations by using objects inheriting from the Observation_Translator class (details on the :ref:`Observations` page). In object-oriented programming, inheritance can be thought of as a sophisticated form of templating. The Observation_Translator class is a template that contains instructions to the user on how to write two standardized methods to (1) read observations from file and process them into a Python dictionary formatted for CHEEREIO, and (2) generate simulated observations :math:`y_i^b` from GEOS-Chem output. Users can easily write their own class inheriting from Observation_Translator for a specific use case (like a particular surface or satellite instrument) by implementing these two methods, optionally employing a provided observation toolkit. Any class written with this strict template will then plug in automatically to the rest of the CHEEREIO workflow, and can be activated like usual from the main configuration file. CHEEREIO also comes with some pre-written observation operators (such as one for the TROPOMI satellite instrument). Many different observation operators can be used simultaneously, making it natural to perform multispecies data assimilation or assimilation using both surface and satellite data within the CHEEREIO framework. The HIST_Ens object processes this information into a simulated observation perturbation matrix :math:`Y_i^b`, which it passes to the Assimilator object. The Assimilator object combines all of the data necessary for the assimilation and actually calculates the LETKF calculation. 


The LETKF classes
-------------

In this section, we cover the main objects involved in the LETKF assimilation process, as outlined in the above technical overview. Users will rarely have to work directly with the classes described on this page (except :ref:`Observations`, which have their own page) -- this page is aimed at developers and users who are seeking to understand CHEEREIO bugs.

.. _GC Translator:

The GC Translator class
~~~~~~~~~~~~~

CHEEREIO uses the GC_Translator object to handle data from (1) GEOS-Chem restart files, which are used by CHEEREIO to represent the concentration components of the state vector, and (2) gridded NetCDF files representing emissions scaling factors, which are fed to GEOS-Chem via HEMCO and are used by CHEEREIO to represent the emissions component of the state vector. The GC_Translator is a two-way translator: it contains tools for converting gridded NetCDF files to Python arrays expected by CHEEREIO, and converts assimilated data fields that CHEEREIO calculates in Python arrays into NetCDF files that feed back into GEOS-Chem for the next simulation cycle. At each assimilation cycle, CHEEREIO creates one GC_Translator object for each GEOS-Chem ensemble run directory. The GC_Translator objects are stored in an array within the Assimilator (:ref:`Assimilator`) object.

The GC_Translator object contains (1) an object of class DataBundle (more below); (2) an object of class StateVector (more below); (3) a method to reconstruct NetCDF restart and emissions scaling factor files when given an assimilated ensemble of LETKF state vectors, which includes separation of concentrations and emissions that LETKF requires to be concatenated in a single state vector; and (4) methods to overwrite restart and emissions scaling factor files for the next GEOS-Chem cycle.

The DataBundle is a simple database class that stores all of the restart and emissions scaling factor data that the GC_Translator object uses. This data is stored separately from the GC_Translator object because the StateVector object will also be referencing and modifying data within DataBundle; a dedicated object ensures that all information is synchronized. Beyond storing data, the DataBundle contains a suite of functions for getting and setting data. DataBundle data should *never* be referenced directly; the get and set methods should be used instead, in case additional checks need to be added later. The DataBundle also is able to add assimilated emissions scaling factor data at a new timestamp, which is used by CHEEREIO to record emissions timeseries output.

The StateVector class converts the data in DataBundle into a one dimensional vector used by the LETKF routine. This includes (1) mapping 3D concentrations from the restart file into other representations, like a column sum, and (2) concatenating scaling factors together and combining them with the concentration representation. The class also contains methods for rapidly subsetting sections of the state vector for assimilation calculations, such as localizing the state vector within x kilometers of a given pixel. 

The DataBundle and StateVector classes are only used by the GC_Translator object, not anywhere else in the CHEEREIO codebase.


.. _HIST Translator:

The HIST Translator class
~~~~~~~~~~~~~

CHEEREIO uses the HIST_Translator object to read in GEOS-Chem output data, specifically the time-evolution of chemical species of interest in 3D, and pass them on to :ref:`HIST Ensemble` where they will be used to form the simulated observation perturbation matrix :math:`Y_i^b`. One HIST_Translator object is produced for each ensemble member at assimilation time and reads data from the OutputDir folder in the run directory. This class is very simple; it mainly just concatenates files using xarray and contains some utilities for subsetting data. :ref:`HIST Ensemble` does all the hard work.

.. _HIST Ensemble:

The HIST Ensemble class
~~~~~~~~~~~~~

CHEEREIO uses HIST_Ensemble object to (1) gather the raw GEOS-Chem output data as aggregated by :ref:`HIST Translator`, (2) gather raw observation data, via the objects detailed in the :ref:`Observations` article, and (3) use all the input data and objects supplied in the first two steps to transform GEOS-Chem output data into observation space. The HIST_Ensemble object then uses this data aggregated to form a variety of quantities that are necessary for the LETKF calculation, including the simulated observation mean :math:`\bar{y}^b`, the simulated observation perturbation matrix :math:`Y^b`, the difference between actual and (mean) simulated observations :math:`y-\bar{y}^b`, and the observational error covariance matrix :math:`R`. HIST_Ensemble also takes care of localizing these quantities in the LETKF calculation and applying user-supplied corrections as they are relevant, such as the regularization factor :math:`\gamma`. 

All this data is calculated and supplied as requested by :ref:`Assimilator`. There is only one HIST_Ensemble class produced per assimilator; it holds and processes the data for the entire ensemble.

The Observation Translator class type
~~~~~~~~~~~~~

Please see :ref:`Observations` for more information. 

.. _GT Container:

The GT Container class
~~~~~~~~~~~~~

This section is under construction, check back later!

.. _Assimilator:

The Assimilator class
~~~~~~~~~~~~~

This section is under construction, check back later!
