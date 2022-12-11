.. _LETKF modules:

The CHEEREIO LETKF implementation
==========

The 4D localized ensemble transform Kalman filter (4D-LETKF) algorithm lies at the heart of CHEEREIO. Everything else merely supports the assimilation operations conducted by this workflow.

Technical overview of the 4D-LETKF implementation
-------------

LETKF assimilation is implemented in CHEEREIO using a structure of nested Python objects, designed primarily to ensure that new observation operators can immediately plug into CHEEREIO and work automatically, without requiring users to have deep knowledge of the CHEEREIO code structure. We use Python because of its familiarity to a broad user base, because of its ease of use, and because the object-oriented structure of the language makes it well suited to the modular design of CHEEREIO.

CHEEREIO works by creating a suite of objects called translators, which load data from gridded NetCDF files used by GEOS-Chem runs, form one-dimensional ensemble state vectors :math:`x_i^b`and background vectors of simulated observations :math:`y_i^b` that are acceptable to the CHEEREIO LETKF routine, and convert assimilated state vectors back into a format acceptable to GEOS-Chem. Translator objects are assembled in a nested structure, with low-level translators performing IO operations and basic calculations to form vectors, which are then passed to objects that operate at a higher level of abstraction. Abstract objects do the actual LETKF calculations without any knowledge of the GEOS-Chem simulation or even the user-defined rules on how to construct the state vector, enabled by the fully general nature of the LETKF. Because all the details of a specific simulation are handled by low-level translators, which are designed to easily expand to include new capabilities added by the community, users are able to modify only one small part of CHEEREIO and remain assured that workflow will continue to operate as expected.

Now for the details. State vectors are handled by an array of GC_Translator objects, one per ensemble member, which wrap around the GEOS-Chem restart file and scaling factors. These objects translate between gridded NetCDF files used by GEOS-Chem and the one-dimensional vectors needed by the LETKF process. The Assimilator object, which we will discuss more soon, contains an array of GC_Translator objects and uses them to form the background perturbation matrix :math:`X_i^b`, a key input into the LETKF algorithm.

To facilitate the 4D-LETKF calculation, the HIST_Translator object (one per ensemble member) wraps around GEOS-Chem species concentration output files along with other files necessary for calculating simulated observations, such as meteorological information. Because the 4D version of LETKF requires aligning observations with the model state in time, the HIST_Translator object accesses model data from throughout the assimilation window. HIST_Translator objects are stored in an array within a HIST_Ens object (one per Assimilator object), which is responsible for aligning observations with simulated observations obtained by applying an observation operator to model output. The HIST_Ens object handles observations by using objects inheriting from the Observation_Translator class. In object-oriented programming, inheritance can be thought of as a sophisticated form of templating. As shown in Figure 4, the Observation_Translator class is a template that contains instructions to the user on how to write two standardized methods to (1) read observations from file and process them into a Python dictionary formatted for CHEEREIO, and (2) generate simulated observations :math:`y_i^b` from GEOS-Chem output. Users can easily write their own class inheriting from Observation_Translator for a specific use case (like a particular surface or satellite instrument) by implementing these two methods, optionally employing a provided observation toolkit. Any class written with this strict template will then plug in automatically to the rest of the CHEEREIO workflow, and can be activated like usual from the main configuration file. CHEEREIO also comes with some pre-written observation operators (such as one for the TROPOMI satellite). Many different observation operators can be used simultaneously, making it natural to perform multispecies data assimilation or assimilation using both surface and satellite data within the CHEEREIO framework. The HIST_Ens object processes this information into a simulated observation perturbation matrix :math:`Y_i^b`, which it passes to the Assimilator object. The Assimilator object combines all of the data necessary for the assimilation and actually calculates the LETKF calculation. 


The LETKF classes
-------------

This section is under construction, check back later!

.. _GC Translator:

The GC Translator class
~~~~~~~~~~~~~

This section is under construction, check back later!

.. _HIST Translator:

The HIST Translator class
~~~~~~~~~~~~~

This section is under construction, check back later!

.. _HIST Ensemble:

The HIST Ensemble class
~~~~~~~~~~~~~

This section is under construction, check back later!

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