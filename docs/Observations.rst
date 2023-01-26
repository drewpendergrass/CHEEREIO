.. _Observations:

Observations
==========

Observations in CHEEREIO are handled by user-generated observation operators. Because users are most likely to need to modify CHEEREIO in order to introduce new observation types (such as a currently unsupported satellite product), special care has been taken to make it as easy as possible to add new observation operators *without modifying existing CHEEREIO source code*. In this article, we outline how observation operators are implemented in CHEEREIO and describe how users can create their own custom observation operators and plug them in to the CHEEREIO workflow.

What is a CHEEREIO observation operator?
-------------

As described in the :ref:`LETKF technical` section, CHEEREIO works by combining a nested series of objects called Translators. Translators successively abstract raw observational data and GEOS-Chem output and convert it into the quantities needed for LETKF assimilation. Assimilation results are then passed back down through the translators and saved out in a form compatible with GEOS-Chem. CHEEREIO observation operators take (1) GEOS-Chem output data, supplied as a 4D NumPy array (time, level, latitude, longitude) from other CHEEREIO translators; (2) observation data from file; and maps the GEOS-Chem output onto observation space. This often looks like applying a satellite averaging kernel or air mass factors to create a vector of simulated observations, each corresponding to a real observation matched in time and space. This output is then passed to other parts of CHEEREIO for LETKF assimilation. The :ref:`New observation` section walks through in detail how to write a CHEEREIO observation operator.

The observation operator toolkit
-------------

Technically speaking, CHEEREIO observation operators are objects that inherit from :ref:`Observation_Translator`. In object-oriented programming, inheritance can be thought of as a sophisticated form of templating. Indeed, the Observation_Translator class itself is mostly empty, and contains instructions to the user on how to write two standardized methods to (1) read observations from file and process them into a Python dictionary formatted for CHEEREIO, and (2) generate simulated observations :math:`y_i^b` from GEOS-Chem output, returning the results as an ObsData object (:ref:`ObsData`). Users can easily write their own class inheriting from Observation_Translator for a specific use case (like a particular surface or satellite instrument) by implementing these two methods, optionally employing the tools provided in the observation toolkit (detailed below) to make coding easier. The full instructions for adding a new observation are given in :ref:`New observation`.

Any class written with this strict Observation_Translator template will then plug in automatically to the rest of the CHEEREIO workflow; it can then be activated from the main configuration file by any user, not just the original author. CHEEREIO also comes with some pre-written observation operators (such as for the TROPOMI and OMI satellite instruments). Many different observation operators can be used simultaneously, making it natural to perform multispecies data assimilation or assimilation using both surface and satellite data within the CHEEREIO framework. Again, because Observation_Translators handle the details of interpreting a specific observation type, the rest of CHEEREIO can remain ignorant of specifics and operate in a fully abstract environment that can be reused for all simulations.

In the below entries, we will briefly describe the functions and classes supplied in the observation operator toolkit (observation_operators.py).

The produceSuperObservationFunction function
~~~~~~~~~~~~~

This section is under construction, check back later!

In CHEEREIO, users can opt to aggregate observations together to the GEOS-Chem grid, rather than ingesting observations individually even if several separate observations are made in a GEOS-Chem grid cell for a given point in time. Using "super-observations" is helpful because it reduces the complexity of the LETKF calculation and can help control ensemble behavior in areas where there are dense observations. The one tricky bit is handling error -- errors are reduced when aggregating observations. That's where the produceSuperObservationFunction function comes in.

.. py:function:: produceSuperObservationFunction(fname)

   Takes as input a string for the function name. Users supply this string in ``ens_config.json`` via the ``SUPER_OBSERVATION_FUNCTION`` entry. Returns a function as output which will then be used to calculate errors after super-observation aggregation.

   :param str fname: The name of the function. Currently supported values are "default", "sqrt", and "constant". The details of these functions are described in the ``AV_TO_GC_GRID`` entry on the :ref:`Configuration` page.
   :return: The super observation function ``super_obs()``
   :rtype: function
   :raises ValueError: if the function name is unrecognized

Users never use the super observation function that is output by the ``produceSuperObservationFunction`` function directly, but CHEEREIO does. Therefore it is important that the super observation function has a standardized call signature and output. Details on how to right your own function are given in the :ref:`New superobservation` section.

.. py:function:: super_obs(mean_error,num_obs,errorCorr=0,min_error=0,[transportError=0])

   A function which takes the mean error for a set of observations, the number of observations that are averaged, and a few other parameters, and outputs the new error resulting from the aggregation (which is always less than or equal to the mean error input).

   :param mean_error: The mean error for the observations which have been aggregated together. Like all parameters, CHEEREIO supplies this number, calculated based on user settings.
   :ptype: float
   :param num_obs: The number of observations which have been aggregated together.
   :ptype: int
   :param errorCorr: The correlation (between 0 and 1) between individual observations. The function must always have this argument in the call signature even if it is not used.
   :ptype: float
   :param min_error: The error floor for the super observation. Error will never fall below this number. The function must always have this argument in the call signature even if it is not used.
   :ptype: float
   :param transportError: Irreducible error attributable to model transport. Only some super observation functions use this quantity.
   :ptype: float
   :return: The reduced error which will be associated with the super observation in the LETKF calculation
   :rtype: float

The apply_filters function
~~~~~~~~~~~~~

This section is under construction, check back later!

The nearest_loc function
~~~~~~~~~~~~~

This section is under construction, check back later!

The getGCCols function
~~~~~~~~~~~~~

This section is under construction, check back later!

The averageByGC function
~~~~~~~~~~~~~

This section is under construction, check back later!

.. _Observation_Translator:

The Observation_Translator class
~~~~~~~~~~~~~

This section is under construction, check back later!

.. _ObsData:

The ObsData class
~~~~~~~~~~~~~

This section is under construction, check back later!

Existing observation toolkits
-------------

This section is under construction, check back later!

.. _TROPOMI tools:

TROPOMI tools
~~~~~~~~~~~~~

This section is under construction, check back later!

.. _OMI tools:

OMI tools
~~~~~~~~~~~~~

This section is under construction, check back later!


Supplementing an existing observation type
~~~~~~~~~~~~~

This section is under construction, check back later!

.. _New observation:

Workflow to add a new observation operator
-------------

This section is under construction, check back later!

(1) Create a class inheriting from Observation_Translator 
~~~~~~~~~~~~~

This section is under construction, check back later!

(2) Implement getObservations() function 
~~~~~~~~~~~~~

This section is under construction, check back later!

(3) Implement gcCompare() function 
~~~~~~~~~~~~~

This section is under construction, check back later!

(4) Update operators.json
~~~~~~~~~~~~~

This section is under construction, check back later!

(5) Link observational files from ens_config.json
~~~~~~~~~~~~~

This section is under construction, check back later!

(6) [optional] Add observation filters via an extension
~~~~~~~~~~~~~

This section is under construction, check back later!

.. _New superobservation:

(7) [optional] Add a new super observation error function
~~~~~~~~~~~~~

This section is under construction, check back later!

