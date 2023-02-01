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

   :param float mean_error: The mean error for the observations which have been aggregated together. Like all parameters, CHEEREIO supplies this number, calculated based on user settings. 
   :param int num_obs: The number of observations which have been aggregated together.
   :param float errorCorr: The correlation (between 0 and 1) between individual observations. The function must always have this argument in the call signature even if it is not used.
   :param float min_error: The error floor for the super observation. Error will never fall below this number. The function must always have this argument in the call signature even if it is not used.
   :param float transportError: Irreducible error attributable to model transport. Only some super observation functions use this quantity.
   :return: The reduced error which will be associated with the super observation in the LETKF calculation
   :rtype: float

The apply_filters function
~~~~~~~~~~~~~

CHEEREIO supports real-time filtering of input observations based on user settings; for example, removing observations with high albedo. The ``apply_filters()`` function is the best way to perform this filtering, because it is able to connect seamlessly with the rest of the CHEEREIO codebase. For a tutorial on how to add filters for your own observation operator, or on how to add additional filters to existing observation operators, see :ref:`Observation filters`.

.. py:function:: apply_filters(OBSDATA,filterinfo)

   Takes raw observations as read from file and a dictionary-based description of how to filter out bad observations, and returns the filtered raw observations with bad data removed.

   :param dict OBSDATA: Raw observation data in dictionary form. The OBSDATA dictionary must have a numpy array labeled "longitude", one labeled "latitude" and one labeled "utctime". The "utctime" entry must be formatted by ISO 8601 data time format. ISO 8601 represents date and time by starting with the year, followed by the month, the day, the hour, the minutes, seconds and milliseconds. For example, 2020-07-10 15:00:00.000, represents the 10th of July 2020 at 3 p.m. Timezone is assumed to be UTC. For a good example on how to format OBSDATA, see the ``read_tropomi()`` or ``read_omi()`` functions in tropomi_tools and omi_tools respectively. 
   :param dict filterinfo: A dictionary describing the filters to be applied. The keys of the dictionary are called filter families and describe the observation type (e.g. OMI_NO2), while the value is a list with filter values specific to that observation type. For example, in the OMI_NO2 filter family, the solar zenith angle filter value is the first entry in the list.
   :return: The post-filtering version of the input dictionary ``OBSDATA``.
   :rtype: dict

The nearest_loc function
~~~~~~~~~~~~~

CHEEREIO uses the ``nearest_loc()`` function to match observation data with the GEOS-Chem grid. The corresponding index lists are used to (1) get GEOS-Chem columns corresponding with observations, and (2) aggregate observation data to the GEOS-Chem grid.

.. py:function:: nearest_loc(GC,OBSDATA)

   Find the GEOS-Chem grid box and time indices which best correspond with real observation data.

   :param DataSet GC: An xarray dataset which contains the combined GEOS-Chem model output. GC is provided by other CHEEREIO translators; users creating new observation operators can take it as a given. 
   :param dict OBSDATA: Observation data in dictionary form. See :py:func:`apply_filters` for more details.
   :return: Three NumPy arrays ``iGC``, ``jGC``, and ``tGC`` containing the spatial and temporal indices on the GEOS-Chem grid which match observations.
   :rtype: List of NumPy arrays

The getGCCols function
~~~~~~~~~~~~~

Observation operators commonly use the getGCCols function to grab GEOS-Chem columns corresponding in space and time with observations. It returns data as a set of 2D arrays (with dimension index of observation, column) which can be used by subsequent observation operator functions. 

.. py:function:: getGCCols(GC,OBSDATA,species,spc_config,returninds=False,returnStateMet=False,GC_area=None)

   Takes aggregated GEOS-Chem data and grabs columns corresponding with observations in space and time, stored as 2D arrays (with dimension index of observation, column). These GEOS-Chem columns are commonly used by observation operators for calculating simulated observations.

   :param DataSet GC: An xarray dataset which contains the combined GEOS-Chem model output. GC is provided by other CHEEREIO translators; users creating new observation operators can take it as a given. 
   :param dict OBSDATA: Observation data in dictionary form. See :py:func:`apply_filters` for more details.
   :param str species: The species of interest, e.g. "CH4"
   :param dict spc_config: The CHEEREIO ensemble configuration data, stored as a dictionary. This is provided by CHEEREIO.
   :param bool returninds: True or False, should we return the three index fectors from :py:func:`nearest_loc` in addition to the GEOS-Chem columns.
   :param bool returnStateMet: True or False, are we also going to subset meteorological data accompanying GEOS-Chem.
   :param array GC_area: If we are using the grid cell areas, we supply them here (rare)
   :return: A dictionary containing (1) GC_SPC, a 2D array containing GEOS-Chem columns of the species of interest corresponding with observations; (2) GC_P, a 2D array with GEOS-Chem pressure levels; (3) additional entries depending on if returnStateMet is True and/or GC_area is supplied. Indices can also be returned as a list if returninds is True.
   :rtype: dict


The averageByGC function
~~~~~~~~~~~~~

CHEEREIO allows users to specify if they want observations to be aggregated to the GEOS-Chem grid (specified by setting ``AV_TO_GC_GRID`` to ``True`` in ``ens_config.json``). The :py:func:`averageByGC` function ensures that observational data are aggregated onto the GEOS-Chem grid; it also returns results in an ObsData object, which is the expected return type for the gcCompare() function in an observation operator.  

The :py:func:`averageByGC` function accepts two kinds of errors --- prescribed errors (relative or absolute) or individual observational errors supplied with observation files. 

.. py:function:: averageByGC(iGC, jGC, tGC, GC,GCmappedtoobs,obsvals,doSuperObs,superObsFunction=None,albedo_swir=None,albedo_nir=None,blended_albedo=None, prescribed_error=None,prescribed_error_type=None, obsInstrumentError = None, modelTransportError = None, errorCorr = None,minError=None)

   Average observational data and other parameters onto the GEOS-Chem grid. Returns an ObsData object, which is compatible with the rest of the CHEEREIO workflow (i.e. expected output of gcCompare()).

   :param array iGC: Index array output by :py:func:`nearest_loc`
   :param array jGC: Index array output by :py:func:`nearest_loc`
   :param array tGC: Index array output by :py:func:`nearest_loc`
   :param DataSet GC: An xarray dataset which contains the combined GEOS-Chem model output. GC is provided by other CHEEREIO translators; users creating new observation operators can take it as a given.
   :param array GCmappedtoobs: An array of simulated observations, output from an observation operator, of the same length as the observation data.
   :param array obsvals: An array of observation data.
   :param bool doSuperObs: True or False, should we reduce error when we aggregate observations together. Supplied by user settings. 
   :param superObsFunction str: Name of the super observation function, fed into :py:func:`produceSuperObservationFunction`. Supplied by user settings.
   :param array albedo_swir: An optional array of short wave infrared albedo from observations, of the same length as obsvals.
   :param array albedo_nir: As with albedo_swir, but for near wave infrared albedo.
   :param array blended_albedo: As with albedo_swir, but for blended albedo.
   :param float prescribed_error: If working with prescribed errors, then this is either the percent or absolute error associated with the observational data.
   :param str prescribed_error_type: If using prescribed errors, a string for "relative" or "absolute" denoting how ``prescribed_error`` should be interpreted (as percent or an absolute value). Supplied by user configuration settings.
   :param array obsInstrumentError: If working with observational errors, an array of errors associated with each observation.
   :param float modelTransportError: If using a super-observation function that accounts for model transport error, the transport error. Supplied by the user configuration settings.
   :param array errorCorr: : If using a super-observation function that accounts for correlation between errors, the error correlation. Supplied by the user configuration settings.
   :param array minError: If using a super-observation function that accounts for minimum error, the minimum error allowed for a specific observation. Supplied by the user configuration settings.
   :return: An ObsData object containing the aggregated observations.
   :rtype: ObsData
   :raises ValueError: if the error information is specified incorrectly.

.. _ObsData:

The ObsData class
~~~~~~~~~~~~~

This section is under construction, check back later!

.. py:class:: ObsData(gccol,obscol,obslat,obslon,obstime,**additional_data)

   A simple class for storing labeled data output by observation operators. 

   :var array gccol: A NumPy array which contains GEOS-Chem model output mapped onto corresponding observations (i.e. passed through an operation operator. This is sometimes a 2D array containing all ensemble columns. 
   :var array obscol: An array containing observational data.
   :var array obslat: An array containing latitude of observations.
   :var array obslon: An array containing longitude of observations.
   :var array obstime: An array containing times of observations.
   :var dict additional_data: Additional data supplied to the constructor as keyword arguments are stored in dictionary form in this variable.

   .. py:method:: ObsData.getGCCol()

      Returns the gccol attribute

      :return: The gccol attribute of ObsData.
      :rtype: array

   .. py:method:: ObsData.setGCCol(gccol)

      Sets the gccol attribute to the input argument

      :param array gccol: Set ``gccol`` attribute to this value.

   .. py:method:: ObsData.getObsCol()

      Returns the obscol attribute

      :return: The obscol attribute of ObsData.
      :rtype: array

   .. py:method:: ObsData.setObsCol(gccol)

      Sets the obscol attribute to the input argument

      :param array obscol: Set ``obscol`` attribute to this value.

   .. py:method:: ObsData.getCols()

      Returns the gccol and obscol attributes as a list

      :return: The gccol and obscol attributes of ObsData as a list in that order.
      :rtype: list

   .. py:method:: ObsData.getLatLon()

      Returns the obslat and obslon attributes as a list

      :return: The obslat and obslon attributes of ObsData as a list in that order.
      :rtype: list

   .. py:method:: ObsData.getObsTime()

      Returns the obstime attribute

      :return: The obstime attribute of ObsData.
      :rtype: array

   .. py:method:: ObsData.addData(**data_to_add)

      Add custom data fields to ObsData

      :param \**data_to_add: Add data by keyword argument to the ``additional_data`` dictionary attribute as key-value pairs.

   .. py:method:: ObsData.getDataByKey(key)

      Add custom data fields to ObsData

      :param str,list key: Get data from the ``additional_data`` dictionary attribute by key. If supplying a list of keys, it will return a list of values corresponding to each key.
      :return: The array requested by key, or a list of arrays requested by key, from the ``additional_data`` dictionary attribute.
      :rtype: array,list



.. _Observation_Translator:

The Observation_Translator class
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

.. _Observation filters:

(6) [optional] Add observation filters via an extension
~~~~~~~~~~~~~

This section is under construction, check back later!

.. _New superobservation:

(7) [optional] Add a new super observation error function
~~~~~~~~~~~~~

This section is under construction, check back later!

