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

.. _super observation function:

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

.. _The apply filters function:

The apply_filters function
~~~~~~~~~~~~~

CHEEREIO supports real-time filtering of input observations based on user settings; for example, removing observations with high albedo. The ``apply_filters()`` function is the best way to perform this filtering, because it is able to connect seamlessly with the rest of the CHEEREIO codebase. For a more information, see :ref:`Observation filters`.

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

.. py:function:: averageByGC(iGC, jGC, tGC, GC,GCmappedtoobs,obsvals,doSuperObs,superObsFunction=None,other_fields_to_avg=None, prescribed_error=None,prescribed_error_type=None, obsInstrumentError = None, modelTransportError = None, errorCorr = None,minError=None)

   Average observational data and other parameters onto the GEOS-Chem grid. Returns an ObsData object, which is compatible with the rest of the CHEEREIO workflow (i.e. expected output of gcCompare()).

   :param array iGC: Index array output by :py:func:`nearest_loc`
   :param array jGC: Index array output by :py:func:`nearest_loc`
   :param array tGC: Index array output by :py:func:`nearest_loc`
   :param DataSet GC: An xarray dataset which contains the combined GEOS-Chem model output. GC is provided by other CHEEREIO translators; users creating new observation operators can take it as a given.
   :param array GCmappedtoobs: An array of simulated observations, output from an observation operator, of the same length as the observation data.
   :param array obsvals: An array of observation data.
   :param bool doSuperObs: True or False, should we reduce error when we aggregate observations together. Supplied by user settings. 
   :param superObsFunction str: Name of the super observation function, fed into :py:func:`produceSuperObservationFunction`. Supplied by user settings.
   :param dict other_fields_to_avg: Other fields present in the ObsData object can be averaged to the GC grid (e.g. albedo from TROPOMI). Users wishing to do this can supply a dictionary with keys naming the fields in question (must be present in ObsData) and values storing an array of these fields. Note that CHEEREIO automatically will calculate and provide all fields provided in the ``EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y`` setting in ``ens_config.json`` by default without any user code or input required (including for new observation operators).
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

The ObsData class is a simple data storage class, used by CHEEREIO to handled data output from observation operators. It operates very similarly to a dictionary, but forces expected data to be present and allows users to flexibly add additional data (like albedo) which is present for some but not other operators. 

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

All observation operators that are compatible with CHEEREIO must inherit from the Observation_Translator class.  This class itself is basically empty and functions as a template for users to follow in building their own observation operators. 

.. py:class:: Observation_Translator(verbose)

   All observation operators compatible with CHEEREIO inherit from this abstract class and implement the two required methods ``getObservations()`` and ``gcCompare()``. 

   :attr int verbose: Verbosity of output; 1 is default, higher values print more statements during calculation. 
   :attr dict spc_config: Ensemble configuration data, from ``ens_config.json``.
   :attr str scratch: Path to the ensemble scratch folder.

   .. py:method:: Observation_Translator.initialReadDate()

      This is an **optional** function that observation operators can implement, and as such is not present in the abstract class. For a sorted list of all the observation files, indicate in a dictionary the start and end datetimes of the data in each file and save as a pickle file into the ``scratch/`` directory. If an observation operator implements this function, then it should be indicated in operators.json (see :ref:`operators_json`).

      :return: Start and end dates for each observation data file used in assimilation.
      :rtype: dict

   .. py:method:: Observation_Translator.getObservations(specieskey,timeperiod, interval=None, includeObsError=False))

      For a given species and time period of interest, provide a dictionary with all relevant observations and observation metadata (latitude, longitude, time, and others as user needs).

      :param str specieskey: The species of interest. Please note that the "specieskey" variable MUST be a key in the dictionary OBSERVED_SPECIES in ens_config.
      :param list timeperiod: A list of two datetime objects indicating the start and end times of the observations which we need to aggregate.
      :param int interval: If this value is not None, only read observations at a given interval (e.g. only a certain hour a day). This is only used to speed up certain kinds of plots; while users are required to accept this argument as an argument, they can opt to raise an error if interval is not None rather than implement this functionality.
      :param bool includeObsError: True or False, should we load in error data which accompanies observation data? If your data set does not have error data that is useful, you can ignore this or opt to raise an error if this is set to True.
      :return: Observation data formatted as a dictionary. The returned dictionary must have keys for "latitude", "longitude", and "utctime", where UTC time is an ISO 8601 date time string. Actual observation data can be named however the user would like, so long as ``gcCompare()`` can handle it.
      :rtype: dict
      :raises NotImplementedError: if the user fails to implement this function.

   .. py:method:: Observation_Translator.gcCompare(specieskey,OBSDATA,GC,GC_area=None,saveAlbedo=False,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None))

      THE BIG DADDY. This function takes as input observation data (formatted as a dictionary) and GEOS-Chem data from the ensemble (an xarray DataSet), and returns as output an ObsData object with GEOS-Chem data mapped into observation space, along with relevant metadata. This is the heart and soul of the observation operator.

      :param str specieskey: The species of interest. Please note that the "specieskey" variable MUST be a key in the dictionary OBSERVED_SPECIES in ens_config.
      :param dict OBSDATA: Observation data in dictionary form. See :py:func:`apply_filters` for more details.
      :param DataSet GC: An xarray dataset which contains the combined GEOS-Chem model output. GC is provided by other CHEEREIO translators; users creating new observation operators can take it as a given.
      :param array GC_area: If we are using the grid cell areas, we supply them here (rare)
      :param bool saveAlbedo: True or False, should we save out albedo data into the ObsData object we return. Ignore if you do not use.
      :param bool doErrCalc: In the case where ``AV_TO_GC_GRID`` is set to True, should we calculate the aggregated error and save it out in the returned ObsData object? This is always ignored if ``AV_TO_GC_GRID`` is set to False. THIS AND ALL SUBSEQUENT PARAMETERS ARE ONLY USED IF ``AV_TO_GC_GRID`` is set to True.
      :param bool useObserverError: True or False, are we using error accompanying observation files or (if False) are we using prescribed (absolute/relative) errors.
      :param float prescribed_error: If working with prescribed errors, then this is either the percent or absolute error associated with the observational data.
      :param str prescribed_error_type: If using prescribed errors, a string for "relative" or "absolute" denoting how ``prescribed_error`` should be interpreted (as percent or an absolute value). Supplied by user configuration settings.
      :param float transportError: If using a super-observation function that accounts for model transport error, the transport error. Supplied by the user configuration settings.
      :param array errorCorr: : If using a super-observation function that accounts for correlation between errors, the error correlation. Supplied by the user configuration settings.
      :param array minError: If using a super-observation function that accounts for minimum error, the minimum error allowed for a specific observation. Supplied by the user configuration settings.
      :return: ObsData type object containing observation data, relevant metadata (lat/lon/time/etc), and GEOS-Chem data mapped via this function onto observation space.
      :rtype: ObsData
      :raises NotImplementedError: if the user fails to implement this function.

.. _New observation:

Workflow to add a new observation operator
-------------

Because of CHEEREIO's modular design, adding a new observation operator for an arbitrary new observation type (like a new satellite instrument) is straightforward. Here we walk through the process to write a new CHEEREIO observation operator step-by-step. You can always look at existing tools for an additional model to follow (e.g. :ref:`TROPOMI tools` in ``core/tropomi_tools.py`` and :ref:`OMI tools` in ``core/omi_tools.py``).

(1) Create a class inheriting from Observation_Translator 
~~~~~~~~~~~~~

Create a new file in the ``core/`` folder to contain your observation operator and whatever support methods you need to write. At a minimum, you will need to import the ``observation_operators`` module, because all observation operators compatible with CHEEREIO need to inherit from :py:class:`Observation_Translator`.

Suppose we were writing an operator for surface NO2 monitors. Here's how your new ``surface_no2_translator.py`` file might start:

.. code-block:: python

	import observation_operators as obsop

	class Surface_NO2_Translator(obsop.Observation_Translator):
		def __init__(self,verbose=1):
			super().__init__(verbose)

The init function here is run when the Surface_NO2_Translator object is created. The ``super()`` function here just means that we run the default initialization included in the :py:class:`Observation_Translator` class; unless you have a very good reason to do so, you should probably use this initialization.


(2) Implement getObservations() function 
~~~~~~~~~~~~~

The ``getObservations()`` method takes as input (1) a string representing a key from ``OBSERVED_SPECIES`` in ``ens_config.json``, and (2) a pair of datetime objects representing the start and end of the period of interest. This function should return a dictionary of all the relevant observations of the species of interest in this timeperiod. The dictionary should have (1) the observation data stored as a 1D NumPy array, (2) required metadata fields of latitude, longitude, and utc time of the observations each stored as NumPy arrays within the dictionary (see :py:class:`Observation_Translator` for details), and (3) additional metadata fields as necessary for your observation operator, such as albedo.

Suppose all of our observation data is in a CSV file. Our ``getObservations()`` function might look like this (pseudocode, may not be exactly right):


.. code-block:: python

	import observation_operators as obsop
	import pandas as pd

	class Surface_NO2_Translator(obsop.Observation_Translator):

		def __init__(self,verbose=1):
			super().__init__(verbose)

		def getObservations(self,specieskey,timeperiod, interval=None, includeObsError=False):

			#Get the name of the species we are observaing
			species_of_interest = self.spc_config['OBSERVED_SPECIES'][specieskey] 

			#Here we imagine the user specifies the csv file path in the ens_config.json file.
			data_path = self.spc_config['Surface_dirs'][species_of_interest] 

			#In this implementation, the NO2 data must be named NO2.csv.
			data_file = f'{data_path}/NO2.csv

			#Load the data
			data = pd.read_csv(data_file) 

			#subset the data to the right timespan.
			data = data[(data['date'].dt>=timeperiod[0]) & (data['date'].dt<timeperiod[1])] 

			#make an empty dictionary to return, fill it
			to_return = {} 
			to_return['NO2'] = data['NO2'] #store data in the dictionary
			to_return['latitude'] = data['latitude']
			to_return['longitude'] = data['longitude']
			to_return['utctime'] = data['utctime']

			#we only include the error associated with the measurement if requested.
			if includeObsError: 
				to_return['error'] = data['error']

			return to_return #return the data
			
Note that this function is designed to (1) load data from file, (2) subset to the timeperiod of interest, and (3) return in a standardized dictionary form. Notice that the specific CSV file location is loaded from ``ens_config.json`` and stored in the ``self.spc_config`` object; you can always get user settings from ``ens_config.json`` via this object. Although the ``Surface_dirs`` is not in our current ``ens_config.json`` format, new observation operators will require us to grow the configuration file. See :ref:`observation_link` for more details on how to add new setting fields to ``ens_config.json`` for your observation operator.

You may want to allow users to filter out bad observational data, using filters they can set via an extension. You should implement that filtering in this function, by following the :ref:`Observation filters` procedure.

(3) Implement gcCompare() function 
~~~~~~~~~~~~~

The ``gcCompare()`` method takes as input (1) a string representing a key from ``OBSERVED_SPECIES`` in ``ens_config.json``, (2) observational data in dictionary form, as output by getObservations(), (3) GEOS-Chem data as an xarray DataSet, and (4) additional input data, such as that specifying how to handle errors when aggregating to the GEOS-Chem grid. See :py:class:`Observation_Translator` for details. We return data as an ObsData object. Here we show a very simplified version of a gcCompare function (pseudocode, may not be exactly right).

.. code-block:: python

	def gcCompare(self,specieskey,OBSDATA,GC,GC_area=None,saveAlbedo=False,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):

		#Get the name of the species we are observaing
		species = self.spc_config['OBSERVED_SPECIES'][specieskey] 

		#Use one of the utilities included in the Observation Operator toolkit to get GEOS-Chem data matching observations
		GC_col_data = obsop.getGCCols(GC,OBSDATA,species,self.spc_config,returnStateMet=returnStateMet,GC_area=GC_area)

		#The GEOS-Chem 2D array, where the first dimension matches the length of observations and second is the column height.
		GC_SPC = GC_col_data['GC_SPC']

		#Get the surface data
		GC_surf = GC_SPC[:,0]

		toreturn = obsop.ObsData(GC_surf,OBSDATA[species],OBSDATA['latitude'],OBSDATA['longitude'],OBSDATA['utctime'])

		return toreturn

Here all we do is use the :py:func:`getGCCols` function to get GEOS-Chem columns that line up in space and time with our observational data. Then we get the surface data from this column. We create a :py:class:`ObsData` object with all the relevant data and return it.

Note that a real observation operator will need to handle errors in the case that users choose to aggregate to the GEOS-Chem grid. Consult :py:func:`produceSuperObservationFunction`, and existing observation operator toolkits, for guidance on how to implement this.

.. _operators_json:

(4) Update operators.json
~~~~~~~~~~~~~

Now that we have written our observation operator, we have to let CHEEREIO know about it! In the top level directory of CHEEREIO, there is a file called ``operators.json``. Open it and add a new entry. In our Surface NO2 example, we would add the following:
::

	"SURFACE_NO2" : {
		"module_name" : "surface_no2_translator",
		"translator_name" : "Surface_NO2_Translator",
		"implements_initialReadDate" : "False"
	}

Here is the description of each of the three options:


.. option:: module_name
	
	File name of the Python script where the observation operator and support methods lives, without the ``.py`` extension. 

.. option:: translator_name
	
	Name of the class inheriting from ``Observation_Translator``, where the main observation operator lives.

.. option:: implements_initialReadDate
	
	"True" or "False", do you implement the initial read date function in your observation operator? This is fully optional.

Now that you have added your operator to ``operators.json``, users can activate the observation operator by setting values in the key-value pairs within ``OBS_TYPE`` to the name of your operator; in this case, SURFACE_NO2. Here is an example of what ``ens_config.json`` could look like:
::

	"OBS_TYPE" : {
		"NO2_OMI":"OMI",
		"NO2_MONITORS":"SURFACE_NO2"
	},

In this example, we are using NO2 from OMI and from the SURFACE_NO2 operators. CHEEREIO will now use your operator to handle surface NO2!

.. _observation_link:

(5) Link observational files from ens_config.json
~~~~~~~~~~~~~

CHEEREIO needs to know where to look for your observational data files; indeed, when you wrote your ``getObservations()`` function you needed to get a file path from ``ens_config.json``. All we have to do is add an additional dictionary in ``ens_config.json`` which is compatible with how we wrote our observation operator. Here is an example:
::

	"OMI_dirs" : {
		"NO2" : "/n/holylfs05/LABS/jacob_lab/dpendergrass/omi/NO2"
	},
	"Surface_dirs" : {
		"NO2" : "/n/holylfs05/LABS/jacob_lab/dpendergrass/surface/NO2"
	},

Now future users can point CHEEREIO towards their observational data without modifying any code!

.. _Observation filters:

(6) [optional] Add observation filters via an extension
~~~~~~~~~~~~~

CHEEREIO supports real-time filtering of input observations based on user settings; for example, removing observations with high albedo. The ``apply_filters()`` function (documented at :ref:`The apply filters function`)is the best way to perform this filtering, because it is able to connect seamlessly with the rest of the CHEEREIO codebase. 

The ``apply_filters`` function takes (1) a dictionary of observation data ``OBSDATA`` formatted for the standard gcCompare functions in observation operators, and (2) a dictionary of filter information called ``filterinfo``. Here we focus on ``filterinfo``.

``filter_info`` is a dictionary. Keys are specific to a given observation type. By default, these are ``MAIN``, ``OMI_NO2``, and ``TROPOMI_CH4``. Values are a list of filter values, distinct for each observation type. The ``apply_filters`` function checks to see if a given key is present in ``filter_info``; if it is, it parses the list of filter values accordingly. In the OMI NO2 case, it looks for solar zenith angle, cloud radiance fraction, and surface albedo and removes data that do not match these filters.

Filter thresholds are usually supplied by ensemble configuration extensions, discussed here: :ref:`Extensions`. Users supply filter values in an extension, which observation operators then load and pass to the ``apply_filters`` function. In the case of TROPOMI CH4, this works as follows:

.. code-block:: console

    if (self.spc_config['Extensions']['TROPOMI_CH4']=="True") and (self.spc_config['TROPOMI_CH4_FILTERS']=="True"): #Check first if extension is on before doing the TROPOMI filtering
            filterinfo["TROPOMI_CH4"] = [float(self.spc_config['TROPOMI_CH4_filter_blended_albedo']),float(self.spc_config['TROPOMI_CH4_filter_swir_albedo_low']),float(self.spc_config['TROPOMI_CH4_filter_swir_albedo_high']),float(self.spc_config['TROPOMI_CH4_filter_winter_lat']),float(self.spc_config['TROPOMI_CH4_filter_roughness']),float(self.spc_config['TROPOMI_CH4_filter_swir_aot'])]

First, the operator checks that filters are activated, then creates adds an entry to the filterinfo dictionary listing thresholds supplied by the user. Follow this pattern to add your own filters.


.. _postfilter:

(7) [optional] Add ability to filter observations after the operator is applied
~~~~~~~~~~~~~

Sometimes you may be want to filter observations after the observation operator is applied. This occurs, for example, in the :ref:`iasi`; in that case, we want to filter out retrievals outside a certain range only after we replace the prior column with our GC columns. In the case of :ref:`ObsPack tools`, for nitrous oxide we support removing observations exceeding a certain difference with the GC simulation as being likely in error. The issue with this form of filtering is that it can result in different ensemble members having different observation counts. A marginal observation might count as valid for one ensemble member and invalid for another. LETKF requires that all observations across each ensemble member be equal in length (we cannot have ragged arrays). The way CHEEREIO handles this is through **postfiltering**.

The implementation of postfiltering in your gcCompare() function is very simple. **Do not drop observations that do not meet your post-operator filtering condition,** as this may lead to ragged arrays and thus LETKF failure. Instead, calculate a boolean array with the same length as your observation vector. Entries with value True meet your filtering conditions and will be preserved for LETKF calculations, while entries with value False should be removed from the calculation. To your ObsData object that you are returning, use the addData function to add an entry called "postfilter." Here is an example from IASI_tools:

.. code-block:: python

      valid_after_postfilter = (SFm_inv_abs<1.5e16) & np.logical_not((np.abs(toreturn.getDataByKey('HRI'))>1.5) & (toreturn.getObsCol() < 0)) #These points are all valid after threshholding.
      toreturn.addData(postfilter=valid_after_postfilter) 

If you add the postfiltering like so, CHEEREIO will automatically preserve only those entries which have value True across all ensemble members, avoiding the ragged array problem.

.. _New superobservation:

(8) [optional] Add a new super observation error function
~~~~~~~~~~~~~

In CHEEREIO, users can opt to average observations to the GEOS-Chem grid by setting the ``AV_TO_GC_GRID`` entry to be True. Your observation operator should consider supporting this functionality. The use of "super observations" is a useful technique to balance prior and observational errors while also reducing the computational complexity of the optimization (by reducing the size of the observational vectors and matrices in the LETKF calculation). The main subtlety that needs to be handled for super observation aggregation is the adjustment of observational error. By default, users can specify one of several error reduction functions in the ``SUPER_OBSERVATION_FUNCTION`` section:

.. option:: sqrt

   A modified version of the familiar square root law, where if we aggregate :math:`n` observations (indexed by :math:`i`) with errors :math:`\sigma_i` together, the new error is :math:`\bar{\sigma}/\sqrt{n}` where :math:`\bar{\sigma}` is the mean of the :math:`\sigma_i`. The modification accounts for correlations :math:`c` between errors (e.g. due to correlated retrieval errors from shared surface type or similar albedo), and for a user-specified minimum error :math:`\sigma_{\min}`. Thus the equation that is actually applied is given by :math:`\max\left[\left(\bar{\sigma}\cdot\sqrt{\frac{1-c}{n}+c}\right),\sigma_{\min}\right]`. The correlation :math:`c` is taken from ``OBS_ERROR_SELF_CORRELATION`` with default value 0, and the minimum error :math:`\sigma_{\min}` is taken from ``MIN_OBS_ERROR`` with default value 0 (i.e. the normal square root law).

.. option:: default

   As with "sqrt", but with an additional term accounting for the fact that GEOS-Chem transport errors are perfectly correlated. Because perfectly correlated errors are irriducible no matter how many realizations are averaged, the resulting equation is given by :math:`\max\left[\sqrt{\bar{\sigma}^2\cdot\left(\frac{1-c}{n}+c\right)+\sigma_t^2},\sigma_{\min}\right]` where :math:`\sigma_t` is transport error supplied by the "transport_error" entry from ``OTHER_OBS_ERROR_PARAMETERS`` etnry.

.. option:: constant

   No error reduction applied. In other words, no matter how many observations are averaged, this function just returns :math:`\bar{\sigma}!

To add a new super observation function, add an additional entry to ``produceSuperObservationFunction`` in the ``observation_operator.py`` file. This function is described on the :ref:`super observation function` entry. At a minimum, all super observation functions need to accept the following four arguments: (1) ``mean_error``, an array of individual observation error values generated by CHEEREIO; (2) ``num_obs``, an array showing the number of observations in a given grid cell; (3) ``errorCorr``, a float showing the correlation of observations to one another (set as a default of zero); and (4) ``min_error``, a float telling the function the minimum error acceptable in a grid cell (set as a default of zero). You do not need to use all four inputs (see the ``constant`` function) but you do need to accept all four for compatability. Additional inputs can be supplied from ``ens_config.json``, as in the case of the ``default`` function and its need for ``transport_error``.  Additional arguments arguments are supplied by users through the ``OTHER_OBS_ERROR_PARAMETERS`` parameter in ``ens_config.json``; just accept these parameters according to their key. See :ref:`LETKF settings` for an example for ``transport_error``.  

(9) Share your operator with the CHEEREIO community
~~~~~~~~~~~~~

Once you are done with your observation operator and have tested it, make a pull request to the main CHEEREIO git repository. That way the community can make use of your tool and advance science faster!

Document your operator by adding an entry in the :ref:`Existing obs tools` section of the documentation. In your documentation, provide the paper that people should cite if they use your operator, so that you get appropriate credit for your work. Also add your paper to the :ref:`appropriate citations` page.


