.. _Observations:

.. _Existing obs tools:

Existing observation operators
==========

CHEEREIO ships with observation operators that have produced by the community. Users can add their own operators by following the instruction on the :ref:`New observation` page, or they can use the operators listed below. 


.. _TROPOMI tools:

TROPOMI tools
-------------

This section is under construction, check back later!

.. _OMI tools:

OMI tools
-------------

NASA's Ozone Monitoring Instrument (OMI) on board the Aura satellite measures criteria air pollutants and other trace gases of interest to CHEEREIO users. Currently, NO\ :sub:`2`\ is the only OMI operator written for CHEEREIO, but users can follow the pattern in ``omi_tools.py`` to add support for additional species.

To activate OMI observations, list "OMI" as an observation type in the ``OBS_TYPE`` setting, as described on the :ref:`Observation settings` page.

OMI fields available to save and plot in CHEEREIO
~~~~~~~~~~~~~

All observation operators save standard data from observations and GEOS-Chem and pass it on to CHEEREIO (including latitude, longitude, observation values, GC simulated observation values, and timestamps). However, individual operators might also be able to save additional data, such as albedo in remote sensing cases or site IDs in surface observation cases. Users can list the additional data they would like to save in the ``EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y`` setting in ``ens_config.json``, following the instructions on the :ref:`postprocessing settings` page. 

Currently, **no additional ObsData fields are available to be saved and/or plotted with the OMI operator** as written. Additional fields can be added by saving more metadata into the ObsData object with the ``addData`` function within the ``gcCompare`` method of the ``OMI_Translator`` class. See the :ref:`Observations` page for more details.

OMI operator support functions
~~~~~~~~~~~~~

The OMI observation operator calls two utility functions to process observations from file and format them in such a way that the ``gcCompare`` method of the ``OMI_Translator`` class can perform relevant computations. They are documented below:

.. py:function:: read_omi(filename, species, filterinfo=None, includeObsError = False)

   A utility function which loads OMI observations from file, filters them, and returns a dictionary of important data formatted for input into the `gcCompare`` method of the ``OMI_Translator`` class.

   :param str filename: NetCDF file containing OMI observations to be loaded. Expects standard level 2 data. 
   :param str species: Name of species to be loaded. Currently, only "NO2" is supported.
   :param dict filterinfo: A dictionary of information about data filtering which is passed to a standard observation operator utility function. See :ref:`Observation filters` for more information
   :param bool includeObsError: True or False, read the errors associated with individual observations. 
   :return: A dictionary containing observation values and metadata, ready for input into the `gcCompare`` method of the ``OMI_Translator`` class.
   :rtype: dict

.. py:function:: clearEdgesFilterByQAAndFlatten(met)

   A utility function takes in partially formatted OMI data and does additional processing, outputting a flattened set of arrays which are compatible with CHEEREIO. In the process, the function removes swath edges and bad retrieval values.

   :param dict met: A dictionary with keys naming important observation data and metadata, and values of raw 2D swath data from OMI. 
   :return: A dictionary containing flattened observation values and metadata, with bad data removed, ready for input into the `gcCompare`` method of the ``OMI_Translator`` class.
   :rtype: dict

.. _ObsPack tools:

ObsPack tools
-------------

`ObsPack <https://doi.org/10.5194/essd-6-375-2014>`__is a standardized dataset containing measurements from  surface monitors distributed around the world, aimed at carbon cycle studies. CHEEREIO users often use ObsPack data for validation of CO or CH4 inversions, or as observations for the inversion itself. The ObsPack observation operator, contained in the ``obspack_tools.py`` file, wraps around the the `ObsPack diagnostic <https://geos-chem.readthedocs.io/en/stable/gcclassic-user-guide/obspack.html>`__ produced by GEOS-Chem and translates it into a form acceptable to CHEEREIO.

To activate ObsPack, see the :ref:`Observation settings` page for information on the correct settings for  ``ens_config.json``. 

ObsPack fields available to save and plot in CHEEREIO
~~~~~~~~~~~~~

All observation operators save standard data from observations and GEOS-Chem and pass it on to CHEEREIO (including latitude, longitude, observation values, GC simulated observation values, and timestamps). However, individual operators might also be able to save additional data, such as albedo in remote sensing cases or site IDs in surface observation cases. Users can list the additional data they would like to save in the ``EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y`` setting in ``ens_config.json``, following the instructions on the :ref:`postprocessing settings` page. Below are a list of the supported fields which users can save in ObsPack:

.. option:: utc_conv
   
   A conversion constant that can change UTC timestamps into local time.

.. option:: altitude
   
   Altitude of ObsPack site.

.. option:: pressure
   
   Pressure observed at ObsPack site.

.. option:: obspack_id
   
   Unique identifier of obspack observation.

.. option:: platform
   
   Obspack platform.

.. option:: site_code
   
   Unique identifier of obspack site.

Users wishing to aggregate ObsPack results by site, either for plotting or other analysis, will want to save ``site_code``. Follow the instructions on the :ref:`postprocessing settings` to ensure this is done successfully.

Additional fields can be added by saving more metadata into the ObsData object with the ``addData`` function within the ``gcCompare`` method of the ``ObsPack_Translator`` class. See the :ref:`Observations` page for more details.

ObsPack preprocessing functions
~~~~~~~~~~~~~

CHEEREIO has built in preprocessing functions to translate raw ObsPack data, as downloaded from NOAA, into a form compatible with the GEOS-Chem ObsPack diagnostic. To use this functionality, set ``preprocess_raw_obspack_files`` to ``true`` in ``ens_config.json`` and provide a path to the raw files in the ``raw_obspack_path`` entry. However, some users report that NOAA ObsPack data is not quite standardized. If you run into preprocessing errors, you should set ``preprocess_raw_obspack_files`` to ``false`` and supply an already populated directory of manually preprocessed files. Details for how to do this are provided in the `ObsPack diagnostic <https://geos-chem.readthedocs.io/en/stable/gcclassic-user-guide/obspack.html>` documentation for GEOS-Chem; feel free to use the code provided in CHEEREIO as a model. See the :ref:`Observation settings` page for more information on ensemble configuration settings for ObsPack.

Descriptions of the ObsPack preprocessing functions are below.

.. py:function:: make_filter_fxn(start_date,end_date,lat_bounds=None,lon_bounds=None)

   Generate a function that will filter raw ObsPack data (i.e. downloaded directly from NOAA) and keep only data within certain date and location bounds. The output filter function also does additional filtering and reformatting regardless of these bounds.

   :param datetime start_date: Date of earliest ObsPack data to include
   :param datetime end_date: Date of latest ObsPack data to include
   :param list lat_bounds: If filtering by latitude, a list of two latitudes representing minimum and maximum latitude to be kept. If None, ignore. 
   :param list lon_bounds: If filtering by longitude, a list of two longitudes representing minimum and maximum longitudes to be kept. If None, ignore. 
   :return: A filter function for filtering and formatting raw ObsPack data.
   :rtype: Function



.. py:function:: prep_obspack(raw_obspack_dir,gc_obspack_dir,filename_format,start_date,end_date)

   This is a preprocessing function, designed to take raw ObsPack files as downloaded from NOAA and process them into files compatible with CHEEREIO and the GEOS-Chem ObsPack diagnostic.

   :param str raw_obspack_dir: Directory where raw ObsPack data as downloaded from NOAA is stored. CHEEREIO takes this by default from the ``raw_obspack_path`` path in ``ens_config.json.``
   :param str gc_obspack_dir: Directory where processed ObsPack compatible with the GEOS-Chem ObsPack diagnostic will be saved. CHEEREIO takes this by default from the ``gc_obspack_path`` path in ``ens_config.json.``
   :param str filename_format: File format which CHEEREIO will use to save the preprocessed ObsPack data. CHEEREIO takes this by default from the ``obspack_gc_input_file`` entry in ``ens_config.json.``
   :param datetime start_date: Date of earliest ObsPack data to include. CHEEREIO takes this by default from the ``START_DATE`` entry in ``ens_config.json.``
   :param datetime end_date: Date of latest ObsPack data to include. CHEEREIO takes this by default from the ``END_DATE`` entry in ``ens_config.json.``

