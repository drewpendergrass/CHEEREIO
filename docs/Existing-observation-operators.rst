.. _Observations:

.. _Existing obs tools:

Existing observation operators
==========

CHEEREIO ships with observation operators that have produced by the community. Users can add their own operators by following the instruction on the :ref:`New observation` page, or they can use the operators listed below. 


.. _TROPOMI tools:

TROPOMI tools
-------------

The TROPOspheric Monitoring Instrument (TROPOMI) onboard Sentinel-5 Precursor satellite satellite measures criteria air pollutants and other trace gases of interest to CHEEREIO users. Currently, there are two TROPOMI operators written for CHEEREIO: CH\ :sub:`4`\ and CO (the latter contributed by Sina Voshtani). Users can follow the pattern in ``tropomi_tools.py`` to add support for additional species. An NO\ :sub:`2`\ operator is partially written but not yet functional as of version 1.2.

To activate TROPOMI observations, list "TROPOMI" as an observation type in the ``OBS_TYPE`` setting, as described on the :ref:`Observation settings` page.

TROPOMI fields available to save and plot in CHEEREIO
~~~~~~~~~~~~~

All observation operators save standard data from observations and GEOS-Chem and pass it on to CHEEREIO (including latitude, longitude, observation values, GC simulated observation values, and timestamps). However, individual operators might also be able to save additional data, such as albedo in remote sensing cases or site IDs in surface observation cases. Users can list the additional data they would like to save in the ``EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y`` setting in ``ens_config.json``, following the instructions on the :ref:`postprocessing settings` page. **With the TROPOMI operator, all data read in by the read_tropomi function are available to be saved and/or plotted.** Below is a subset of the supported fields which users can save in TROPOMI:

.. option:: qa_value
   
   Quality assurance value.

.. option:: column_AK
   
   Satellite averaging kernel.

.. option:: albedo_swir
   
   Shortwave infrared albedo.

.. option:: albedo_nir
   
   Near infrared albedo.

.. option:: blended_albedo
   
   Blended albedo.

.. option:: methane_profile_apriori
   
   A priori methane profile. (Methane only)


TROPOMI operator support functions
~~~~~~~~~~~~~

The TROPOMI observation operator calls a handful of utility functions to process observations from file and format them in such a way that the ``gcCompare`` method of the ``TROPOMI_Translator`` class can perform relevant computations. 

The first two utility functions are designed to remap GEOS-Chem pressure levels to satellite pressure levels and apply the averaging kernel.

.. py:function:: GC_to_sat_levels(GC_SPC, GC_edges, sat_edges, species, chunk_size=10000)

   Takes as input GEOS-Chem data and pressure level edges, as well as satellite pressure levels, and calculate GEOS-Chem data values on satellite pressure levels

   :param array GC_SPC: A NumPy array containing GEOS-Chem columns. 
   :param array GC_edges: A NumPy array containing GEOS-Chem pressure level edges.
   :param array sat_edges: A NumPy array containing TROPOMI pressure level edges
   :param str species: Species to be processed.
   :param int chunk_size: For CO, the number of observations to be processed at once. This is to save memory for TROPOMI observations with high vertical resolution.
   :return: A NumPy array containing GEOS-Chem columns remapped to be on the TROPOMI pressure levels.
   :rtype: array

.. py:function:: apply_avker(sat_avker, sat_pressure_weight, GC_SPC, sat_prior=None,filt=None)

   Apply the averaging kernel

   :param array sat_avker: TROPOMI averaging kernel. 
   :param array sat_prior: The satellite prior profile in ppb, optional (used for CH4).
   :param array sat_pressure_weight: The relative pressure weights for each level
   :param array GC_SPC: The GC species on the satellite levels, output by GC_to_sat_levels
   :param array filt: A filter, optional
   :return: A NumPy array containing simulated GEOS-Chem values such that they are directly comparable to TROPOMI (i.e. because the averaging kernel has been applied).
   :rtype: array

The remaining three utility functions are very similar. They read TROPOMI level 2 observations from file, but for methane there are a variety of TROPOMI observations with a variety of formattings (operational, science product, Harvard-specific standard). Each function is designed to read a different formatting. Users can select which function they would like to use by specifying the ``WHICH_TROPOMI_PRODUCT`` setting in the ``TROPOMI_CH4_extension.json`` file. ``DEFAULT`` selects the TROPOMI operational product, ``ACMG`` for the ACMG/Harvard TROPOMI product, and ``BLENDED`` for Belasus et al., 2023 which also works for the TROPOMI science product. For CO, users must select the ``DEFAULT`` option.

.. py:function:: read_tropomi(filename, species, filterinfo=None, includeObsError = False)

   **Designed for the TROPOMI operational level 2 product.** A utility function which loads TROPOMI observations from file, filters them, and returns a dictionary of important data formatted for input into the ``gcCompare`` method of the ``TROPOMI_Translator`` class. This function is selected when users supply the ``DEFAULT`` value to the ``WHICH_TROPOMI_PRODUCT`` setting in the ``TROPOMI_CH4_extension.json`` file.

   :param str filename: NetCDF file containing TROPOMI observations to be loaded. Expects standard level 2 data. 
   :param str species: Name of species to be loaded. Currently, only "CH4" is supported, though an "NO2" operator is partially written.
   :param dict filterinfo: A dictionary of information about data filtering which is passed to a standard observation operator utility function. See :ref:`Observation filters` for more information
   :param bool includeObsError: True or False, read the errors associated with individual observations. 
   :return: A dictionary containing observation values and metadata, ready for input into the ``gcCompare`` method of the ``TROPOMI_Translator`` class.
   :rtype: dict

   .. py:function:: read_tropomi_acmg(filename, species, filterinfo=None, includeObsError = False)

   **Designed for the Harvard/ACMG version of the TROPOMI operational level 2 product.** A utility function which loads TROPOMI observations from file, filters them, and returns a dictionary of important data formatted for input into the ``gcCompare`` method of the ``TROPOMI_Translator`` class. This function is selected when users supply the ``ACMG`` value to the ``WHICH_TROPOMI_PRODUCT`` setting in the ``TROPOMI_CH4_extension.json`` file.

   :param str filename: NetCDF file containing TROPOMI observations to be loaded. Expects standard level 2 data. 
   :param str species: Name of species to be loaded. Currently, only "CH4" is supported, though an "NO2" operator is partially written.
   :param dict filterinfo: A dictionary of information about data filtering which is passed to a standard observation operator utility function. See :ref:`Observation filters` for more information
   :param bool includeObsError: True or False, read the errors associated with individual observations. 
   :return: A dictionary containing observation values and metadata, ready for input into the ``gcCompare`` method of the ``TROPOMI_Translator`` class.
   :rtype: dict

   .. py:function:: read_tropomi_gosat_corrected(filename, species, filterinfo=None, includeObsError = False)

   **Designed for the Belasus et al., 2023 version of the TROPOMI operational level 2 product,** but also works for the TROPOMI science product. A utility function which loads TROPOMI observations from file, filters them, and returns a dictionary of important data formatted for input into the ``gcCompare`` method of the ``TROPOMI_Translator`` class. This function is selected when users supply the ``BLENDED`` value to the ``WHICH_TROPOMI_PRODUCT`` setting in the ``TROPOMI_CH4_extension.json`` file.

   :param str filename: NetCDF file containing TROPOMI observations to be loaded. Expects standard level 2 data. 
   :param str species: Name of species to be loaded. Currently, only "CH4" is supported, though an "NO2" operator is partially written.
   :param dict filterinfo: A dictionary of information about data filtering which is passed to a standard observation operator utility function. See :ref:`Observation filters` for more information
   :param bool includeObsError: True or False, read the errors associated with individual observations. 
   :return: A dictionary containing observation values and metadata, ready for input into the ``gcCompare`` method of the ``TROPOMI_Translator`` class.
   :rtype: dict

.. _OMI tools:

OMI tools
-------------

NASA's Ozone Monitoring Instrument (OMI) onboard the Aura satellite measures criteria air pollutants and other trace gases of interest to CHEEREIO users. Currently, NO\ :sub:`2`\ is the only OMI operator written for CHEEREIO, but users can follow the pattern in ``omi_tools.py`` to add support for additional species.

To activate OMI observations, list "OMI" as an observation type in the ``OBS_TYPE`` setting, as described on the :ref:`Observation settings` page.

OMI fields available to save and plot in CHEEREIO
~~~~~~~~~~~~~

All observation operators save standard data from observations and GEOS-Chem and pass it on to CHEEREIO (including latitude, longitude, observation values, GC simulated observation values, and timestamps). However, individual operators might also be able to save additional data, such as albedo in remote sensing cases or site IDs in surface observation cases. Users can list the additional data they would like to save in the ``EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y`` setting in ``ens_config.json``, following the instructions on the :ref:`postprocessing settings` page. 

Currently, **no additional ObsData fields are available to be saved and/or plotted with the OMI operator** as written. Additional fields can be added by saving more metadata into the ObsData object with the ``addData`` function within the ``gcCompare`` method of the ``OMI_Translator`` class. See the :ref:`Observations` page for more details.

OMI operator support functions
~~~~~~~~~~~~~

The OMI observation operator calls two utility functions to process observations from file and format them in such a way that the ``gcCompare`` method of the ``OMI_Translator`` class can perform relevant computations. They are documented below:

.. py:function:: read_omi(filename, species, filterinfo=None, includeObsError = False)

   A utility function which loads OMI observations from file, filters them, and returns a dictionary of important data formatted for input into the ``gcCompare`` method of the ``OMI_Translator`` class.

   :param str filename: NetCDF file containing OMI observations to be loaded. Expects standard level 2 data. 
   :param str species: Name of species to be loaded. Currently, only "NO2" is supported.
   :param dict filterinfo: A dictionary of information about data filtering which is passed to a standard observation operator utility function. See :ref:`Observation filters` for more information
   :param bool includeObsError: True or False, read the errors associated with individual observations. 
   :return: A dictionary containing observation values and metadata, ready for input into the ``gcCompare`` method of the ``OMI_Translator`` class.
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

