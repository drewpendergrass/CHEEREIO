.. _Observations:

.. _Existing obs tools:

Existing observation operators
==========

CHEEREIO ships with observation operators that have produced by the community. Users can add their own operators by following the instruction on the :ref:`New observation` page, or they can use the operators listed below. 


.. _TROPOMI tools:

TROPOMI tools
-------------

The TROPOspheric Monitoring Instrument (TROPOMI) onboard Sentinel-5 Precursor satellite satellite measures criteria air pollutants and other trace gases of interest to CHEEREIO users. Currently, there are four TROPOMI operators written for CHEEREIO: CH\ :sub:`4`\, CO (contributed by Sina Voshtani), SO\ :sub:`2`\, and NO\ :sub:`2`\ . Users can follow the pattern in ``tropomi_tools.py`` to add support for additional species.

To activate TROPOMI observations, list "TROPOMI" as an observation type in the ``OBS_TYPE`` setting, as described on the :ref:`Observation settings` page.

TROPOMI fields available to save and plot in CHEEREIO
~~~~~~~~~~~~~

All observation operators save standard data from observations and GEOS-Chem and pass it on to CHEEREIO (including latitude, longitude, observation values, GC simulated observation values, and timestamps). However, individual operators might also be able to save additional data, such as albedo in remote sensing cases or site IDs in surface observation cases. Users can list the additional data they would like to save in the ``EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y`` setting in ``ens_config.json``, following the instructions on the :ref:`postprocessing settings` page. **With the TROPOMI operator, all data read in by the read_tropomi function into the met dictionary are available to be saved and/or plotted** and users can save out additional fields by modifying the ``met`` dictionary. Below is a subset of the supported fields which users can save in TROPOMI:

.. option:: qa_value
   
   Quality assurance value.

.. option:: column_AK
   
   Averaging kernel

.. option:: albedo_swir 
   
   Shortwave infrared albedo (methane only).

.. option:: albedo_nir
   
   Near infrared albedo (methane only).

.. option:: blended_albedo
   
   Blended albedo (methane only).

.. option:: swir_aot
   
   Shortwave infrared aerosol optical thickness (methane only).

.. option:: methane_profile_apriori
   
   A priori methane profile (methane only).


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
   :param int chunk_size: For CO, NO2, and SO2, the number of observations to be processed at once. This is to save memory for TROPOMI observations with high vertical resolution.
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

The remaining two utility functions are very similar. They read TROPOMI level 2 observations from file either from the operational produce format (``read_tropomi``) or from the CHEEREIO-produced compression script (``read_tropomi_compressed``; see :ref:`Compress TROPOMI` for details). Users can select which function they would like to use by specifying the ``WHICH_TROPOMI_**SPC**_PRODUCT`` setting in the ``TROPOMI_**SPC**_extension.json`` file (where ``**SPC**`` is replaced by CH4, CO, NO2, or SO2). ``DEFAULT`` selects the TROPOMI operational product and ``COMPRESSED`` the compressed workflow.

.. py:function:: read_tropomi OR read_tropomi_compressed(filename, species, filterinfo=None, includeObsError = False)

   **Designed for the TROPOMI operational level 2 product.** A utility function which loads TROPOMI observations from file, filters them, and returns a dictionary of important data formatted for input into the ``gcCompare`` method of the ``TROPOMI_Translator`` class. This function is selected when users supply the ``DEFAULT`` value to the ``WHICH_TROPOMI_PRODUCT`` setting in the ``TROPOMI_CH4_extension.json`` file.

   :param str filename: NetCDF file containing TROPOMI observations to be loaded. Expects standard level 2 data. 
   :param str species: Name of species to be loaded. CH4, NO2, SO2, and CO are supported.
   :param dict filterinfo: A dictionary of information about data filtering which is passed to a standard observation operator utility function. See :ref:`Observation filters` for more information
   :param bool includeObsError: True or False, read the errors associated with individual observations. 
   :return: A dictionary containing observation values and metadata, ready for input into the ``gcCompare`` method of the ``TROPOMI_Translator`` class.
   :rtype: dict


.. _Compress TROPOMI:

Utility script to compress TROPOMI files
~~~~~~~~~~~~~

TROPOMI operational files can be very large, and CHEEREIO typically only reads a small fraction of saved fields. Users can optionally run the ``core/compress_tropomi.py`` script from the command line to save out only those fields which CHEEREIO uses. Run the script as follows: 

Within your CHEEREIO conda environment (and probably inside a SLURM batch script, as this process is slow), execute the script (``python compress_tropomi.py``) with the following command line options:

.. option:: --path_to_tropomi_input (-i) [required]
   
   Path to directory containing operational TROPOMI input files which you would like to process. The directory can have a nested structure and can also include other files; CHEEREIO processes all files matching ``S5P_*.nc`` structure. An example entry: ``/hpc/group/shindell/ap851/CO/TROPOMI``

.. option:: --path_to_tropomi_output (-o) [required]
   
   Path to directory where you would like to store your output compressed TROPOMI files. Files saved into the output directory will have an identical nested structure as the input directory and will also have identical filenames. An example entry: ``/hpc/group/shindell/ap851/CO/compressedTROPOMI``

.. option:: --species (-s) [required]
   
   TROPOMI species observed. Supported as of this version: CO, SO2, CH4, and NO2.

.. option:: --delete_original (-d)
   
   Delete original TROPOMI file after compressing? True or False (default: False). Use caution with this option, as CHEEREIO will only save out fields used in the inversion. You may be interested in other fields.

Once the TROPOMI files are compressed, you need to tell CHEEREIO to expect this format. In your inversion run directory, change the ``WHICH_TROPOMI_**SPC**_PRODUCT`` inside the ``extensions/TROPOMI_**SPC**_extension.json`` file to ``COMPRESSED`` (where ``**SPC**`` is replaced by CH4, CO, NO2, or SO2) and ensure that your relevant entry in ``TROPOMI_dirs`` within ``ens_config.json`` points to your compressed output directory (e.g. ``"CO" : "/hpc/group/shindell/ap851/CO/compressedTROPOMI"``). 


A note on TROPOMI NO2 and SO2 autofiltering
~~~~~~~~~~~~~

As of version 1.5.0, certain filters are hard-coded in to the read_tropomi function for NO2 and SO2:

#. Land mask: Only data over land and with no snow or ice are retained (both products)
#. Cloud flag: only keep if below 0.3 (minimal cloud; both products).
#. Solar zenith angle: only keep if below 70 (both products).
#. Volcano flag: only keep if likely not the result of a volcanic eruption (SO2 only)

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

.. _iasi:

IASI tools
-------------

The Infrared Atmospheric Sounding Interferometer (IASI) instrument measures trace gases of interest to CHEEREIO users. Currently, NH\ :sub:`3`\  (as retrieved by the ANNI v4 algorithm from the Université Libre de Bruxelles group) is the only IASI operator written for CHEEREIO, but users can follow the pattern in ``iasi_tools.py`` to add support for additional species.

To activate IASI observations, list "IASI" as an observation type in the ``OBS_TYPE`` setting, as described on the :ref:`Observation settings` page.

.. warning::
   CHEEREIO saves out IASI data in μmol per mol (both from GEOS-Chem and adjusted from the operational product; multiplying each by :math:`10^6`). This is to avoid some floating point errors.


IASI fields available to save and plot in CHEEREIO
~~~~~~~~~~~~~

Users can list additional fields (beyond the minimum observations and spatiotemporal location) they would like to save from IASI in the ``EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y`` setting in ``ens_config.json``, following the instructions on the :ref:`postprocessing settings` page. **With the IASI operator, all data read in by the read_iasi function are available to be saved and/or plotted.**

Note that the IASI operator calculates whether an observation should be discarded *after* the operator is applied, as we can only evaluate whether the observation is in appropriate error bounds recommended by ULB once we replace the prior column with the GC column. This is implemented through the postfilter functionality, as described on the :ref:`postfilter` entry.

IASI operator support functions
~~~~~~~~~~~~~

The IASI observation operator calls two utility functions to process observations from file and format them in such a way that the ``gcCompare`` method of the ``IASI_Translator`` class can perform relevant computations. They are documented below:

.. py:function:: read_iasi(filename, species, filterinfo=None, includeObsError = False)

   A utility function which loads IASI observations from file, filters them, and returns a dictionary of important data formatted for input into the ``gcCompare`` method of the ``IASI_Translator`` class.

   :param str filename: NetCDF file containing IASI observations to be loaded. Expects standard level 2 data. 
   :param str species: Name of species to be loaded. Currently, only "NH3" is supported.
   :param dict filterinfo: A dictionary of information about data filtering which is passed to a standard observation operator utility function. See :ref:`Observation filters` for more information
   :param bool includeObsError: True or False, read the errors associated with individual observations. 
   :return: A dictionary containing observation values and metadata, ready for input into the ``gcCompare`` method of the ``IASI_Translator`` class.
   :rtype: dict

.. py:function:: GC_to_sat_levels(GC_SPC, GC_bxheight, sat_edges)

   See the function from TROPOMI_tools. Note that for the IASI product we only have height from surface, so the function uses GC Boxheight diagnostic to do a pressure level regrid approximation.

.. _iasi error:

IASI observational error
~~~~~~~~~~~~~

The IASI operational product separates observational error into systematic and random components (see `the IASI v4 paper <https://doi.org/10.5194/amt-16-5009-2023>`__ for details). If you would like to use this information in your inversion (recommended) when aggregating to the GEOS-Chem grid (super-observations) then be sure to do the following in ``ens_config.json``:

#. Set the NH3 IASI observation key under ``OBS_ERROR_TYPE`` to ``product``
#. Set the NH3 IASI observation key under ``AV_TO_GC_GRID`` to ``True``
#. Set the NH3 IASI observation key under ``SUPER_OBSERVATION_FUNCTION`` to ``sys_rand_err``
#. Provide an estimate for systematic transport error under ``OTHER_OBS_ERROR_PARAMETERS`` for a NH3 IASI observation subdictionary by including the following key-value pair: ``"transport_error":"TKTKT"`` where you would replace TKTKT with your estimate. See CHEEREIO papers for estimating this value.

The transport and systematic error components are taken to be fully correlated and thus irreducible from averaging, while the random component follows the central limit theorem and reduces as the observation count increases; the total super-observation error covariance is obtained by adding these components in quadrature. See the :ref:`New superobservation` entry for more information

.. _tccon:

TCCON tools
-------------

The Total Carbon Column Observing Network (TCCON) measures trace gases of interest to CHEEREIO users. Currently, CO and nitrous oxide are supported for CHEEREIO, but users can follow the pattern in ``tccon_tools.py`` to add support for additional species. To activate TCCON observations, list "TCCON" as an observation type in the ``OBS_TYPE`` setting, as described on the :ref:`Observation settings` page.

The TCCON operator was originally built by Sina Voshtani. See :ref:`appropriate citations` for a recommended citation.

TCCON fields available to save and plot in CHEEREIO
~~~~~~~~~~~~~

Users can list additional fields (beyond the minimum observations and spatiotemporal location) they would like to save from TCCON in the ``EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y`` setting in ``ens_config.json``, following the instructions on the :ref:`postprocessing settings` page. **With the TCCON operator, all data read in by the read_tccon function are available to be saved and/or plotted.**

TCCON operator support functions
~~~~~~~~~~~~~

The TCCON observation operator calls two utility functions to process observations from file and format them in such a way that the ``gcCompare`` method of the ``TCCON_Translator`` class can perform relevant computations. They are documented below:

.. py:function:: read_tccon(filename, species, filterinfo=None, includeObsError = False,doN2OCorrectionPT700=False)

   A utility function which loads TCCON observations from file, filters them, and returns a dictionary of important data formatted for input into the ``gcCompare`` method of the ``TCCON_Translator`` class.

   :param str filename: NetCDF file containing IASI observations to be loaded. Expects data produced through the prep_tccon_aggregated.py script, as described below. 
   :param str species: Name of species to be loaded. Currently, only "N2O" and "CO" are supported.
   :param dict filterinfo: A dictionary of information about data filtering which is passed to a standard observation operator utility function. See :ref:`Observation filters` for more information
   :param bool includeObsError: True or False, read the errors associated with individual observations.
   :param bool doN2OCorrectionPT700: True or False, do the temperature correction for N2O in the GGG2020.0 product set? This won't be necessary after 2020.1 is released.
   :return: A dictionary containing observation values and metadata, ready for input into the ``gcCompare`` method of the ``TCCON_Translator`` class.
   :rtype: dict

.. py:function:: correct_xn2o_from_pt700(xn2o,prior_temperature,prior_pressure,xn2o_error=None,n2o_aicf=0.9821, m=0.000626, b=0.787)

   Handle temperature-dependent bias in TCCON N2O. Based on Josh Laughner's code for GGG2020: py_tccon_netcdf/write_tccon_netcdf/bias_corrections.py. Will no longer be needed after GGG2020.1 is released.

.. py:function:: _compute_pt700(prior_temperature,prior_pressure)

   Handle temperature-dependent bias in TCCON N2O. Based on Josh Laughner's code for GGG2020: py_tccon_netcdf/write_tccon_netcdf/bias_corrections.py. Will no longer be needed after GGG2020.1 is released.

.. py:function:: GC_to_sat_levels(GC_SPC, GC_edges, sat_edges)

   See the function from TROPOMI_tools.

.. py:function:: gravity(altitudes, latitudes)

   Compute g at vertical layers for each observation site.

.. py:function:: integrate_column(gas_profile,h2o_profile,obh2o_profile,obpout,obpressure_profile,altitude_profile,ensemble_profile,oblat,AK)

   This is the main function for TCCON column intergation, using its a-priori and averaging kernels.

   :param array gas_profile: The model gas profile of interest. 
   :param array h2o_profile: The model h2o profile.
   :param array obh2o_profile: The observation h2o profile. 
   :param array pressure_profile: The model pressure profile that corresponds with the gas profile in hPa.
   :param array ensemble_profile: The ensemble profile from equation 25 in Rodgers and Connor 2000 - this will likely be the a priori profile from GFIT; most often multiplied by the scaling factors (VSF) from GFIT for the spectra near the aircraft overpass
   :param float obalt: The altitude of the ground-based site in m - geometric altitude
   :param float oblat: The latitude (in degrees) of the ground-based site
   :param array AK: The averaging kernels for all windows of the molecule of interest, in a structure


TCCON preprocessing for CHEEREIO
~~~~~~~~~~~~~

CHEEREIO has built in preprocessing functions to translate raw GGG2020 TCCON data, as downloaded from the CalTech site, into a form compatible with the CHEEREIO TCCON operator. To use this, execute from the command line (in CHEEREIO/core) before installing CHEEREIO: ``python prep_tccon_aggregated.py ARGUMENTS``

Arguments are as follows:

#. **-i or --input_path (required)**: Path to your input TCCON files, e.g. downloaded from tccondata.org.
#. **-o or --output_path (required)**: Path to where to save your pre-processed TCCON files, ready for CHEEREIO.
#. **-s or --start_time (required)**: Start time, in format YYYY-MM-DDTHH:MM:SS. For example, 2023-01-01T00:00:00. Will only process data from after this time.
#. **-e or --end_time (required)**: End time, in format YYYY-MM-DDTHH:MM:SS. For example, 2023-02-01T00:00:00. Will only process data from before this time.
#. **-s2k or --species_to_keep (optional)**: Species to keep in your pre-processed TCCON files (e.g. just "co" or "co,co2"). If multiple, comma separated.
#. **-p or --input_file_pattern (optional)**: File pattern for input tccon files. If you have the public files, you do not need to modify (default is ``*public.qc.nc``).

.. _ObsPack tools:

ObsPack tools
-------------

`ObsPack <https://doi.org/10.5194/essd-6-375-2014>`__ is a standardized dataset containing measurements from  surface monitors distributed around the world, aimed at carbon cycle studies. CHEEREIO users often use ObsPack data for validation of CO or CH4 inversions, or as observations for the inversion itself. The ObsPack observation operator, contained in the ``obspack_tools.py`` file, wraps around the Obspack `diagnostic <https://geos-chem.readthedocs.io/en/stable/gcclassic-user-guide/obspack.html>`__ produced by GEOS-Chem and translates it into a form acceptable to CHEEREIO.

To activate ObsPack, see the :ref:`Observation settings` page for information on the correct settings for  ``ens_config.json``. 

ObsPack fields available to save and plot in CHEEREIO
~~~~~~~~~~~~~

All observation operators save standard data from observations and GEOS-Chem and pass it on to CHEEREIO (including latitude, longitude, observation values, GC simulated observation values, and timestamps). However, individual operators might also be able to save additional data, such as albedo in remote sensing cases or site IDs in surface observation cases. Users can list the additional data they would like to save in the ``EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y`` setting in ``ens_config.json``, following the instructions on the :ref:`postprocessing settings` page. Below are a list of the supported fields which users can save in ObsPack:

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

CHEEREIO has built in preprocessing functions to translate raw ObsPack data, as downloaded from NOAA, into a form compatible with the GEOS-Chem ObsPack diagnostic. To use this functionality, set ``preprocess_raw_obspack_files`` to ``true`` in ``ens_config.json`` and provide a path to the raw files in the ``raw_obspack_path`` entry. However, some users report that NOAA ObsPack data is not quite standardized. If you run into preprocessing errors, you should set ``preprocess_raw_obspack_files`` to ``false`` and supply an already populated directory of manually preprocessed files. Details for how to do this are provided in the `ObsPack diagnostic <https://geos-chem.readthedocs.io/en/stable/gcclassic-user-guide/obspack.html>`__ documentation for GEOS-Chem; feel free to use the code provided in CHEEREIO as a model. See the :ref:`Observation settings` page for more information on ensemble configuration settings for ObsPack.

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

