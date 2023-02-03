.. _Postprocessing workflow:

The postprocessing workflow  
==========

This section is under construction, check back later!

How to run the default postprocessing workflow
-------------

This section is under construction, check back later!

The postprocessed figures
~~~~~~~~~~~~~

This section is under construction, check back later!

The postprocessed data files
~~~~~~~~~~~~~

This section is under construction, check back later!

The postprocessing API
-------------

This section is under construction, check back later!

Functions in postprocess_tools
~~~~~~~~~~~~~

Here is the documentation for the postprocessing toolkit. Users are especially likely to use :py:func:`combineScaleFactors`. TKTKTK

.. py:function:: globDirs(ensemble_dir,removeNature=False,includeOutputDir=False)

   For a given ensemble directory, get all of the ensemble member run directory paths, their directory names, and their numeric labels and returns them in sorted order.

   :param str ensemble_dir: Path to ensemble run directory. 
   :param bool removeNature: True or False, should we remove ensemble run directory 0 (control run).
   :param bool includeOutputDir: True or False, write paths to go to OutputDir (where the GEOS-Chem model history is stored) or the individual top level ensemble run directory.
   :return: List of (1) list of paths to ensemble run directory members; (2) list of ensemble run directory names; and (3) list of numeric directory labels.
   :rtype: list


.. py:function:: globSubDir(hist_dir,timeperiod=None,hourlysub = 6)

   For a given GEOS-Chem output directory, get all the SpeciesConc files in order for a given time period and return the filenames as a list.

   :param str hist_dir: Path to GEOS-Chem output directory. 
   :param list timeperiod: A list of two datetime objects indicating the start and end of the time period of interest. Leave as None to request the entire time period.
   :param int hourlysub: Only grab files whose hour timestamp is divisible by this number; i.e. 6 means that we grab data every six hours.
   :return: List of filenames for SpeciesConc files.
   :rtype: list

.. py:function:: globSubDirLevelEdge(hist_dir,timeperiod=None,hourlysub = 6)

   As with :py:func:`globSubDir`, but for LevelEdgeDiag files.

   :param str hist_dir: Path to GEOS-Chem output directory. 
   :param list timeperiod: A list of two datetime objects indicating the start and end of the time period of interest. Leave as None to request the entire time period.
   :param int hourlysub: Only grab files whose hour timestamp is divisible by this number; i.e. 6 means that we grab data every six hours.
   :return: List of filenames for LevelEdgeDiag files.
   :rtype: list

.. py:function:: combineScaleFactors(ensemble_dir,output_dir,timeperiod=None,flag_snapshot=False,return_not_write=False)

   Combine emissions scaling factors from across the ensemble and save (or return) them as a single NetCDF or xarray DataSet, with a new dimension called "Ensemble" representing ensemble number. One dataset is saved or returned for each scale factor type. 

   :param str ensemble_dir: Path to CHEEREIO ensemble directory. 
   :param str output_dir: Path to where the combined scaling factor NetCDF should be saved. 
   :param list timeperiod: A list of two datetime objects indicating the start and end of the time period of interest. Leave as None to request the entire time period.
   :param bool flag_snapshot: Flag the output file as a snapshot (True only by the CHEEREIO snapshot script).
   :param bool return_not_write: Return the combined dataset rather than writing it as a NetCDF file.
   :return: If return_not_write is True, a dictionary containing the scale factor names as keys and xarray DataSets with the combined scaling factors as values.
   :rtype: dict

.. py:function:: combineHemcoDiag(ensemble_dir,output_dir,timeperiod=None)

   Combine HEMCO Diagnostics (e.g. emissions) from across the ensemble and save them as a single NetCDF, with a new dimension called "Ensemble" representing ensemble number.

   :param str ensemble_dir: Path to CHEEREIO ensemble directory. 
   :param str output_dir: Path to where the combined HEMCO diagnostic NetCDF should be saved. 
   :param list timeperiod: A list of two datetime objects indicating the start and end of the time period of interest. Leave as None to request the entire time period.

.. py:function:: combineHemcoDiagControl(ensemble_dir,output_dir,timeperiod=None)

   Combine HEMCO Diagnostics (e.g. emissions) from the control run as a single NetCDF.

   :param str ensemble_dir: Path to CHEEREIO ensemble directory. 
   :param str output_dir: Path to where the combined HEMCO diagnostic NetCDF should be saved. 
   :param list timeperiod: A list of two datetime objects indicating the start and end of the time period of interest. Leave as None to request the entire time period.

.. py:function:: makeDatasetForDirectory(hist_dir,species_names,timeperiod=None,hourlysub = 6,subset_rule = 'SURFACE', fullpath_output_name = None)

   Combine GEOS-Chem species concentration output from a single ensemble member as a single dataset and either save to a NetCDF file or return.

   :param str hist_dir: Path to GEOS-Chem output directory. 
   :param list species_names: List of species that we would like to process. 
   :param list timeperiod: A list of two datetime objects indicating the start and end of the time period of interest. Leave as None to request the entire time period.
   :param int hourlysub: Only grab files whose hour timestamp is divisible by this number; i.e. 6 means that we grab data every six hours.
   :param str subset_rule: Which vertical level(s) to save data from. SURFACE is the surface, 850 is the 850hPa level, and ALL is all vertical data.
   :param str fullpath_output_name: Path and filename of the NetCDF file to which we should save the combined data. If ``None``, return the data instead
   :return: If fullpath_output_name is ``None``, an xarray DataSet with the combined concentrations.
   :rtype: DataSet


.. py:function:: makeDatasetForEnsemble(ensemble_dir,species_names,timeperiod=None,hourlysub = 6,subset_rule = 'SURFACE', fullpath_output_name = None)

   Combine GEOS-Chem species concentration output from across the ensemble as a single dataset, with a new Ensemble dimension denoting ensemble member number, and either save to a NetCDF file or return.

   :param str ensemble_dir: Path to CHEEREIO ensemble directory. 
   :param list species_names: List of species that we would like to process. 
   :param list timeperiod: A list of two datetime objects indicating the start and end of the time period of interest. Leave as None to request the entire time period.
   :param int hourlysub: Only grab files whose hour timestamp is divisible by this number; i.e. 6 means that we grab data every six hours.
   :param str subset_rule: Which vertical level(s) to save data from. SURFACE is the surface, 850 is the 850hPa level, and ALL is all vertical data.
   :param str fullpath_output_name: Path and filename of the NetCDF file to which we should save the combined data. If ``None``, return the data instead
   :return: If fullpath_output_name is ``None``, an xarray DataSet with the combined concentrations.
   :rtype: DataSet

.. py:function:: makeYEachAssimPeriod(path_to_bigy_subsets,startdate=None,enddate=None,fullpath_output_name = None)

   Combine the intermediate Y datasets, as output by :ref:`HIST Ensemble`, into a dictionary for the entire ensemble run. These Y datasets contain simulated observations from GEOS-Chem aligned with actual observations, along with other metadata. Returns a combined dictionary, where the keys are timestamps.

   :param str path_to_bigy_subsets: Path to where the intermediate values of Y are saved out by the ensemble (postprocess/bigy in your ensemble installation). 
   :param datetime startdate: Start date for when we should process Y datasets. 
   :param datetime enddate: End date for when we should process Y datasets. 
   :param str fullpath_output_name: Path and filename of the pickle file to which we should save the combined Y dataset. If ``None``, return the data instead.
   :return: If fullpath_output_name is ``None``, an dictionary with the combined Y datasets. The keys are timestamps and the values are the bigY data dictionary for the given day -- entries in these include 'Latitude', 'Longitude', 'Observations', and the species name. See :ref:`HIST Ensemble` for details.
   :rtype: dict


.. py:function:: plotSurfaceCell(ds,species_name,latind,lonind,outfile=None,unit='ppt',includesNature=False)

   Plot a timeseries of the surface concentrations for a single grid cell.

   :param DataSet ds: DataSet of the combined ensemble species concentrations, output by :py:func:`makeDatasetForEnsemble`. 
   :param list species_name: Species to plot. 
   :param int latind: Index of latitude of cell we will plot. 
   :param int lonind: Index of longitude of cell we will plot.
   :param str outfile: Name of image file containing plot. If None, display the plot instead.
   :param str unit: either 'ppm', 'ppb', or 'ppt'; CHEEREIO will multiply the GEOS-Chem mole-mole ratio to get these units.
   :param bool includesNature: True or False, include the no-assimilation control run in plot.


plotSurfaceMean(ds,species_name,outfile=None,unit='ppt',includesNature=False):

tsPlotTotalEmissions(ds_ensemble,ds_prior,collectionName,useLognormal = False, timeslice=None,outfile=None):

tsPlotSatCompare(bigY,species,numens,unit='ppb',observer_name='Observations',useControl=False,outfile=None):

tsPlot(time,ensmean,enssd,species_name,unit,nature=None,priortime=None,prior=None,outfile=None):

makeBigYArrays(bigy,gclat,gclon,nEnsemble,postprocess_save_albedo=False,useControl=False):
return arraysbase

Functions in map_tools
~~~~~~~~~~~~~

This section is under construction, check back later!

plotMap(m,lat,lon,flat,labelname,outfile,clim=None,cmap=None,useLog=False,minval = None):

plotEmissions(m,lat,lon,ppdir, hemco_diags_to_process,plotWithLogScale=True, min_emis=None,min_emis_std=None, plotcontrol=True,useLognormal = False, aggToMonthly=True):

plotScaleFactor(m,lat,lon,ppdir, useLognormal = False, aggToMonthly=True):

agg_to_monthly(dates, to_agg)
return [dates,to_return]

.. _New field in postprocessing:

Adding a new observation field to the postprocessing workflow
-------------

This section is under construction, check back later!

