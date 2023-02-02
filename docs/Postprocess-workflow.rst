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

This section is under construction, check back later!

.. py:function:: globDirs(ensemble_dir,removeNature=False,includeOutputDir=False)

   For a given ensemble directory, get all of the ensemble member run directory paths, their directory names, and their numeric labels and returns them in sorted order.

   :param str ensemble_dir: Path to ensemble run directory. 
   :param bool removeNature: True or False, should we remove ensemble run directory 0 (control run).
   :param bool includeOutputDir: True or False, write paths to go to OutputDir (where the GEOS-Chem model history is stored) or the individual top level ensemble run directory.
   :return: List of (1) list of paths to ensemble run directory members; (2) list of ensemble run directory names; and (3) list of numeric directory labels.
   :rtype: list

globSubDir(hist_dir,timeperiod=None,hourlysub = 6)
return specconc_list

globSubDirLevelEdge(hist_dir,timeperiod=None,hourlysub = 6)
return edgeconc_list

combineScaleFactors(ensemble_dir,output_dir,timeperiod=None,flag_snapshot=False,return_not_write=False)
to_return

combineHemcoDiag(ensemble_dir,output_dir,timeperiod=None)

combineHemcoDiagControl(control_dir,output_dir,timeperiod=None)

makeDatasetForDirectory(hist_dir,species_names,timeperiod=None,hourlysub = 6,subset_rule = 'SURFACE', fullpath_output_name = None)
return ds

makeDatasetForEnsemble(ensemble_dir,species_names,timeperiod=None,hourlysub = 6,subset_rule = 'SURFACE',fullpath_output_name = None)
return ds

makeYEachAssimPeriod(path_to_bigy_subsets,startdate=None,enddate=None,fullpath_output_name = None):
return masterY

plotSurfaceCell(ds,species_name,latind,lonind,outfile=None,unit='ppt',includesNature=False):

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

