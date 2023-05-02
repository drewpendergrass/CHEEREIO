import xarray as xr
import postprocess_tools as pt
from glob import glob
from datetime import datetime,timedelta
import pickle
import pandas as pd
import numpy as np
from os.path import exists
import sys 
sys.path.append('../core')
import settings_interface as si 

data = si.getSpeciesConfig()

gclat,gclon = si.getLatLonVals(data)
gclat = np.array(gclat)
gclon = np.array(gclon)

hemco_diags_to_process = data['hemco_diags_to_process']
pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"
ens_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs"
useControl=data['DO_CONTROL_RUN']=="true"
controlInEns = data['DO_CONTROL_WITHIN_ENSEMBLE_RUNS']=="true"

if useControl and controlInEns:
	control_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs/{data['RUN_NAME']}_0000"
elif useControl and not controlInEns:
	control_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/control_run"

savelevel = data['SaveLevelEdgeDiags']
controlvec = data['CONTROL_VECTOR_CONC']
observed_species = data['OBSERVED_SPECIES']
obs_units = data['OBSERVATION_UNITS']

nEnsemble = int(data['nEnsemble'])
statevec = data['STATE_VECTOR_CONC']
emisvec = list(data['CONTROL_VECTOR_EMIS'].keys())
lognormalErrors=data['lognormalErrors']=="true"
POSTPROCESS_START_DATE=datetime.strptime(data['POSTPROCESS_START_DATE'], "%Y%m%d")
POSTPROCESS_END_DATE=datetime.strptime(data['POSTPROCESS_END_DATE'], "%Y%m%d")
ASSIM_TIME=data['ASSIM_TIME']
delta = timedelta(hours=int(ASSIM_TIME))
timeperiod = (POSTPROCESS_START_DATE,POSTPROCESS_END_DATE)
av_to_gc_grid = data['AV_TO_GC_GRID']
useLevelEdge=data['SaveLevelEdgeDiags']=="True"
useStateMet=data['SaveStateMet']=="True"
useArea=data['SaveArea']=="True"

observers_to_plot_as_points=data['OBSERVERS_TO_PLOT_AS_POINTS']
extra_obsdata_fields=data['EXTRA_OBSDATA_FIELDS_TO_REGRID_AND_PLOT']

print('Starting scale factor postprocessing.')
if len(emisvec) > 0:
	if not exists(f'{pp_dir}/{emisvec[0]}_SCALEFACTOR.nc'):
		print('Scale factor postprocessing file not detected; generating now.')
		pt.combineScaleFactors(ens_dir,pp_dir,timeperiod)
	else:
		print('Detected existing scale factor postprocessing file.')

print('Scale factor postprocessing complete.')
print('Starting HEMCO diagnostic (e.g. emissions) postprocessing.')

try:
	area = np.load(f'{pp_dir}/area.npy')
except OSError:
	print('Area postprocessing file not detected; generating now.')
	pt.getArea(ens_dir,pp_dir)
	area = np.load(f'{pp_dir}/area.npy')

try:
	hemcodiag = xr.open_dataset(f'{pp_dir}/combined_HEMCO_diagnostics.nc')
except FileNotFoundError:
	print('HEMCO diagnostic postprocessing file not detected; generating now.')
	pt.combineHemcoDiag(ens_dir,pp_dir,timeperiod)
	hemcodiag = xr.open_dataset(f'{pp_dir}/combined_HEMCO_diagnostics.nc')

print('HEMCO diagnostic (e.g. emissions) postprocessed and loaded.')
print('Starting concentration (surface) postprocessing.')

try:
	ds = xr.open_dataset(f'{pp_dir}/controlvar_pp.nc')
except FileNotFoundError:
	print('Surface concentration postprocessing file not detected; generating now.')
	_ = pt.makeDatasetForEnsemble(ens_dir,controlvec,timeperiod,fullpath_output_name=f'{pp_dir}/controlvar_pp.nc')
	ds = xr.open_dataset(f'{pp_dir}/controlvar_pp.nc')	

print('Surface concentration data postprocessed and loaded.')


if "histprocess" in sys.argv:
	print('Starting simulated observation vs actual observation (Y) postprocessing.')
	try:
		with open(f"{pp_dir}/bigY.pkl",'rb') as f:
			bigy=pickle.load(f)
	except FileNotFoundError:
		print('Observation postprocessing files not detected; generating now.')
		bigy = pt.makeYEachAssimPeriod(path_to_bigy_subsets=f"{pp_dir}/bigy",assim_time=int(ASSIM_TIME),startdate=POSTPROCESS_START_DATE,enddate=POSTPROCESS_END_DATE, fullpath_output_name=f"{pp_dir}/bigY.pkl")
	print('Simulated observation vs actual observation (Y) postprocessed and loaded.')

	print('Beginning to regrid simulated observation vs actual observation (Y) for plotting and analysis.')
	try: 		
		with open(f"{pp_dir}/bigy_arrays_for_plotting.pkl",'rb') as f:
			arraysbase=pickle.load(f)
		except FileNotFoundError:
			print('Gridded observation postprocessing (Y) files not detected; generating now.')
			arraysbase = pt.makeBigYArrays(bigy,gclat,gclon,nEnsemble,av_to_grid=av_to_gc_grid, observers_to_plot_as_points=observers_to_plot_as_points,extra_obsdata_fields=extra_obsdata_fields,useControl=useControl)
			file = open(f'{pp_dir}/bigy_arrays_for_plotting.pkl',"wb")
			pickle.dump(arraysbase,file)
			file.close()
		print('Observation postprocessing (Y) files regridded and loaded.')


if useControl:
	print('Starting control run HEMCO diagnostic (e.g. emissions) postprocessing.')
	try:
		hemcocontroldiag = xr.open_dataset(f'{pp_dir}/control_HEMCO_diagnostics.nc')
	except FileNotFoundError:
		print('Control run HEMCO diagnostic postprocessing file not detected; generating now.')
		pt.combineHemcoDiagControl(control_dir,pp_dir,timeperiod)
		hemcocontroldiag = xr.open_dataset(f'{pp_dir}/control_HEMCO_diagnostics.nc')
	for collection in hemco_diags_to_process:
		pt.tsPlotTotalEmissions(ds_ensemble=hemcodiag,ds_prior=hemcocontroldiag,area=area,collectionName=collection,useLognormal = lognormalErrors,timeslice=[POSTPROCESS_START_DATE,POSTPROCESS_END_DATE], outfile=f'{pp_dir}/timeseries_totalemissions_{collection}_against_prior.png')
	print('Control run HEMCO diagnostic (e.g. emissions) postprocessed and loaded.')

if "calc850" in sys.argv:
	print('Calculating/loading 850hPa pressure level')
	try:
		ds850 = xr.open_dataset(f'{pp_dir}/controlvar_pp_850hPa.nc')
	except FileNotFoundError:
		_ = pt.makeDatasetForEnsemble(ens_dir,controlvec,timeperiod,subset_rule="850",fullpath_output_name=f'{pp_dir}/controlvar_pp_850hPa.nc')
		ds850 = xr.open_dataset(f'{pp_dir}/controlvar_pp_850hPa.nc')


for spec in controlvec:
	pt.plotSurfaceMean(ds,spec,outfile=f'{pp_dir}/surfmean_ts_{spec}.png',includesNature=False)
	if "calc850" in sys.argv:
		 pt.plotSurfaceMean(ds850,spec,outfile=f'{pp_dir}/mean850hPa_ts_{spec}.png',includesNature=False)

if "histprocess" in sys.argv:
	for spec in observed_species:
		pt.tsPlotSatCompare(bigy,spec,nEnsemble,unit=obs_units[spec],observer_name=data['OBS_TYPE'][spec],useControl=False,outfile=f'{pp_dir}/observations_ts_compare_{spec}.png')
		if useControl:
			pt.tsPlotSatCompare(bigy,spec,nEnsemble,unit=obs_units[spec],observer_name=data['OBS_TYPE'][spec],useControl=True,outfile=f'{pp_dir}/observations_ts_compare_{spec}_w_control.png')
