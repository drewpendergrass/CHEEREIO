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

hemco_diags_to_process = data['hemco_diags_to_process']
pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"
ens_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs"
control_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/control_run"
savelevel = data['SaveLevelEdgeDiags']
controlvec = data['CONTROL_VECTOR_CONC']
postprocess_save_albedo = data['postprocess_save_albedo']=="True"
nEnsemble = int(data['nEnsemble'])
statevec = data['STATE_VECTOR_CONC']
emisvec = list(data['CONTROL_VECTOR_EMIS'].keys())
START_DATE=datetime.strptime(data['START_DATE'], "%Y%m%d")
ASSIM_START_DATE=datetime.strptime(data['ASSIM_START_DATE'], "%Y%m%d")
endtime=datetime.strptime(data['END_DATE'], "%Y%m%d")
ASSIM_TIME=data['ASSIM_TIME']
delta = timedelta(hours=int(ASSIM_TIME))
starttime = START_DATE
timeperiod = (starttime,endtime)
avtogcgrid = data['AV_TO_GC_GRID']=="True"
useLevelEdge=data['SaveLevelEdgeDiags']=="True"
useStateMet=data['SaveStateMet']=="True"
useArea=data['SaveArea']=="True"
useControl=data['DO_CONTROL_RUN']=="true"

print('Starting scale factor postprocessing.')
if len(emisvec) > 0:
	if not exists(f'{pp_dir}/{emisvec[0]}_SCALEFACTOR.nc'):
		pt.combineScaleFactors(ens_dir,pp_dir)
		print('Scale factor postprocessing file not detected; generating now.')
	else:
		print('Detected existing scale factor postprocessing file. Loading now.')
	scalefactor_files = glob(f'{pp_dir}/*_SCALEFACTOR.nc')
	for scalefactor in scalefactor_files:
		sf_name = '_'.join(scalefactor.split('/')[-1].split('_')[0:-1])

print('Scale factor postprocessing complete.')
print('Starting HEMCO diagnostic (e.g. emissions) postprocessing.')

try:
	hemcodiag = xr.open_dataset(f'{pp_dir}/combined_HEMCO_diagnostics.nc')
except FileNotFoundError:
	print('HEMCO diagnostic postprocessing file not detected; generating now.')
	pt.combineHemcoDiag(ens_dir,pp_dir)
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
	daterange = np.arange(ASSIM_START_DATE,endtime+timedelta(hours=1),delta).astype(datetime)
	dates_string_array = [dateval.strftime("%Y%m%d_%H%M") for dateval in daterange]
	try:
		with open(f"{pp_dir}/bigY.pkl",'rb') as f:
			bigy=pickle.load(f)
	except FileNotFoundError:
		print('Observation postprocessing files not detected; generating now.')
		bigy = pt.makeYEachAssimPeriod(dates_string_array,useLevelEdge=useLevelEdge,useStateMet = useStateMet,useArea=useArea,use_numav=avtogcgrid,use_albedo=postprocess_save_albedo,useControl = useControl, fullpath_output_name=f"{pp_dir}/bigY.pkl")
	print('Simulated observation vs actual observation (Y) postprocessed and loaded.')

if useControl:
	print('Starting control run HEMCO diagnostic (e.g. emissions) postprocessing.')
	try:
		hemcocontroldiag = xr.open_dataset(f'{pp_dir}/control_HEMCO_diagnostics.nc')
	except FileNotFoundError:
		print('Control run HEMCO diagnostic postprocessing file not detected; generating now.')
		pt.combineHemcoDiagControl(control_dir,pp_dir)
		hemcocontroldiag = xr.open_dataset(f'{pp_dir}/control_HEMCO_diagnostics.nc')
	for collection in hemco_diags_to_process:
		pt.tsPlotTotalEmissions(ds_ensemble=hemcodiag,ds_prior=hemcocontroldiag,collectionName=collection,timeslice=[START_DATE,endtime], outfile=f'{pp_dir}/timeseries_totalemissions_{collection}_against_prior.png')
	print('Control run HEMCO diagnostic (e.g. emissions) postprocessed and loaded.')

if "calc850" in sys.argv:
	print('Calculating/loading 850hPa pressure level')
	try:
		ds850 = xr.open_dataset(f'{pp_dir}/controlvar_pp_850hPa.nc')
	except FileNotFoundError:
		_ = pt.makeDatasetForEnsemble(ens_dir,controlvec,timeperiod,subset_rule="850",fullpath_output_name=f'{pp_dir}/controlvar_pp_850hPa.nc')
		ds850 = xr.open_dataset(f'{pp_dir}/controlvar_pp_850hPa.nc')


for spec in controlvec:
	if "histprocess" in sys.argv:
		pt.tsPlotSatCompare(bigy,spec,nEnsemble,unit='ppb',satellite_name='TROPOMI',outfile=f'{pp_dir}/satellite_ts_compare_{spec}.png')
	pt.plotSurfaceMean(ds,spec,outfile=f'{pp_dir}/surfmean_ts_{spec}.png',includesNature=False)
	if "calc850" in sys.argv:
		 pt.plotSurfaceMean(ds850,spec,outfile=f'{pp_dir}/mean850hPa_ts_{spec}.png',includesNature=False)

