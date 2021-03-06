import json
import xarray as xr
import postprocess_tools as pt
from glob import glob
from datetime import datetime,timedelta
import sys 
import pickle
import pandas as pd
import numpy as np

with open('../ens_config.json') as f:
	data = json.load(f)

pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"
ens_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs"
savelevel = data['SaveLevelEdgeDiags']
controlvec = data['CONTROL_VECTOR_CONC']
postprocess_save_albedo = data['postprocess_save_albedo']=="True"
nEnsemble = int(data['nEnsemble'])
statevec = data['STATE_VECTOR_CONC']
emisvec = data['CONTROL_VECTOR_EMIS']
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

pt.combineScaleFactors(ens_dir,pp_dir)
scalefactor_files = glob(f'{pp_dir}/*_SCALEFACTOR.nc')
for scalefactor in scalefactor_files:
	sf_name = '_'.join(scalefactor.split('/')[-1].split('_')[0:-1])
	pt.plotEmissionsCell(scalefactor,30,59,outfile=f'{pp_dir}/wuhan_cell_emis_{sf_name}.png')

try:
	ds = xr.open_dataset(f'{pp_dir}/controlvar_pp.nc')
except FileNotFoundError:
	_ = pt.makeDatasetForEnsemble(ens_dir,controlvec,timeperiod,fullpath_output_name=f'{pp_dir}/controlvar_pp.nc')
	ds = xr.open_dataset(f'{pp_dir}/controlvar_pp.nc')

if "histprocess" in sys.argv:
	daterange = np.arange(ASSIM_START_DATE,endtime+timedelta(hours=1),delta).astype(datetime)
	dates_string_array = [dateval.strftime("%Y%m%d_%H%M") for dateval in daterange]
	try:
		with open(f"{pp_dir}/bigY.pkl",'rb') as f:
			bigy=pickle.load(f)
	except FileNotFoundError:
		bigy = pt.makeYEachAssimPeriod(dates_string_array,useLevelEdge=useLevelEdge,useStateMet = useStateMet,useArea=useArea,use_numav=avtogcgrid,use_albedo=postprocess_save_albedo,fullpath_output_name=f"{pp_dir}/bigY.pkl")

if "calccol" in sys.argv:
	try:
		ds_level = xr.open_dataset(f'{pp_dir}/leveledge_pp.nc')
		ds_all = xr.open_dataset(f'{pp_dir}/controlvar_pp_all.nc')
	except FileNotFoundError:
		_ = pt.makeDatasetForEnsemble(ens_dir,controlvec,timeperiod,subset_rule = 'ALL',fullpath_output_name=f'{pp_dir}/controlvar_pp_all.nc')
		subdirs,_,_ = pt.globDirs(ens_dir,includeOutputDir=True)
		_ = pt.makeDatasetForDirectoryLevelEdge(subdirs[0],timeperiod,fullpath_output_name=f'{pp_dir}/controlvar_pp_all.nc')
		ds_level = xr.open_dataset(f'{pp_dir}/leveledge_pp.nc')
		ds_all = xr.open_dataset(f'{pp_dir}/controlvar_pp_all.nc')
	pressure_edges = ds_level[Met_PEDGE].values
	pressure_weights = (pressure_edges[:,1:,:,:]+pressure_edges[:,:-1,:,:])/2
	

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
	pt.plotSurfaceCell(ds,spec,30,59,outfile=f'{pp_dir}/wuhan_cell_ts_{spec}.png',includesNature=False)
	pt.plotSurfaceMean(ds,spec,outfile=f'{pp_dir}/surfmean_ts_{spec}.png',includesNature=False)
	if "calc850" in sys.argv:
		 pt.plotSurfaceMean(ds850,spec,outfile=f'{pp_dir}/mean850hPa_ts_{spec}.png',includesNature=False)

