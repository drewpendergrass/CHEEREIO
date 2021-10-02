import json
import xarray as xr
import postprocess_tools as pt
from glob import glob
from datetime import datetime,timedelta
import sys 

with open('../ens_config.json') as f:
	data = json.load(f)

pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"
ens_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs"
controlvec = data['CONTROL_VECTOR_CONC']
statevec = data['STATE_VECTOR_CONC']
emisvec = data['CONTROL_VECTOR_EMIS']
ASSIM_START_DATE=datetime.strptime(data['ASSIM_START_DATE'], "%Y%m%d")
endtime=datetime.strptime(data['END_DATE'], "%Y%m%d")
ASSIM_TIME=data['ASSIM_TIME']
delta = timedelta(hours=int(ASSIM_TIME))
starttime = ASSIM_START_DATE-(2*delta)
timeperiod = (starttime,endtime)

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

if "calc850" in sys.argv:
	print('Calculating/loading 850hPa pressure level')
	try:
		ds850 = xr.open_dataset(f'{pp_dir}/controlvar_pp_850hPa.nc')
	except FileNotFoundError:
		_ = pt.makeDatasetForEnsemble(ens_dir,controlvec,timeperiod,subset_rule="850",fullpath_output_name=f'{pp_dir}/controlvar_pp_850hPa.nc')
		ds850 = xr.open_dataset(f'{pp_dir}/controlvar_pp_850hPa.nc')

for spec in controlvec:
	pt.plotSurfaceCell(ds,spec,30,59,outfile=f'{pp_dir}/wuhan_cell_ts_{spec}.png',includesNature=True)
	pt.plotSurfaceMean(ds,spec,outfile=f'{pp_dir}/surfmean_ts_{spec}.png',includesNature=True)
	if "calc850" in sys.argv:
		 pt.plotSurfaceMean(ds850,spec,outfile=f'{pp_dir}/mean850hPa_ts_{spec}.png',includesNature=True)

