import json
import xarray as xr
import postprocess_tools as pt
from glob import glob
from datetime import datetime,timedelta

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
starttime = ASSIM_START_DATE-delta
endtime = ASSIM_START_DATE+delta
timeperiod = (starttime,endtime)

pt.combineScaleFactors(ens_dir,pp_dir)
scalefactor_files = glob(f'{pp_dir}/*_SCALEFACTOR.nc')
for scalefactor in scalefactor_files:
	sf_name = '_'.join(scalefactor.split('/')[-1].split('_')[0:-1])
	pt.plotEmissionsCell(scalefactor,30,59,outfile=f'{pp_dir}/wuhan_cell_emis_{sf_name}.png')

try:
	ds = xr.open_dataset(f'{pp_dir}/controlvar_pp.nc')
except FileNotFoundError:
	_ = pt.makeDatasetForEnsemble(ens_dir,controlvec,timeperiod,hourlysub=1,fullpath_output_name=f'{pp_dir}/controlvar_pp.nc')
	ds = xr.open_dataset(f'{pp_dir}/controlvar_pp.nc')

for spec in controlvec:
	#pt.plotSurfaceCellEnsMeanNorm(ds,spec,30,59,outfile=f'{pp_dir}/wuhan_cell_ts_{spec}_zeromean.png',unit='ppm')
	pt.plotSurfaceCell(ds,spec,30,59,outfile=f'{pp_dir}/wuhan_cell_ts_{spec}.png',unit='ppm',includesNature=True)

