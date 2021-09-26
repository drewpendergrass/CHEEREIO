import json
import xarray as xr
import postprocess_tools as pt

with open('../ens_config.json') as f:
	data = json.load(f)

pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"
ens_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs"
controlvec = data['CONTROL_VECTOR_CONC']
statevec = data['STATE_VECTOR_CONC']

try:
	ds = xr.open_dataset(f'{pp_dir}/controlvar_pp.nc')
except FileNotFoundError:
	ds = pt.makeDatasetForEnsemble(ens_dir,controlvec,fullpath_output_name=f'{pp_dir}/controlvar_pp.nc')

controlvec = ['NO','NO2']
for spec in controlvec:
	pt.plotSurfaceCellEnsMeanNorm(ds,spec,30,59,outfile=f'{pp_dir}/wuhan_cell_ts_{spec}_zeromean.png',unit='ppm')
	#pt.plotSurfaceCell(ds,spec,30,59,outfile=f'{pp_dir}/wuhan_cell_ts_{spec}.png',unit='ppm',includesNature=True,nature_error=0.2,natureErrType='relative')
