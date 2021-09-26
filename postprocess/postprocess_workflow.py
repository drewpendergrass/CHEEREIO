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
except OSError:
	ds = pt.makeDatasetForEnsemble(ens_dir,controlvec,fullpath_output_name=f'{pp_dir}/controlvar_pp.nc')

for spec in controlvec:
	pt.plotSurfaceCell(ds,spec,30,59,outfile=f'{pp_dir}/wuhan_cell_ts_{spec}.png')