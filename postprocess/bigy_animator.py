import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
import pickle
import sys 
sys.path.append('../core')
import settings_interface as si 

data = si.getSpeciesConfig()
gclat,gclon = si.getLatLonVals(data)
gclat = np.array(gclat)
gclon = np.array(gclon)
pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"

anim_fps = int(data['animation_fps_scalingfactor'])

with open(f'{pp_dir}/bigy_arrays_for_plotting.pkl','rb') as f:
	bigy_arrays=pickle.load(f)

filename_helper = {
	'sim_obs':'simulated_observations',
	'obs':'actual_observations',
	'obscount':'total_raw_observation_counts',
	'obscount_avg':'total_aggregated_observation_counts'
}

for species in bigy_arrays:
	#We only plot gridded data here.
	if bigy_arrays[species]['interpret_as'] == 'points':
		continue
	fields_to_plot = ['sim_obs','obs','obscount','obscount_avg'] + data['EXTRA_OBSDATA_FIELDS_TO_REGRID_AND_PLOT'][species]
	labels = [f"{data['OBSERVED_SPECIES'][species]} ({data['OBSERVATION_UNITS'][species]})"]*2+['Count']*2
	for bonus_field in data['EXTRA_OBSDATA_FIELDS_TO_REGRID_AND_PLOT'][species]:
		labels.append(data['extra_plot_field_units'][bonus_field])
	m = Basemap(projection='cyl', resolution='l',llcrnrlat=-90, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180)
	for field,label in zip(fields_to_plot,labels):
		if field not in bigy_arrays[species]:
			continue #Skip if field is not present
		arrayval = bigy_arrays[species][field]
		dates = bigy_arrays[species]['dates']
		if field in filename_helper:
			outfile = f'{pp_dir}/{species}_{filename_helper[field]}.mp4'
		else:
			outfile = f'{pp_dir}/{species}_{field}.mp4'
        animateData(m, data=arrayval, file_out=outfile, lon=gclon, lat=gclat, anim_fps = anim_fps, timestr = dates)
