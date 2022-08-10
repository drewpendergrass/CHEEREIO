import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from map_tools import *
import pickle
import json

with open('../ens_config.json') as f:
	data = json.load(f)

pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"

anim_fps = int(data['animation_fps_scalingfactor'])
postprocess_save_albedo = data['postprocess_save_albedo']=="True"

with open(f'{pp_dir}/bigy_arrays_for_plotting.pkl','rb') as f:
	pickledata=pickle.load(f)

dates = pickledata[2]
specieslist = pickledata[3]
total_satellite_obs=pickledata[0]
total_averaged_obs=pickledata[1]
true_obs = pickledata[4]
sim_obs = pickledata[5]

with open(f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch/latlon_vals.json") as f:
	ll_data = json.load(f)

gclat = np.array(ll_data['lat'])
gclon = np.array(ll_data['lon'])

m = Basemap(projection='cyl', resolution='l',llcrnrlat=-90, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180)

total_obs_in_period = np.sum(total_satellite_obs,axis=0)
total_weighted_mean_true_obs = np.zeros(np.shape(total_obs_in_period))

for i,species in enumerate(specieslist):
	for j in range(len(gclat)):
		for k in range(len(gclon)):
			total_weighted_mean_true_obs[i,j,k] = np.average(true_obs[:,i,j,k],weights=total_satellite_obs[:,i,j,k])

#Plot observation means and counts

for i,species in enumerate(specieslist):
	plotmap(m,gclat,gclon,total_obs_in_period[i,:,:],species,f'total_obs_count_{species}.png',useLog=True)
	plotmap(m,gclat,gclon,total_weighted_mean_true_obs[i,:,:],species,f'weighted_mean_obs_{species}.png') 

