import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from map_tools import *
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
postprocess_save_albedo = data['postprocess_save_albedo']=="True"

with open(f'{pp_dir}/bigy_arrays_for_plotting.pkl','rb') as f:
	pickledata=pickle.load(f)

dates = pickledata["dates"]
specieslist = pickledata["species"]
total_satellite_obs=pickledata["obscount"]
total_averaged_obs=pickledata["obscount_avg"]
true_obs = pickledata["obs"]
sim_obs = pickledata["sim_obs"]
ctrl_obs = pickledata["control"]

m = Basemap(projection='cyl', resolution='l',llcrnrlat=-90, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180)

total_obs_in_period = np.sum(total_satellite_obs,axis=0)
total_weighted_mean_true_obs = np.zeros(np.shape(total_obs_in_period))
assim_minus_obs = np.zeros(np.shape(total_obs_in_period))*np.nan
ctrl_minus_obs = np.zeros(np.shape(total_obs_in_period))*np.nan

for i,species in enumerate(specieslist):
	for j in range(len(gclat)):
		for k in range(len(gclon)):
			if np.sum(~np.isnan(true_obs[:,i,j,k]))>0:
				assim_minus_obs[i,j,k] = np.nanmean(sim_obs[:,i,j,k]-true_obs[:,i,j,k])
				ctrl_minus_obs[i,j,k] = np.nanmean(ctrl_obs[:,i,j,k]-true_obs[:,i,j,k])
			if np.sum(total_satellite_obs[:,i,j,k]) == 0:
				total_weighted_mean_true_obs[i,j,k] = np.nan
			else:
				total_weighted_mean_true_obs[i,j,k] = np.average(true_obs[:,i,j,k],weights=total_satellite_obs[:,i,j,k])

#Plot observation means and counts
#Plot assimilation minus obs and ctrl minus obs
for i,species in enumerate(specieslist):
	plotMap(m,gclat,gclon,total_obs_in_period[i,:,:],species,f'{pp_dir}/total_obs_count_{species}.png',useLog=True)
	plotMap(m,gclat,gclon,total_weighted_mean_true_obs[i,:,:],species,f'{pp_dir}/weighted_mean_obs_{species}.png') 
	clim_abs = np.max([np.nanmax(np.abs(assim_minus_obs[i,:,:])),np.nanmax(np.abs(ctrl_minus_obs[i,:,:]))])
	plotMap(m,gclat,gclon,assim_minus_obs[i,:,:],species,f'{pp_dir}/assim_minus_obs_{species}.png',cmap=plt.cm.seismic,clim = [-1*clim_abs,clim_abs])
	print(f'For species {species} we have, for assimilation minus observations, a mean of {np.nanmean(assim_minus_obs[i,:,:])} and a standard deviation of {np.nanstd(assim_minus_obs[i,:,:])}')
	plotMap(m,gclat,gclon,ctrl_minus_obs[i,:,:],species,f'{pp_dir}/ctrl_minus_obs_{species}.png',cmap=plt.cm.seismic,clim = [-1*clim_abs,clim_abs]) 
	print(f'For species {species} we have, for control minus observations, a mean of {np.nanmean(ctrl_minus_obs[i,:,:])} and a standard deviation of {np.nanstd(ctrl_minus_obs[i,:,:])}')


#Plot scale factor slices

plotScaleFactor(m,gclat,gclon,pp_dir, plotMonthStartOnly=True)
