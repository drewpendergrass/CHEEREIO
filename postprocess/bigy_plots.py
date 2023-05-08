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
scalefactor_plot_freq = data['scalefactor_plot_freq']
hemco_diags_to_process = data['hemco_diags_to_process']
min_emis_value_to_plot = data['min_emis_value_to_plot']
omit_diff_cells_with_fewer_than_n_observations = int(data['omit_diff_cells_with_fewer_than_n_observations'])
if len(min_emis_value_to_plot) > 0: #If it's not none
	min_emis_value_to_plot = float(min_emis_value_to_plot) #parse to float
	if np.isnan(min_emis_value_to_plot): #If it's a nan, store as None for plotting function
		min_emis_value_to_plot = None
else:
	min_emis_value_to_plot = None

min_emis_std_value_to_plot = data['min_emis_std_value_to_plot'] 
if len(min_emis_std_value_to_plot) > 0: #If it's not none
	min_emis_std_value_to_plot = float(min_emis_std_value_to_plot) #parse to float
	if np.isnan(min_emis_std_value_to_plot): #If it's a nan, store as None for plotting function
		min_emis_std_value_to_plot = None
else:
	min_emis_std_value_to_plot = None

useControl=data['DO_CONTROL_RUN']=="true"
lognormalErrors=data['lognormalErrors']=="true"
useLogScaleForEmissionsMaps = data['useLogScaleForEmissionsMaps']=="True"
anim_fps = int(data['animation_fps_scalingfactor'])

with open(f'{pp_dir}/bigy_arrays_for_plotting.pkl','rb') as f:
	pickledata=pickle.load(f)

regridded_bigy = regridBigYdata(pickledata,gclat,gclon)

m = Basemap(projection='cyl', resolution='l',llcrnrlat=-90, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180)

#Plot observation means and counts
#Plot assimilation minus obs and ctrl minus obs
for species in regridded_bigy:
	assim_minus_obs = regridded_bigy[species]['assim_minus_obs']
	ctrl_minus_obs = regridded_bigy[species]['ctrl_minus_obs']
	total_obs_in_period=regridded_bigy[species]['total_obs_in_period']
	if regridded_bigy[species]['interpret_as'] == 'map':
		plotMap(m,gclat,gclon,total_obs_in_period,species,f'{pp_dir}/total_obs_count_{species}.png',useLog=True)
		plotMap(m,gclat,gclon,regridded_bigy[species]['total_weighted_mean_true_obs'],species,f'{pp_dir}/weighted_mean_obs_{species}.png') 
		#Remove pixels that are too low in observation (tend to be noisy)
		pixels_to_remove = np.where(total_obs_in_period<omit_diff_cells_with_fewer_than_n_observations)
		assim_minus_obs[pixels_to_remove[0],pixels_to_remove[1]]=np.nan 
		ctrl_minus_obs[pixels_to_remove[0],pixels_to_remove[1]]=np.nan 
		clim_abs = np.max([np.nanmax(np.abs(assim_minus_obs)),np.nanmax(np.abs(ctrl_minus_obs))])
		plotMap(m,gclat,gclon,assim_minus_obs,species,f'{pp_dir}/assim_minus_obs_{species}.png',cmap=plt.cm.seismic,clim = [-1*clim_abs,clim_abs])
		plotMap(m,gclat,gclon,ctrl_minus_obs,species,f'{pp_dir}/ctrl_minus_obs_{species}.png',cmap=plt.cm.seismic,clim = [-1*clim_abs,clim_abs]) 
	elif regridded_bigy[species]['interpret_as'] == 'points':
		clim_abs = np.max([np.nanmax(np.abs(assim_minus_obs)),np.nanmax(np.abs(ctrl_minus_obs))])
		plotMapPoints(m, regridded_bigy[species]['lat'], regridded_bigy[species]['lon'], total_obs_in_period, species,f'{pp_dir}/total_obs_count_{species}.png',useLog=True)
		plotMapPoints(m, regridded_bigy[species]['lat'], regridded_bigy[species]['lon'], regridded_bigy[species]['mean_obs'], species,f'{pp_dir}/mean_obs_{species}.png')
		plotMapPoints(m, regridded_bigy[species]['lat'], regridded_bigy[species]['lon'], assim_minus_obs, species,f'{pp_dir}/assim_minus_obs_{species}.png',cmap=plt.cm.seismic,clim = [-1*clim_abs,clim_abs])
		plotMapPoints(m, regridded_bigy[species]['lat'], regridded_bigy[species]['lon'], ctrl_minus_obs, species,f'{pp_dir}/ctrl_minus_obs_{species}.png',cmap=plt.cm.seismic,clim = [-1*clim_abs,clim_abs])
	print(f'For species {species} we have, for assimilation minus observations, a mean of {np.nanmean(assim_minus_obs)} and a standard deviation of {np.nanstd(assim_minus_obs)}')
	print(f'For species {species} we have, for control minus observations, a mean of {np.nanmean(ctrl_minus_obs)} and a standard deviation of {np.nanstd(ctrl_minus_obs)}')


if scalefactor_plot_freq == 'monthly':
	aggToMonthly=True
else:
	aggToMonthly=False 

#Plot scale factor slices
plotScaleFactor(m,gclat,gclon,pp_dir, useLognormal = lognormalErrors, aggToMonthly=aggToMonthly)
#Plot emission slices
plotEmissions(m,gclat,gclon,pp_dir,hemco_diags_to_process=hemco_diags_to_process,plotWithLogScale=useLogScaleForEmissionsMaps, min_emis=min_emis_value_to_plot,min_emis_std=min_emis_std_value_to_plot, useLognormal = lognormalErrors, plotcontrol=useControl, aggToMonthly=aggToMonthly)

