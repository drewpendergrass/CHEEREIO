import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm,LinearSegmentedColormap
import pickle
from datetime import datetime
from glob import glob

#Plot 2d field over a map
def plotMap(m,lat,lon,flat,labelname,outfile,clim=None,cmap=None,useLog=False,minval = None):
	fig = plt.figure(figsize=(10, 6))
	m.drawcountries(color='lightgray')
	m.drawcoastlines(color='lightgray')
	if cmap is None:
		cmap = plt.cm.jet
	if useLog:
		if minval is not None:
			flat[flat<=minval] = np.nan #user can optionally supply minimum value for log plots; anything below is not shown
		else:
			flat[flat<=0] = np.nan
		mesh = m.pcolormesh(lon, lat, flat,latlon=True,cmap=cmap,norm=LogNorm())
	else:
		mesh = m.pcolormesh(lon, lat, flat,latlon=True,cmap=cmap)
	if clim is not None:
		plt.clim(clim[0],clim[1])
	else:
		plt.clim(np.nanmin(flat), np.nanmax(flat))
	plt.colorbar(label=labelname)
	fig.savefig(outfile)

#Plot points on a map, with color codes
def plotMapPoints(m, lat, lon, zvals, labelname,outfile,clim=None,cmap=None,useLog=False,minval = None):
	fig = plt.figure(figsize=(10, 6))
	m.drawcountries(color='lightgray')
	m.drawcoastlines(color='lightgray')
	if cmap is None:
		cmap = plt.cm.jet
	if useLog:
		if minval is not None:
			zvals[zvals<=minval] = np.nan #user can optionally supply minimum value for log plots; anything below is not shown
		else:
			zvals[zvals<=0] = np.nan
		mesh = m.scatter(lon, lat, c=zvals,latlon=True,cmap=cmap,norm=LogNorm(),edgecolors='black')
	else:
		mesh = m.scatter(lon, lat, c=zvals,latlon=True,cmap=cmap,edgecolors='black')
	if clim is not None:
		plt.clim(clim[0],clim[1])
	else:
		plt.clim(np.nanmin(zvals), np.nanmax(zvals))
	plt.colorbar(label=labelname)
	fig.savefig(outfile)


def plotEmissions(m,lat,lon,ppdir, hemco_diags_to_process,plotWithLogScale=True, min_emis=None,min_emis_std=None, clim=None, clim_std=None,plotcontrol=True,useLognormal = False, aggToMonthly=True,conversion_factor=None):
	hemcodiag = xr.open_dataset(f'{ppdir}/combined_HEMCO_diagnostics.nc')
	if plotcontrol:
		hemcocontroldiag = xr.open_dataset(f'{ppdir}/control_HEMCO_diagnostics.nc')
	dates = hemcodiag['time'].values
	for diag in hemco_diags_to_process:
		hemcofield = hemcodiag[diag].values
		if conversion_factor is not None:
			hemcofield*=conversion_factor
		if plotcontrol:
			ctrlfield = hemcocontroldiag[diag].values
			if conversion_factor is not None:
				ctrlfield*=conversion_factor
		if len(np.shape(hemcofield)) == 5: #we have emissions at higher levels (e.g. aircraft)
			hemcofield = np.sum(hemcofield,axis=2) #ensemble gone, dim 0 is ens, dim 1 is time, dim 2 is lev. Sum up.
			if plotcontrol:
				ctrlfield = np.sum(ctrlfield,axis=1) #dim 0 is time, dim 1 is lev.
		if aggToMonthly:
			dates,scalar = agg_to_monthly(dates, hemcofield)
			_,ctrlfield = agg_to_monthly(dates, ctrlfield)
		if useLognormal: 
			hemcofield_std = np.exp(np.std(np.log(hemcofield),axis=0)) #std of log will give the lognormal shape parameter; exponentiate back into emissions space. 
			hemcofield = np.exp(np.mean(np.log(hemcofield),axis=0)) #average across ensemble
		else:
			hemcofield_std = np.std(hemcofield,axis=0) #standard deviation across ensemble
			hemcofield = np.mean(hemcofield,axis=0) #average across ensemble
		#Now hemcofield is of dim time, lat, lon
		timelabels = [str(timeval)[0:13] for timeval in dates]
		#Do the plotting.
		#Get colormaps
		if clim_std is None:
			if min_emis_std is not None:
				clim_std = [min_emis_std,np.max(hemcofield_std)]
			else:
				clim_std = [np.min(hemcofield_std),np.max(hemcofield_std)]
		if clim is None:
			if plotcontrol:
				if min_emis is not None:
					clim  = [min_emis, np.max([np.max(hemcofield),np.max(ctrlfield)])]
				else:
					clim  = [np.min([np.min(hemcofield),np.min(ctrlfield)]), np.max([np.max(hemcofield),np.max(ctrlfield)])]
			else:
				if min_emis is not None:
					clim  = [min_emis, np.max(hemcofield)]
				else:
					clim  = [np.min(hemcofield), np.max(hemcofield)]
		cmap = plt.cm.jet
		for i,dateval in enumerate(timelabels):
			plotMap(m,lat,lon,hemcofield[i,:,:],diag,f'{ppdir}/{diag}_{dateval}_ensemble_mean.png',clim = clim, useLog=plotWithLogScale,minval = min_emis)
			plotMap(m,lat,lon,hemcofield_std[i,:,:],diag,f'{ppdir}/{diag}_{dateval}_ensemble_std.png',clim = clim_std, useLog=plotWithLogScale,minval = min_emis_std)
			if plotcontrol:
				plotMap(m,lat,lon,ctrlfield[i,:,:],diag,f'{ppdir}/{diag}_{dateval}_control.png',clim = clim, useLog=plotWithLogScale,minval = min_emis)


def plotScaleFactor(m,lat,lon,ppdir, useLognormal = False, aggToMonthly=True,plot_on_log_scale=False,clim = None):
	files = glob(f'{ppdir}/*_SCALEFACTOR.nc')
	files.sort()
	sf_names = [pts.split('/')[-1][0:-15] for pts in files]
	for file,name in zip(files,sf_names):
		ds = xr.load_dataset(file)
		dates = ds['time'].values
		scalar = ds['Scalar'].values
		if aggToMonthly:
			dates,scalar = agg_to_monthly(dates, scalar)
		if useLognormal:
			scalar = np.exp(np.mean(np.log(scalar),axis=0)) #average across ensemble
		else:
			scalar = np.mean(scalar,axis=0) #average across ensemble
		timelabels = [str(timeval)[0:13] for timeval in dates]
		if plot_on_log_scale:
			cmap = plt.cm.bwr #Let user supplied clim and plotMap do the rest.
		else:
			#Make custom blue-white-red colorbar centered at one
			cvals  = [0.0, 1.0, np.max([np.max(scalar),1.1])]
			colors = ["blue","white","red"]
			pltnorm=plt.Normalize(min(cvals),max(cvals))
			tuples = list(zip(map(pltnorm,cvals), colors))
			cmap = LinearSegmentedColormap.from_list("", tuples)
			if clim is None:
				clim = [0,np.max([np.max(scalar),1.1])]
		for i,dateval in enumerate(timelabels):
			plotMap(m,lat,lon,scalar[i,:,:],'Scaling factor',f'{ppdir}/{name}_{dateval}_scalefactor.png',clim=clim,cmap=cmap,useLog=plot_on_log_scale)


#Takes regridded bigy data; goes through and calculates key fields to plot (OmF etc)
#Handles both gridded and station data, though treatments of each are distinct
#Control is required for this function
def regridBigYdata(bigy,gclat,gclon,timeperiod=None):
	to_return = {}
	for species in bigy:
		to_return[species] = {}
		spec_dict = bigy[species]
		to_return[species]['interpret_as'] = spec_dict['interpret_as']
		#Handle point data
		if spec_dict['interpret_as'] == 'points':
			num_stations = len(list(spec_dict.keys()))-1 #Number of stations in dataset
			#Create empty arrays to fill later
			to_return[species]['stations'] = [] #List of stations
			to_return[species]['lat'] = np.zeros(num_stations)*np.nan #station lat
			to_return[species]['lon'] = np.zeros(num_stations)*np.nan #station lon
			to_return[species]['assim_minus_obs'] = np.zeros(num_stations)*np.nan #posterior minus observations
			to_return[species]['ctrl_minus_obs'] = np.zeros(num_stations)*np.nan #prior minus observations
			to_return[species]['total_obs_in_period'] = np.zeros(num_stations)*np.nan #total number of observations
			to_return[species]['mean_obs'] = np.zeros(num_stations)*np.nan #mean obs in period
			counter = 0
			#Loop through stations
			for station in spec_dict:
				if station == 'interpret_as':
					continue
				else:
					station_dict = spec_dict[station]
					times = station_dict['time']
					if timeperiod is not None: #Slice data down to timeperiod
						inds = np.where((times >= timeperiod[0]) & (times < timeperiod[1]))[0]
						sim_obs = station_dict["sim_obs"][inds]
						true_obs = station_dict["obs"][inds]
						ctrl_obs = station_dict["control"][inds]
						to_return[species]['total_obs_in_period'][counter]  = len(inds)
					else:
						sim_obs = station_dict["sim_obs"]
						true_obs = station_dict["obs"]
						ctrl_obs = station_dict["control"]
						to_return[species]['total_obs_in_period'][counter]  = len(ctrl_obs)
					to_return[species]['stations'].append(station)
					to_return[species]['lat'][counter] = station_dict['lat'][0] #Should be all identical
					to_return[species]['lon'][counter] = station_dict['lon'][0] #Should be all identical
					to_return[species]['assim_minus_obs'][counter]  = np.nanmean(sim_obs-true_obs)
					to_return[species]['ctrl_minus_obs'][counter]  = np.nanmean(ctrl_obs-true_obs)
					to_return[species]['mean_obs'][counter]  = np.nanmean(true_obs)
					counter+=1
			to_return[species]['stations'] = np.array(to_return[species]['stations'])
		elif spec_dict['interpret_as'] == 'map':
			dates = spec_dict["dates"]
			datevals = [datetime.strptime(dateval,'%Y%m%d_%H%M') for dateval in dates]
			total_satellite_obs=spec_dict["obscount"]
			true_obs = spec_dict["obs"]
			sim_obs = spec_dict["sim_obs"]
			ctrl_obs = spec_dict["control"]
			if timeperiod is not None: #Slice data down to timeperiod
				inds = [i for i, e in enumerate(datevals) if (e >= timeperiod[0]) & (e < timeperiod[1])]
				total_satellite_obs=total_satellite_obs[inds,:,:]
				true_obs = true_obs[inds,:,:]
				sim_obs = sim_obs[inds,:,:]
				ctrl_obs = ctrl_obs[inds,:,:]
			total_obs_in_period = np.sum(total_satellite_obs,axis=0)
			to_return[species]['total_obs_in_period'] = total_obs_in_period
			to_return[species]['total_weighted_mean_true_obs'] = np.zeros(np.shape(total_obs_in_period))
			to_return[species]['assim_minus_obs'] = np.zeros(np.shape(total_obs_in_period))*np.nan
			to_return[species]['ctrl_minus_obs'] = np.zeros(np.shape(total_obs_in_period))*np.nan
			for j in range(len(gclat)):
				for k in range(len(gclon)):
					if np.sum(~np.isnan(true_obs[:,j,k]))>0:
						to_return[species]['assim_minus_obs'][j,k] = np.nanmean(sim_obs[:,j,k]-true_obs[:,j,k])
						to_return[species]['ctrl_minus_obs'][j,k] = np.nanmean(ctrl_obs[:,j,k]-true_obs[:,j,k])
					if np.sum(total_satellite_obs[:,j,k]) == 0:
						to_return[species]['total_weighted_mean_true_obs'][j,k] = np.nan
					else:
						indices = np.where(np.logical_not(np.isnan(true_obs[:,j,k])))[0]
						to_return[species]['total_weighted_mean_true_obs'][j,k] = np.average(true_obs[indices,j,k],weights=total_satellite_obs[indices,j,k])
	return to_return

def agg_to_monthly(dates, to_agg):
	agg_dim = len(np.shape(to_agg))
	years = dates.astype('datetime64[Y]').astype(int) + 1970
	months = dates.astype('datetime64[M]').astype(int) % 12 + 1
	ym = (years*100) + months
	_, ind = np.unique(ym,return_index = True)
	dates = dates[ind]
	shape = np.array(np.shape(to_agg))
	if agg_dim == 3:
		shape[0] = len(ind)
	else:
		shape[1] = len(ind)
	to_return = np.zeros(shape)
	for i in range(len(ind)):
		if agg_dim == 3:
			if i == len(ind)-1:
				subset = to_agg[ind[i]:,:,:]
			else:
				subset = to_agg[ind[i]:ind[i+1],:,:]
			to_return[i,:,:] = np.mean(subset,axis=0)
		else:
			if i == len(ind)-1:
				subset = to_agg[:,ind[i]:,:,:]
			else:
				subset = to_agg[:,ind[i]:ind[i+1],:,:]
			to_return[:,i,:,:] = np.mean(subset,axis=1)
	return [dates,to_return]
