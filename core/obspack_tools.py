#Some of this code is adapted from Hannah Nesser.

from datetime import datetime,timedelta
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np
import pandas as pd
import observation_operators as obsop

data_vars = ['time', 'start_time', 'midpoint_time', 'time_components', 'value',
'latitude', 'longitude', 'altitude', 'assimilation_concerns',
'obspack_id']

#Make a filter function (output function) with start and end date and 2-value (min,max) tuples for lat and lon bounds
def make_filter_fxn(start_date,end_date,lat_bounds=None,lon_bounds=None):
	# Define a filtering function
	def filter_obspack(data):
	# Define the variables to keep
		# Subset variables
		data = data[data_vars]
		# Subset for time and location
		data = data.where((data['time'].dt.date >= start_date.date()) & (data['time'].dt.date <= end_date.date()), drop=True)
		if lat_bounds is not None:
			data = data.where((data['latitude'] >= lat_bounds[0]) & (data['latitude'] <= lat_bounds[1]), drop=True)
		if lon_bounds is not None:
			data = data.where((data['longitude'] >= lon_bounds[0]) & (data['longitude'] <= lon_bounds[1]),drop=True)
		# Save out a platform variable
		platform = data.attrs['dataset_project'].split('-')[0]
		data['platform'] = xr.DataArray([platform]*len(data.obs), dims=('obs'))
		# Correct to local timezone if it's an in situ or surface observation
		if (len(data.obs) > 0) and (platform in ['surface', 'tower']):
			utc_conv = data.attrs['site_utc2lst']
			if int(utc_conv) != utc_conv:
				print('UTC CONVERSION FACTOR IS NOT AN INTEGER : ', data.attrs['dataset_name'])
			data['utc_conv'] = xr.DataArray(utc_conv*np.ones(len(data.obs)),dims=('obs'))
			# data['time_ltc'] = dc(data['time']) + np.timedelta64(int(utc_conv), 'h')
		else:
			data['utc_conv'] = xr.DataArray(np.zeros(len(data.obs)), dims=('obs'))
			# data['time_ltc'] = dc(data['time'])
		return data
	return filter_obspack

#Prepare raw obspack data for input into GEOSChem
#Take raw obspack data from raw_obspack_dir --> process from start_date to end_date --> output to gc_obspack_dir
def prep_obspack(raw_obspack_dir,gc_obspack_dir,filename_format,start_date,end_date):
	# Get a list of the files
	files = glob(f'{raw_obspack_dir}/*.nc')
	files = [f for f in files if f.split('/')[-1][:11] != 'obspack_ch4']
	files.sort()
	#Make filter function
	filter_obspack = make_filter_fxn(start_date,end_date)
	## Iterate through the files and see which are relevant to the domain
	filtered_files = []
	for i, f in enumerate(files):
		op = xr.open_dataset(f)
		# Only use files in the needed time, latitude, and longitude
		# ranges
		try:
			op = filter_obspack(op)
		except ValueError:
			continue
		except KeyError:
			print(f)
		# If the file is empty, continue through the loop
		if len(op.obs) == 0:
			continue
		# If the file still has observations, append it to conus_files
		filtered_files.append(op)
	# Now combine all the files
	obspack = xr.concat(filtered_files,dim='obs')
	# Check for the sampling strategy
	## Get the time in hours of each sample
	obspack['obs_length'] = (obspack['time'] - obspack['start_time'])
	obspack['obs_length'] = obspack['obs_length'].dt.seconds*2/(60*60)
	## Convert that to the sampling strategy flag
	## ss = place holder for sampling strategy
	obspack['ss'] = xr.DataArray(999*np.ones(len(obspack.obs)), dims=('obs'))
	## Closest to 4 hours
	obspack['ss'] = obspack['ss'].where(obspack['obs_length'] > 5.25, 1)
	## Closest to 90 minutes
	obspack['ss'] = obspack['ss'].where(obspack['obs_length'] > 2.75, 3)
	## Closest to 1 hour
	obspack['ss'] = obspack['ss'].where(obspack['obs_length'] > 1.25, 2)
	## Closest to instantaneous
	obspack['ss'] = obspack['ss'].where(obspack['obs_length'] > 0.5, 4)
	## Cast to int
	obspack['ss'] = obspack['ss'].astype(int)
	# Rename and add attributes
	obspack = obspack.rename({'ss' : 'CT_sampling_strategy'})
	obspack['CT_sampling_strategy'].attrs = {'_FillValue' : -9,'long_name' : 'model sampling strategy','values' : 'How to sample model. 1=4-hour avg; 2=1-hour avg; 3=90-min avg; 4=instantaneous'}
	# Other clean up
	obspack.attrs = {}
	obspack = obspack.drop(['obs_length', 'start_time', 'midpoint_time'])
	#Drop data without sampling strategy
	inds_to_drop = np.where(obspack['ss']==999)[0]
	obspack=obspack.drop_isel(obs=inds_to_drop)
	# And iterate through the unique days
	name_str = filename_format.split('YYYYMMDD.nc')[0]
	delta = end_date - start_date   # returns timedelta
	for i in range(delta.days + 1):
		day = start_date + timedelta(days=i)
		# Subset for that day
		daily = obspack.where((obspack['time'].dt.month == day.month) & (obspack['time'].dt.day == day.day) & (obspack['time'].dt.year == day.year), drop=True)
		# If there is no data, continue
		if len(daily.obs) == 0:
			continue
		# Data type fix
		daily['obspack_id'] = daily['obspack_id'].astype('S200')
		daily['platform'] = daily['platform'].astype('S50')
		# Time fix
		daily['time'].encoding['units'] = 'seconds since 1970-01-01 00:00:00 UTC'
		daily['time'].encoding['calendar'] = 'proleptic_gregorian'
		# daily['time_ltc'].encoding['units'] = 'seconds since 1970-01-01 00:00:00 UTC'
		# daily['time_ltc'].encoding['calendar'] = 'proleptic_gregorian'
		# One last rename
		# daily = daily.rename({'value' : 'obs'})
		# Otherwise, save out
		print(f'Saving {day.year:04d}-{day.month:02d}-{day.day:02d}')
		daily.to_netcdf(f'{gc_obspack_dir}/{name_str}{day.year:04d}{day.month:02d}{day.day:02d}.nc',unlimited_dims=['obs'])

def filter_postprocess_obspack_from_file(data):
	return data[['obspack_id', 'value', 'altitude', 'latitude', 'longitude', 'time', 'utc_conv', 'platform']]

class ObsPack_Translator(obsop.Observation_Translator):
	def __init__(self,verbose=1):
		super().__init__(verbose)
	#Save dictionary of dates for later use
	def initialReadDate(self):
		sourcedir = self.spc_config['gc_obspack_path']
		obs_list = glob(f'{sourcedir}/*.nc')
		obs_list.sort()
		obs_dates = [datetime.strptime(obs.split('/')[-1][-11:-3], "%Y%m%d") for obs in obs_list]
		with open(f"{self.scratch}/obspack_dates.pickle", 'wb') as handle:
			pickle.dump(obs_dates, handle)
		return obs_dates
	#Timeperiod is two datetime objects
	def globObs(self,species,timeperiod, interval=None):
		sourcedir = self.spc_config['gc_obspack_path']
		if os.path.exists(f"{self.scratch}/obspack_dates.pickle"):
			with open(f"{self.scratch}/obspack_dates.pickle", 'rb') as handle:
				obs_dates = pickle.load(handle)
		else:
			obs_dates = self.initialReadDate()
		obs_list = glob(f'{sourcedir}/*.nc')
		obs_list.sort()
		if interval:
			obs_list = [obs for obs,t in zip(obs_list,obs_dates) if (t>=timeperiod[0]) and (t<timeperiod[1]) and ((t.hour % interval == 0) or (t.hour % interval == (interval-1)))]
		else:
			obs_list = [obs for obs,t in zip(obs_list,obs_dates) if (t>=timeperiod[0]) and (t<timeperiod[1])]
		return obs_list
	def getObservations(self,specieskey,timeperiod, interval=None, includeObsError=False):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		obs_list = self.globObs(species,timeperiod,interval)
		obspack = xr.open_mfdataset(obs_list, concat_dim='obs', combine='nested', mask_and_scale=False, preprocess=filter_postprocess_obspack_from_file)
		met = {}
		met['value'] = obspack['value'].values
		met['longitude'] = obspack['longitude'].values
		met['latitude'] = obspack['latitude'].values
		met['altitude'] = obspack['altitude'].values
		met['utctime'] = obspack['time'].values
		met['utc_conv'] = obspack['utc_conv'].values
		met['platform'] = obspack['platform'].values
		met['obspack_id'] = obspack['obspack_id'].values
		return met
	def gcCompare(self,specieskey,ObsPack,GC,GC_area=None,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		toreturn = obsop.ObsData(GC[species].values,ObsPack['value'],ObsPack['latitude'],ObsPack['longitude'],ObsPack['utctime'])
		toreturn.addData(utc_conv=ObsPack['utc_conv'],altitude=ObsPack['altitude'],pressure=GC['pressure'].values,obspack_id=ObsPack['obspack_id'],platform=ObsPack['platform'])
		return toreturn

