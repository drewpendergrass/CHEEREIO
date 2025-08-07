#Some of this code is adapted from Hannah Nesser.

from datetime import datetime,timedelta
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np
import pandas as pd
import observation_operators as obsop
import settings_interface as si 

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
		#Save out the site code
		site_code = data.attrs['site_code']
		data['site_code'] = xr.DataArray([site_code]*len(data.obs), dims=('obs'))
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
	inds_to_drop = np.where(obspack['CT_sampling_strategy'].values==999)[0]
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
		daily['site_code'] = daily['site_code'].astype('S50')
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
	return data[['obspack_id', 'value', 'altitude', 'latitude', 'longitude', 'time', 'platform','site_code']]

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
		if len(obs_list)>0:
			obspack = xr.open_mfdataset(obs_list, concat_dim='obs', combine='nested', mask_and_scale=False, preprocess=filter_postprocess_obspack_from_file)
			met = {}
			met['value'] = obspack['value'].values
			met['longitude'] = obspack['longitude'].values
			met['latitude'] = obspack['latitude'].values
			met['altitude'] = obspack['altitude'].values
			met['utctime'] = obspack['time'].values
			met['platform'] = obspack['platform'].values
			met['obspack_id'] = obspack['obspack_id'].values
			met['site_code'] = obspack['site_code'].values
		else: #If no observations, return something empty
			met = {}
			met['value'] = np.array([])
			met['longitude'] = np.array([])
			met['latitude'] = np.array([])
			met['altitude'] = np.array([])
			met['utctime'] = np.array([])
			met['platform'] = np.array([])
			met['obspack_id'] = np.array([])
			met['site_code'] = np.array([])
		return met
	def gcCompare(self,specieskey,ObsPack,GC,GC_area=None,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		#Have to convert to PPB
		if species=='CH4':
			obs_multiplier = 1e9
			gc_multiplier=1e9
		elif species=='N2O':
			obs_multiplier = 1
			gc_multiplier=1e9
		if len(ObsPack['value'])==0:
			toreturn = obsop.ObsData(np.array([]),ObsPack['value']*obs_multiplier,ObsPack['latitude'],ObsPack['longitude'],ObsPack['utctime'])
			toreturn.addData(altitude=ObsPack['altitude'],pressure=np.array([]),obspack_id=ObsPack['obspack_id'],platform=ObsPack['platform'],site_code=ObsPack['site_code'])
		else:
			gc_overrides = si.checkGCSpeciesOverride(self.spc_config) #Here we check to see if GC stores the species under a different name. Only applies for N2O.
			if (len(gc_overrides)>0)&(species in gc_overrides):
				gc_sim = GC[gc_overrides[species]].values*gc_multiplier
			else:
				gc_sim = GC[species].values*gc_multiplier
			# Apply error filter from the extension. Remove observations where error exceeds some (very high) maximum threshhold as being likely in error.
			if (self.spc_config['Extensions']['Obspack_N2O']=="True"):
				maxerr = float(self.spc_config['Max_Obspack_N2O_Error'])
				if ~np.isnan(maxerr):
					err = np.abs(gc_sim-ObsPack['value']*obs_multiplier)
					valid = np.where(err<maxerr)[0]
					gc_sim = gc_sim[valid]
					obsval = ObsPack['value'][valid]*obs_multiplier
					obslat = ObsPack['latitude'][valid]
					obslon = ObsPack['longitude'][valid]
					obstime = ObsPack['utctime'][valid]
					obsalt = ObsPack['altitude'][valid]
					gcpressure = GC['pressure'].values[valid]
					obsid = ObsPack['obspack_id'][valid]
					obsplat = ObsPack['platform'][valid]
					obssite = ObsPack['site_code'][valid]
				else:
					obsval = ObsPack['value']*obs_multiplier
					obslat = ObsPack['latitude']
					obslon = ObsPack['longitude']
					obstime = ObsPack['utctime']
					obsalt = ObsPack['altitude']
					gcpressure = GC['pressure'].values
					obsid = ObsPack['obspack_id']
					obsplat = ObsPack['platform']
					obssite = ObsPack['site_code']
			toreturn = obsop.ObsData(gc_sim,obsval,obslat,obslon,obstime)
			toreturn.addData(altitude=obsalt,pressure=gcpressure,obspack_id=obsid,platform=obsplat,site_code=obssite)
		return toreturn

