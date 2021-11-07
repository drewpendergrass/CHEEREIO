from datetime import datetime
import toolbox as tx
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np

def read_tropomi(filename, species):
	"""
	Read TROPOMI data and save important variables to dictionary.

	Arguments
		filename [str]  : TROPOMI netcdf data file to read
		species [str]   : Species string 

	Returns
		met	  [dict] : Dictionary of important variables from TROPOMI:
							- species
 							- Latitude
							- Longitude
	   						- QA value
 							- UTC time
							- Averaging kernel
							- Prior profile
							- Dry air subcolumns
							- Latitude bounds
 							- Longitude bounds
							- Vertical pressure profile
	"""

	# Initialize dictionary for TROPOMI data
	met = {}
	
	# Store species, QA, lat, lon, time, averaging kernel
	data = xr.open_dataset(filename, group='PRODUCT')
	if species=='NO2':
		met[species] = data['nitrogendioxide_tropospheric_column'].values[0,:,:]
	elif species=='CH4':
		met[species] = data['methane_mixing_ratio_bias_corrected'].values[0,:,:]
	else:
		raise ValueError('Species not supported')

	met['qa_value'] = data['qa_value'].values[0,:,:]
	met['longitude'] = data['longitude'].values[0,:,:]
	met['latitude'] = data['latitude'].values[0,:,:]
	met['utctime'] = data['time_utc'].values[0,:]

	if species=='NO2':
		met['column_AK'] = data['averaging_kernel'].values[0,:,:,::-1]

	if species=='NO2':
		a = data['tm5_constant_a'].values[:,:]
		b = data['tm5_constant_b'].values[:,:]
	
	data.close()

	if species=='CH4':
		data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/DETAILED_RESULTS')
		met['column_AK'] = data['column_averaging_kernel'].values[0,:,:,::-1]
		data.close()

	# Store methane prior profile, dry air subcolumns
	data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
	if species=='CH4':
		pressure_interval = data['pressure_interval'].values[0,:,:]/100

	surface_pressure = data['surface_pressure'].values[0,:,:]/100				# Pa -> hPa
	data.close()

	# Store lat, lon bounds for pixels
	data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/GEOLOCATIONS')
	met['longitude_bounds'] = data['longitude_bounds'].values[0,:,:,:]
	met['latitude_bounds'] = data['latitude_bounds'].values[0,:,:,:]
	data.close()

	if species=='NO2':
		#Pressure levels (hpa) with dimension level, latind, lonind,edge (bottom/top)
		pressures = np.zeros((np.shape(b)[0],np.shape(surface_pressure)[0],np.shape(surface_pressure)[1],np.shape(b)[1]))
		for i in range(np.shape(surface_pressure)[0]):
			for j in range(np.shape(surface_pressure)[1]):
				pressures[:,i,j,:]=a+(b*surface_pressure[i,j])
	elif species=='CH4':
		N1 = met[species].shape[0]
		N2 = met[species].shape[1]
		pressures = np.zeros([N1,N2,13],dtype=np.float)
		pressures.fill(np.nan)
		for i in range(13):
			pressures[:,:,i]=surface_pressure-(i*pressure_interval)

		met['pressures'] = pressures
	
	return met

class TROPOMI_Translator(object):
	def __init__(self,testing=False):
		self.testing = testing
		self.spc_config = tx.getSpeciesConfig(self.testing)
		self.scratch = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/scratch"
	#Save dictionary of dates for later use
	def initialReadDate(self):
		sourcedirs = self.spc_config['TROPOMI_dirs']
		TROPOMI_date_dict = {}
		for key in list(sourcedirs.keys()):
			sourcedir = sourcedirs[key]
			obs_list = glob(f'{sourcedir}/S5P_*.nc')
			obs_list.sort()
			TROPOMI_date_dict[key] = {}
			TROPOMI_date_dict[key]['start'] = [datetime.strptime(obs.split('_')[-6], "%Y%m%dT%H%M%S") for obs in obs_list]
			TROPOMI_date_dict[key]['end'] = [datetime.strptime(obs.split('_')[-5], "%Y%m%dT%H%M%S") for obs in obs_list]
		with open(f"{self.scratch}/tropomi_dates.pickle", 'wb') as handle:
			pickle.dump(TROPOMI_date_dict, handle)
		return TROPOMI_date_dict
	#Timeperiod is two datetime objects
	def globObs(self,species,timeperiod):
		sourcedir = self.spc_config['TROPOMI_dirs'][species]
		if os.path.exists(f"{self.scratch}/tropomi_dates.pickle"):
			with open(f"{self.scratch}/tropomi_dates.pickle", 'rb') as handle:
				TROPOMI_date_dict = pickle.load(handle)
		else:
			TROPOMI_date_dict = self.initialReadDate()
		obs_dates = TROPOMI_date_dict[species]
		obs_list = glob(f'{sourcedir}/S5P_*.nc')
		obs_list.sort()
		obs_list = [obs for obs,t1,t2 in zip(obs_list,obs_dates['start'],obs_dates['end']) if (t1>=timeperiod[0]) and (t2<timeperiod[1])]
		return obs_list
	def getTROPOMI(self,species,timeperiod):
		obs_list = self.globObs(species,timeperiod)
		trop_obs = {}
		for obs in obs_list:
			trop_obs[obs]=read_tropomi(obs,species)
		return trop_obs
