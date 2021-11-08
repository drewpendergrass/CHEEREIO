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

	# Initialize list for TROPOMI data
	met = {}
	
	# Store species, QA, lat, lon, time, averaging kernel
	data = xr.open_dataset(filename, group='PRODUCT')
	qa = data['qa_value'].values[0,:,:] #time,scanline,groundpixel
	sl,gp=np.where(qa>0.5)
	met['qa_value'] = qa[sl,gp]
	if species=='NO2':
		met[species] = data['nitrogendioxide_tropospheric_column'].values[0,:,:]
	elif species=='CH4':
		met[species] = data['methane_mixing_ratio_bias_corrected'].values[0,sl,gp] #time,scanline,groundpixel
	else:
		raise ValueError('Species not supported')

	
	met['longitude'] = data['longitude'].values[0,sl,gp] #time,scanline,groundpixel
	met['latitude'] = data['latitude'].values[0,sl,gp] #time,scanline,groundpixel
	met['utctime'] = data['time_utc'].values[0,sl] #time, scanline

	if species=='NO2':
		met['column_AK'] = data['averaging_kernel'].values[0,:,:,::-1]
		a = data['tm5_constant_a'].values[:,:]
		b = data['tm5_constant_b'].values[:,:]
	
	data.close()

	if species=='CH4':
		data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/DETAILED_RESULTS')
		met['column_AK'] = data['column_averaging_kernel'].values[0,sl,gp,::-1] #time,scanline,groundpixel,layer
		data.close()

	# Store methane prior profile, dry air subcolumns
	data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
	if species=='CH4':
		met['methane_profile_apriori']=data['methane_profile_apriori'].values[0,sl,gp,::-1]
		met['dry_air_subcolumns']=data['dry_air_subcolumns'].values[0,sl,gp,::-1]
		pressure_interval = data['pressure_interval'].values[0,sl,gp]/100 #time,scanline,groundpixel

	surface_pressure = data['surface_pressure'].values[0,sl,gp]/100 #time,scanline,groundpixel				# Pa -> hPa
	data.close()

	# Store lat, lon bounds for pixels
	# data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/GEOLOCATIONS')
	# met['longitude_bounds'] = data['longitude_bounds'].values[0,:,:,:]
	# met['latitude_bounds'] = data['latitude_bounds'].values[0,:,:,:]
	# data.close()

	if species=='NO2':
		#Pressure levels (hpa) with dimension level, latind, lonind,edge (bottom/top)
		pressures = np.zeros((np.shape(b)[0],np.shape(surface_pressure)[0],np.shape(surface_pressure)[1],np.shape(b)[1]))
		for i in range(np.shape(surface_pressure)[0]):
			for j in range(np.shape(surface_pressure)[1]):
				pressures[:,i,j,:]=a+(b*surface_pressure[i,j])
	elif species=='CH4':
		pressures = np.zeros([len(sl),13],dtype=np.float)
		pressures.fill(np.nan)
		for i in range(13):
			pressures[:,i]=surface_pressure-(i*pressure_interval)
		met['pressures'] = pressures
	
	return met

def nearest_loc(GC,TROPOMI):
	# Find the grid box and time indices corresponding to TROPOMI obs
	# i index
	iGC = np.abs(GC.lon.values.reshape((-1, 1))
				 - TROPOMI['longitude'].reshape((1, -1)))
	iGC = iGC.argmin(axis=0)
	# j index
	jGC = np.abs(GC.lat.values.reshape((-1, 1))
				 - TROPOMI['latitude'].reshape((1, -1)))
	jGC = jGC.argmin(axis=0)
	# Time index
	tGC = np.abs(GC.time.values.reshape((-1, 1))
				 - TROPOMI['utctime'].astype('datetime64').reshape((1, -1)))
	tGC = tGC.argmin(axis=0)
	return iGC, jGC, tGC

def getGCCols(GC,TROPOMI,species):
	i,j,t = nearest_loc(GC,TROPOMI)
	return [GC[f'SpeciesConc_{species}'].values[t,:,j,i],GC[f'Met_PEDGE'].values[t,:,j,i]]

def GC_to_sat_levels(GC_CH4, GC_edges, sat_edges):
	'''
	The provided edges for GEOS-Chem and the satellite should
	have dimension number of observations x number of edges
	'''
	# We want to account for the case when the GEOS-Chem surface
	# is above the satellite surface (altitude wise) or the GEOS-Chem
	# top is below the satellite top.. We do this by adjusting the
	# GEOS-Chem surface pressure up to the TROPOMI surface pressure
	idx_bottom = np.less(GC_edges[:, 0], sat_edges[:, 0])
	idx_top = np.greater(GC_edges[:, -1], sat_edges[:, -1])
	GC_edges[idx_bottom, 0] = sat_edges[idx_bottom, 0]
	GC_edges[idx_top, -1] = sat_edges[idx_top, -1]
	# Define vectors that give the "low" and "high" pressure
	# values for each GEOS-Chem and satellite layer.
	GC_lo = GC_edges[:, 1:][:, :, None]
	GC_hi = GC_edges[:, :-1][:, :, None]
	sat_lo = sat_edges[:, 1:][:, None, :]
	sat_hi = sat_edges[:, :-1][:, None, :]
	# Get the indices where the GC-to-satellite mapping, which is
	# a nobs x ngc x nsat matrix, is non-zero
	idx = (np.less_equal(sat_lo, GC_hi) & np.greater_equal(sat_hi, GC_lo))
	# Find the fraction of each GC level that contributes to each
	# TROPOMI level. We should first divide (to normalize) and then
	# multiply (to apply the map to the column) by the GC pressure
	# difference, but we exclude this (since it's the same as x1).
	GC_to_sat = np.minimum(sat_hi, GC_hi) - np.maximum(sat_lo, GC_lo)
	GC_to_sat[~idx] = 0
	# Now map the GC CH4 to the satellite levels
	GC_on_sat = (GC_to_sat*GC_CH4[:, :, None]).sum(axis=1)
	GC_on_sat = GC_on_sat/GC_to_sat.sum(axis=1)
	return GC_on_sat

def apply_avker(sat_avker, sat_prior, sat_pressure_weight, GC_CH4, filt=None):
	'''
	Apply the averaging kernel
	Inputs:
		sat_avker			The averaging kernel for the satellite
		sat_prior			The satellite prior profile in ppb
		sat_pressure_weight  The relative pressure weights for each level
		GC_CH4			   The GC methane on the satellite levels
		filt				 A filter, optional
	'''
	if filt is None:
		filt = np.ones(sat_avker.shape[1])
	else:
		filt = filt.astype(int)
	GC_col = (filt*sat_pressure_weight
			  *(sat_prior + sat_avker*(GC_CH4 - sat_prior)))
	GC_col = GC_col.sum(axis=1)
	return GC_col 

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
	def getSatellite(self,species,timeperiod):
		obs_list = self.globObs(species,timeperiod)
		trop_obs = []
		for obs in obs_list:
			trop_obs.append(read_tropomi(obs,species))
		met = {}
		for key in list(trop_obs[0].keys()):
			met[key] = np.concatenate([metval[key] for metval in trop_obs])
		return met
	def gcCompare(self,species,timeperiod,TROPOMI,GC):
		TROP_CH4 = 1e9*(TROPOMI['methane_profile_apriori']/TROPOMI['dry_air_subcolumns'])
		TROP_PW = (-np.diff(TROPOMI['pressures'])/(TROPOMI['pressures'][:, 0] - TROPOMI['pressures'][:, -1])[:, None])
		GC_CH4,GC_P = getGCCols(GC,TROPOMI,species)
		if species=='CH4':
			GC_CH4*=1e9 #scale to ppb
		GC_on_sat = GC_to_sat_levels(GC_CH4, GC_P, TROPOMI['pressures'])
		GC_on_sat = apply_avker(TROPOMI['column_AK'],TROP_CH4, TROP_PW, GC_on_sat)
		return [GC_on_sat,TROPOMI[species],TROPOMI['latitude'],TROPOMI['longitude']]


