#Some of this code is adapted from the excellent work done by Hannah Nesser, Melissa Sulprizio, Daniel Varon, Lu Shen
#and other contributors to the integrated methane inversion workflow. I am indebted to them and to Elise Penn and Alba Lorente for 
#explaining much of TROPOMI.

#TROPOMI CO operator is added by Sina Voshtani (2023/06/23)

from datetime import datetime
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np
import observation_operators as obsop

def read_tropomi(filename, species, filterinfo=None, includeObsError = False):
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

	if filename=="testNO2":
		filename = "/n/holylfs05/LABS/jacob_lab/dpendergrass/tropomi/NO2/2019/01/S5P_OFFL_L2__NO2____20190130T160807_20190130T174937_06729_01_010202_20190205T180152.nc"
	if filename=="testCH4":
		filename = "/n/holylfs05/LABS/jacob_lab/dpendergrass/tropomi/CH4/2019/01/S5P_OFFL_L2__CH4____20190131T223511_20190201T001641_06747_01_010202_20190207T000450.nc"
	# Initialize list for TROPOMI data
	met = {}
	
	# Store species, QA, lat, lon, time, averaging kernel
	data = xr.open_dataset(filename, group='PRODUCT')
	qa = data['qa_value'].values[0,:,:] #time,scanline,groundpixel
	if species=='NO2':
		sl,gp=np.where(qa>0.75)
		with xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA') as flag:
			with xr.open_dataset(filename, group = 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS') as flag2:
				with xr.open_dataset(filename, group = 'PRODUCT/SUPPORT_DATA/GEOLOCATIONS') as flag3:
					ocean_flag = flag['snow_ice_flag'].values[0,sl,gp]
					cloud_flag = flag2['cloud_radiance_fraction_nitrogendioxide_window'].values[0, sl, gp]
					sza_flag = flag3['solar_zenith_angle'].values[0, sl, gp]
					land_mask = np.where((ocean_flag == 0) & (cloud_flag < 0.3) & (sza_flag < 70))
					sl = sl[land_mask]
					gp = gp[land_mask]
	elif species=="CH4":
		sl,gp=np.where(qa>0.5)
	elif species=="CO":
		sl,gp=np.where(qa>0.75)
	else:
		raise ValueError('Species not supported')
	met['qa_value'] = qa[sl,gp]
	if species=='NO2':
		airmassfactor_trop = data['air_mass_factor_troposphere'].values[0,sl,gp]
		airmassfactor_total = data['air_mass_factor_total'].values[0,sl,gp]
		trop_layer_index = data['tm5_tropopause_layer_index'].values[0,sl,gp].astype(int)
		met[species] = data['nitrogendioxide_tropospheric_column'].values[0,sl,gp]
	elif species=='CH4':
		met[species] = data['methane_mixing_ratio_bias_corrected'].values[0,sl,gp] #time,scanline,groundpixel
	elif species=='CO':
		met[species] = data['carbonmonoxide_total_column_corrected'].values[0,sl,gp] #time,scanline,groundpixel
	else:
		raise ValueError('Species not supported')

	if includeObsError:
		if species=='NO2':
			met['Error'] = data['nitrogendioxide_tropospheric_column_precision'].values[0,sl,gp]
		elif species=='CH4':
			met['Error'] = data['methane_mixing_ratio_precision'].values[0,sl,gp]
		elif species=='CO':
			met['Error'] = data['carbonmonoxide_total_column_precision'].values[0,sl,gp]
	
	met['longitude'] = data['longitude'].values[0,sl,gp] #time,scanline,groundpixel
	met['latitude'] = data['latitude'].values[0,sl,gp] #time,scanline,groundpixel
	met['utctime'] = data['time_utc'].values[0,sl] #time, scanline

	if species=='NO2':
		met['column_AK'] = data['averaging_kernel'].values[0,sl,gp,::-1]
		#Multiply by total airmass over trop airmass, as in PUM 8.8 "Using the averaging kernel"
		met['column_AK'] = met['column_AK'] * (airmassfactor_total/airmassfactor_trop).reshape((len(airmassfactor_total),1))
		#set averaging kernel above tropopause to 0, as in PUM 8.8. Layer is 0-indexed
		for i in range(len(airmassfactor_total)):
			met['column_AK'][i,(trop_layer_index[i]+1):] = 0
		a = np.append(data['tm5_constant_a'].values[:,0],data['tm5_constant_a'].values[-1,1]) #layer, vertices (bottom/top)
		b = np.append(data['tm5_constant_b'].values[:,0],data['tm5_constant_b'].values[-1,1])

	data.close()

	if species=='CH4':
		data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/DETAILED_RESULTS')
		met['column_AK'] = data['column_averaging_kernel'].values[0,sl,gp,::-1] #time,scanline,groundpixel,layer
		met['albedo_swir'] = data['surface_albedo_SWIR'].values[0,sl,gp]
		met['albedo_nir'] = data['surface_albedo_NIR'].values[0,sl,gp]
		met['blended_albedo'] = (met['albedo_nir']*2.4)-(met['albedo_swir']*1.13)
		met['swir_aot'] = data['aerosol_optical_thickness_SWIR'].values[0,sl,gp]
		data.close()
	# Store CO averaging_kernel (missing in older version of TROPOMI CO!)
	if species=='CO':
		data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/DETAILED_RESULTS')
		met['column_AK'] = data['column_averaging_kernel'].values[0,sl,gp,::-1] #time,scanline,groundpixel,layer
		pressure_levels = data['pressure_levels'].values[0,sl,gp,::-1]/100 #time,scanline,groundpixel
		data.close()

	# Store methane prior profile, dry air subcolumns. Not needed for NO2, though surface pressure is
	if species=='CH4':
		data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
		met['methane_profile_apriori']=data['methane_profile_apriori'].values[0,sl,gp,::-1] #in mol/m2
		met['dry_air_subcolumns']=data['dry_air_subcolumns'].values[0,sl,gp,::-1] #in mol/m2
		met['surface_elevation_sd'] = data['surface_altitude_precision'].values[0,sl,gp]
		pressure_interval = data['pressure_interval'].values[0,sl,gp]/100 #time,scanline,groundpixel
		surface_pressure = data['surface_pressure'].values[0,sl,gp]/100 #time,scanline,groundpixel				# Pa -> hPa
		data.close()
	elif species=='NO2':
		data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
		surface_pressure = data['surface_pressure'].values[0,sl,gp] #time,scanline,groundpixel				# Leave Pa
		data.close()
	# Store CO prior profile (missing in older version of TROPOMI CO!)  
	if species=='CO':
		data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
		met['carbonmonoxide_profile_apriori']=data['carbonmonoxide_profile_apriori'].values[0,sl,gp,::-1] # in mol/m2
		met['surface_elevation']=data['surface_altitude'].values[0,sl,gp] #in m
		surface_pressure = data['surface_pressure'].values[0,sl,gp]/100 #time,scanline,groundpixel                              # Pa -> hPa
		data.close()


	# Store lat, lon bounds for pixels
	# data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/GEOLOCATIONS')
	# met['longitude_bounds'] = data['longitude_bounds'].values[0,:,:,:]
	# met['latitude_bounds'] = data['latitude_bounds'].values[0,:,:,:]
	# data.close()

	if species=='NO2':
		#Pressure levels (hpa) with dimension len(sl),levels
		pressures = np.zeros((len(sl),len(b)),dtype=float)
		pressures.fill(np.nan)
		for i in range(len(surface_pressure)):
			pressures[i,:]=(a+(b*surface_pressure[i]))/100 #Pa -> hPa
	elif species=='CH4':
		pressures = np.zeros([len(sl),13],dtype=np.float)
		pressures.fill(np.nan)
		for i in range(13):
			pressures[:,i]=surface_pressure-(i*pressure_interval)
	elif species=='CO':
		pressures = np.zeros([len(sl),51],dtype=float)
		pressures[:, :50] = pressure_levels
		pressures[:, 50]=0
	
	met['pressures'] = pressures
	
	if filterinfo is not None:
		met = obsop.apply_filters(met,filterinfo)

	return met


#Read the ACMG version of TROPOMI CH4
def read_tropomi_acmg(filename, species, filterinfo=None, includeObsError = False):

	if species !='CH4':
		raise ValueError('Species not supported.')
	# Initialize list for TROPOMI data
	met = {}
	
	# Store species, QA, lat, lon, time, averaging kernel
	data = xr.open_dataset(filename)
	qa = data['qa_value'].values
	goodvals=np.where(qa>0.5)[0]

	met['qa_value'] = qa[goodvals]
	met[species] = data['xch4_corrected'].values[goodvals] #nobs

	if includeObsError:
		met['Error'] = data['xch4_precision'].values[goodvals] #nobs
	
	met['longitude'] = data['longitude_center'].values[goodvals] #nobs
	met['latitude'] = data['latitude_center'].values[goodvals] #nobs
	timeraw = data['time'].values[goodvals,:] #nobs, ntime. Seven entries in format year', 'month', 'day', 'hour','minute', 'second'
	#format as CHEEREIO-compliant string
	timestring = [f'{str(int(timestamp[0]))}-{str(int(timestamp[1])).zfill(2)}-{str(int(timestamp[2])).zfill(2)}T{str(int(timestamp[3])).zfill(2)}:{str(int(timestamp[4])).zfill(2)}:{str(int(timestamp[5])).zfill(2)}Z' for timestamp in timeraw]
	met['utctime'] =  np.array(timestring)

	met['column_AK'] = data['xch4_column_averaging_kernel'].values[goodvals,::-1] #nobs,layer
	met['albedo_swir'] = data['surface_albedo'].values[goodvals,1] #nobs, nwin. 0 for nwin is NIR, 1 is swir
	met['albedo_nir'] = data['surface_albedo'].values[goodvals,0]
	met['blended_albedo'] = (met['albedo_nir']*2.4)-(met['albedo_swir']*1.13)
	met['swir_aot'] = data['aerosol_optical_thickness'].values[goodvals,1]

	#no surface elevation std, make sure it is set to nan in ens_config.
	met['methane_profile_apriori']=data['ch4_profile_apriori'].values[goodvals,::-1] #nobs,layer. in molec/cm2, but conversion factor divides out
	met['dry_air_subcolumns']=data['dry_air_subcolumns'].values[goodvals,::-1] #nobs,layer. in molec/cm2
	pressure_interval = data['dp'].values[goodvals] #nobs #already in hPa
	surface_pressure = data['surface_pressure'].values[goodvals] #nobs	#already in hPa
		
	data.close()

	pressures = np.zeros([len(goodvals),13],dtype=np.float) #nobs,layer
	pressures.fill(np.nan)
	for i in range(13):
		pressures[:,i]=surface_pressure-(i*pressure_interval)
	
	met['pressures'] = pressures
	
	if filterinfo is not None:
		met = obsop.apply_filters(met,filterinfo)

	return met



#Read the Belasus et al 2023 version of TROPOMI CH4, corrected using GOSAT
#no swir_aot, make sure it is set to nan in ens_config.
def read_tropomi_gosat_corrected(filename, species, filterinfo=None, includeObsError = False):

	if species !='CH4':
		raise ValueError('Species not supported.')
	# Initialize list for TROPOMI data
	met = {}
	
	# Store species, QA, lat, lon, time, averaging kernel
	data = xr.open_dataset(filename)

	qa = data['qa_value'].values

	landflag = np.bitwise_and(data['surface_classification'].values.astype(int),3) #Get landflag, apply bitwise and of mask value (3) to get 0=land,1=water,2=someland,3=coastline
	goodvals=np.where((qa>0.5) & (landflag==0))[0] #Only include land values (landflag=0); no coasts or ocean retrievals

	met['dry_air_subcolumns']=data['dry_air_subcolumns'].values[goodvals,::-1] #nobs,layer. in molec/cm2
	pressure_interval = data['pressure_interval'].values[goodvals]/100 #nobs #in Pa, convert to hPa
	surface_pressure = data['surface_pressure'].values[goodvals]/100 #nobs	#in Pa, convert to hPa
	met['surface_elevation_sd'] = data['surface_altitude_precision'].values[goodvals]

	met['qa_value'] = qa[goodvals]

	met[species] = data['methane_mixing_ratio_blended'].values[goodvals] #nobs

	if includeObsError:
		met['Error'] = data['methane_mixing_ratio_precision'].values[goodvals] #nobs

	met['methane_profile_apriori']=data['methane_profile_apriori'].values[goodvals,::-1] #nobs,layer. in molec/cm2, but conversion factor divides out
	met['column_AK'] = data['column_averaging_kernel'].values[goodvals,::-1] #nobs,layer

	met['longitude'] = data['longitude'].values[goodvals] #nobs
	met['latitude'] = data['latitude'].values[goodvals] #nobs
	met['utctime'] =  data['time_utc'].values[goodvals]

	met['albedo_swir'] = data['surface_albedo_SWIR'].values[goodvals] #nobs, nwin. 0 for nwin is NIR, 1 is swir
	met['albedo_nir'] = data['surface_albedo_NIR'].values[goodvals]
	met['blended_albedo'] = (met['albedo_nir']*2.4)-(met['albedo_swir']*1.13)

	pressures = np.zeros([len(goodvals),13],dtype=np.float) #nobs,layer
	pressures.fill(np.nan)
	for i in range(13):
		pressures[:,i]=surface_pressure-(i*pressure_interval)
	
	met['pressures'] = pressures
	
	if filterinfo is not None:
		met = obsop.apply_filters(met,filterinfo)

	return met

#This seems to be the memory bottleneck
#This is customized for CO due to the more serious memory bottleneck. It processes a chunk of observations at a time and concatenates the outputs,
#thus keeping the memory usage low depending on the "chunk_size". This can be adaptable to other species and especially in case of memory shortage.
def GC_to_sat_levels(GC_SPC, GC_edges, sat_edges, species, chunk_size=10000):
	'''
	The provided edges for GEOS-Chem and the satellite should
	have dimension number of observations x number of edges
	'''
	if species in ["CO", "NO2"]:
		nobs = GC_SPC.shape[0]
		GC_on_sat_list = []
		for k in range(0, nobs, chunk_size):
			GC_SPC_chunk = GC_SPC[k:k+chunk_size]
			GC_edges_chunk = GC_edges[k:k+chunk_size]
			sat_edges_chunk = sat_edges[k:k+chunk_size]
			# We want to account for the case when the GEOS-Chem surface
			# is above the satellite surface (altitude wise) or the GEOS-Chem
			# top is below the satellite top.. We do this by adjusting the
			# GEOS-Chem surface pressure up to the TROPOMI surface pressure
			idx_bottom = np.less(GC_edges_chunk[:, 0], sat_edges_chunk[:, 0])
			idx_top = np.greater(GC_edges_chunk[:, -1], sat_edges_chunk[:, -1])
			GC_edges_chunk[idx_bottom, 0] = sat_edges_chunk[idx_bottom, 0]
			GC_edges_chunk[idx_top, -1] = sat_edges_chunk[idx_top, -1]
			#Low mem calculation avoids allocating large arrays but runs slower.
			#This is accomplished by iterating over the nobs axis.
			#Default, higher memory calculation exploits vectorization for speed
			#But can result in enormous matrices for some experiments.
			# Define vectors that give the "low" and "high" pressure
			# values for each GEOS-Chem and satellite layer.
			GC_lo = GC_edges_chunk[:, 1:][:, :, None]
			GC_hi = GC_edges_chunk[:, :-1][:, :, None]
			sat_lo = sat_edges_chunk[:, 1:][:, None, :]
			sat_hi = sat_edges_chunk[:, :-1][:, None, :]
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
			GC_on_sat_chunk = (GC_to_sat*GC_SPC_chunk[:, :, None]).sum(axis=1)
			GC_to_sat_sum = GC_to_sat.sum(axis=1)
			mask = GC_to_sat_sum != 0
			GC_on_sat_chunk[mask] = GC_on_sat_chunk[mask] / GC_to_sat_sum[mask]
			#add to list
			GC_on_sat_list.append(GC_on_sat_chunk)
		#Concatenate list
		if len(GC_on_sat_list) > 0: #If there are observations
			GC_on_sat = np.concatenate(GC_on_sat_list, axis=0)
		else: #No observations
			GC_on_sat = np.zeros((nobs,sat_edges.shape[1]-1))
	else:
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
		#Low mem calculation avoids allocating large arrays but runs slower.
		#This is accomplished by iterating over the nobs axis.
		#Default, higher memory calculation exploits vectorization for speed
		#But can result in enormous matrices for some experiments. 
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
		GC_on_sat = (GC_to_sat*GC_SPC[:, :, None]).sum(axis=1)
		GC_on_sat = GC_on_sat/GC_to_sat.sum(axis=1)
	return GC_on_sat

def apply_avker(sat_avker, sat_pressure_weight, GC_SPC, sat_prior=None,filt=None):
	'''
	Apply the averaging kernel
	Inputs:
		sat_avker			The averaging kernel for the satellite
		sat_prior			The satellite prior profile in ppb, optional (used for CH4)
		sat_prior (CO)			The satellite prior profile in mol/m2, converted to ppb (used for CO)
		sat_pressure_weight  The relative pressure weights for each level
		GC_SPC			   The GC species on the satellite levels
		filt				 A filter, optional
	'''
	if filt is None:
		filt = np.ones(sat_avker.shape[1])
	else:
		filt = filt.astype(int)
	if sat_prior is None:
		#GC_col = (filt*sat_pressure_weight*sat_avker*GC_SPC)
		#no sat_pressure_weight version
		GC_col = filt*sat_avker*GC_SPC
	else:
		GC_col = (filt*sat_pressure_weight
				  *(sat_prior + sat_avker*(GC_SPC - sat_prior)))
	GC_col = GC_col.sum(axis=1)
	return GC_col 

class TROPOMI_Translator(obsop.Observation_Translator):
	def __init__(self,verbose=1):
		super().__init__(verbose)
	#Save dictionary of dates for later use
	def initialReadDate(self):
		sourcedirs = self.spc_config['TROPOMI_dirs']
		TROPOMI_date_dict = {}
		for key in list(sourcedirs.keys()):
			sourcedir = sourcedirs[key]
			obs_list = glob(f'{sourcedir}/**/S5P_*.nc', recursive=True)
			obs_list.sort()
			TROPOMI_date_dict[key] = {}
			TROPOMI_date_dict[key]['start'] = [datetime.strptime(obs.split('_')[-6], "%Y%m%dT%H%M%S") for obs in obs_list]
			TROPOMI_date_dict[key]['end'] = [datetime.strptime(obs.split('_')[-5], "%Y%m%dT%H%M%S") for obs in obs_list]
		with open(f"{self.scratch}/tropomi_dates.pickle", 'wb') as handle:
			pickle.dump(TROPOMI_date_dict, handle)
		return TROPOMI_date_dict
	#Timeperiod is two datetime objects
	def globObs(self,species,timeperiod, interval=None):
		sourcedir = self.spc_config['TROPOMI_dirs'][species]
		if os.path.exists(f"{self.scratch}/tropomi_dates.pickle"):
			with open(f"{self.scratch}/tropomi_dates.pickle", 'rb') as handle:
				TROPOMI_date_dict = pickle.load(handle)
		else:
			TROPOMI_date_dict = self.initialReadDate()
		obs_dates = TROPOMI_date_dict[species]
		obs_list = glob(f'{sourcedir}/**/S5P_*.nc', recursive=True)
		obs_list.sort()
		if interval:
			obs_list = [obs for obs,t1,t2 in zip(obs_list,obs_dates['start'],obs_dates['end']) if (t1>=timeperiod[0]) and (t2<timeperiod[1]) and ((t1.hour % interval == 0) or (t1.hour % interval == (interval-1)))]
		else:
			obs_list = [obs for obs,t1,t2 in zip(obs_list,obs_dates['start'],obs_dates['end']) if (t1>=timeperiod[0]) and (t2<timeperiod[1])]
		return obs_list
	def getObservations(self,specieskey,timeperiod, interval=None, includeObsError=False):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		obs_list = self.globObs(species,timeperiod,interval)
		trop_obs = []
		filterinfo = {}
		if species=='CH4':
			if (self.spc_config['Extensions']['TROPOMI_CH4']=="True") and (self.spc_config['TROPOMI_CH4_FILTERS']=="True"): #Check first if extension is on before doing the TROPOMI filtering
				filterinfo["TROPOMI_CH4"] = [float(self.spc_config['TROPOMI_CH4_filter_blended_albedo']),float(self.spc_config['TROPOMI_CH4_filter_swir_albedo_low']),float(self.spc_config['TROPOMI_CH4_filter_swir_albedo_high']),float(self.spc_config['TROPOMI_CH4_filter_winter_lat']),float(self.spc_config['TROPOMI_CH4_filter_roughness']),float(self.spc_config['TROPOMI_CH4_filter_swir_aot'])]
		elif species=='CO':
			if (self.spc_config['Extensions']['TROPOMI_CO']=="True") and (self.spc_config['TROPOMI_CO_FILTERS']=="True"): #Check first if extension is on before doing the TROPOMI filtering
				filterinfo["TROPOMI_CO"] = [float(self.spc_config['TROPOMI_CO_filter_blended_albedo']),float(self.spc_config['TROPOMI_CO_filter_swir_albedo_low']),float(self.spc_config['TROPOMI_CO_filter_swir_albedo_high']),float(self.spc_config['TROPOMI_CO_filter_winter_lat']),float(self.spc_config['TROPOMI_CO_filter_roughness']),float(self.spc_config['TROPOMI_CO_filter_swir_aot'])]
		if specieskey in list(self.spc_config["filter_obs_poleward_of_n_degrees"].keys()):
			filterinfo['MAIN']=[float(self.spc_config["filter_obs_poleward_of_n_degrees"][specieskey])]
		for obs in obs_list:
			if self.spc_config[f'WHICH_TROPOMI_{species}_PRODUCT'] == 'ACMG':
				trop_obs.append(read_tropomi_acmg(obs,species,filterinfo,includeObsError=includeObsError))
			elif self.spc_config[f'WHICH_TROPOMI_{species}_PRODUCT'] == 'BLENDED':
				trop_obs.append(read_tropomi_gosat_corrected(obs,species,filterinfo,includeObsError=includeObsError))
			else:
				trop_obs.append(read_tropomi(obs,species,filterinfo,includeObsError=includeObsError))
		met = {}
		for key in list(trop_obs[0].keys()):
			met[key] = np.concatenate([metval[key] for metval in trop_obs])
		return met
	def gcCompare(self,specieskey,TROPOMI,GC,GC_area=None,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		extra_obsdata_to_save = self.spc_config['EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y'][specieskey]
		TROP_PW = (-np.diff(TROPOMI['pressures'])/(TROPOMI['pressures'][:, 0] - TROPOMI['pressures'][:, -1])[:, None])
		returnStateMet = self.spc_config['SaveStateMet']=='True'
		GC_col_data = obsop.getGCCols(GC,TROPOMI,species,self.spc_config,returninds=True,returnStateMet=returnStateMet,GC_area=GC_area)
		GC_SPC = GC_col_data['GC_SPC']
		GC_P = GC_col_data['GC_P']
		i,j,t = GC_col_data['indices']
		if species=='CH4':
			GC_SPC*=1e9 #scale to mol/mol
			TROP_PRIOR = 1e9*(TROPOMI['methane_profile_apriori']/TROPOMI['dry_air_subcolumns'])
			synthetic_partial_columns = False
		elif species=='NO2':
			TROP_PRIOR=None
			AIRMOL_VOL=GC_col_data['Met_AIRDEN'] / 0.028964 # get model layer dry air density in mol/m3 (air molar mass: 0.028964 kg/mol), Met_AIRDEN(kg/m3)
			AIRMOL_M2=(AIRMOL_VOL*GC_col_data['Met_BXHEIGHT']) # convert to column mol/m2 of dry air
			GC_SPC = GC_SPC*AIRMOL_M2 # convert GEOS-CHEM NO2 from VVdry to mol/m2. Scale up by 1.
			synthetic_partial_columns = True
		elif species=='CO':
			GC_SPC*=1e9 #scale to mol/mol
			TROP_P0=TROPOMI['pressures'][:,0]
			TROP_z0=TROPOMI['surface_elevation'] #  (m)
			TROP_T0=GC_col_data['Met_T'][:,0]
			TROP_z = np.zeros_like(TROPOMI['carbonmonoxide_profile_apriori'])
			for l in range(50):
				TROP_z[:,l]= TROP_z0 + 500 + l*1000
			TROP_P_mid=TROPOMI['pressures']+np.diff(TROPOMI['pressures'],axis=1,append=0)/2.0
			TROP_P_mid=TROP_P_mid[:,0:-1]
			TROP_PRIOR_mol = (TROPOMI['carbonmonoxide_profile_apriori']) # mole/m2
			AIRMOL_VOL=GC_col_data['Met_AIRDEN'] / 0.028964 # get model layer dry air density in mol/m3 (air molar mass: 0.028964 kg/mol), Met_AIRDEN(kg/m3)
			AIRMOL_COL=(AIRMOL_VOL*GC_col_data['Met_BXHEIGHT']).sum(axis=1) # convert to column mol/m2 of dry air
			TROP_PRIOR_ppb = ((TROP_PRIOR_mol) / (np.repeat(AIRMOL_COL[:, np.newaxis], TROP_PW.shape[1], axis=1)*TROP_PW)) *1e9 # convert prior in mol/m2 to ppbv
			TROP_PRIOR = TROP_PRIOR_ppb
			TROPOMI[species]=1e9*TROPOMI[species]/AIRMOL_COL # convert TROPOMI CO from mol/m2 to ppbv
			synthetic_partial_columns = False
		#If super observations, prep the error calculation
		if self.spc_config['AV_TO_GC_GRID'][specieskey]=="True":
			superobs=True
			superObsFunction = self.spc_config['SUPER_OBSERVATION_FUNCTION'][specieskey]
			additional_args_avgGC = {}
			if doErrCalc:
				if useObserverError:
					if species=='CO':
						TROPOMI['Error']=1e9*TROPOMI['Error']/AIRMOL_COL # convert tropomi errro from mol/m2 to ppbv
					additional_args_avgGC['obsInstrumentError'] = TROPOMI['Error']
				elif prescribed_error is not None:
					additional_args_avgGC['prescribed_error'] = prescribed_error
					additional_args_avgGC['prescribed_error_type'] = prescribed_error_type
				if minError is not None:
					additional_args_avgGC['minError'] = minError
				if errorCorr is not None:
					additional_args_avgGC['errorCorr'] = errorCorr
				if transportError is not None:
					additional_args_avgGC['modelTransportError'] = transportError
			#If saving extra fields, add them here
			if len(extra_obsdata_to_save)>0:
				additional_args_avgGC['other_fields_to_avg'] = {}
				for field in extra_obsdata_to_save:
					additional_args_avgGC['other_fields_to_avg'][field] = TROPOMI[field]
		else:
			superobs=False
		#If user is taking super observations and wants to average before applying averaging kernel, do so now
		if superobs and (self.spc_config['SUPER_OBS_BEFORE_Hx'][specieskey]=="True"):
			presuperobs=True
			trop_agg_args = {'GC_SPC':GC_SPC,'GC_P':GC_P, 'T_press':TROPOMI['pressures'], 'ak':TROPOMI['column_AK'],'TROP_PW':TROP_PW}
			if TROP_PRIOR is not None: #only include prior in averaging if there is a prior
				trop_agg_args['TROP_PRIOR']=TROP_PRIOR
			agg_data = obsop.averageFieldsToGC(i,j,t,GC,TROPOMI[species],doSuperObs=doErrCalc,superObsFunction=superObsFunction,**additional_args_avgGC,**trop_agg_args)
			GC_on_sat_l = GC_to_sat_levels(agg_data['GC_SPC'], agg_data['GC_P'], agg_data['T_press'],species)
			#If trop prior is provided, use it in the averaging kernel. otherwise, pass in None.
			if TROP_PRIOR is not None: 
				tp = agg_data['TROP_PRIOR']
			else:
				tp = None
			GC_on_sat = apply_avker(agg_data['ak'],agg_data['TROP_PW'], GC_on_sat_l,tp)
		else:
			presuperobs=False
			GC_on_sat_l = GC_to_sat_levels(GC_SPC, GC_P, TROPOMI['pressures'],species)
			GC_on_sat = apply_avker(TROPOMI['column_AK'],TROP_PW, GC_on_sat_l,TROP_PRIOR)
		nan_indices = np.argwhere(np.isnan(GC_on_sat))
		GC_on_sat = np.nan_to_num(GC_on_sat)
		if superobs:
			if presuperobs:
				toreturn = obsop.ObsData(GC_on_sat,agg_data['obs'],agg_data['lat'],agg_data['lon'], agg_data['time'],num_av=agg_data['num'])
				if doErrCalc:
					toreturn.addData(err_av=agg_data['err'])
				if 'additional_fields' in agg_data:
					toreturn.addData(**agg_data['additional_fields'])
			else:
				toreturn = obsop.averageByGC(i,j,t,GC,GC_on_sat,TROPOMI[species],doSuperObs=doErrCalc,superObsFunction=superObsFunction,**additional_args_avgGC)
		else:
			timevals = GC.time.values[t]
			toreturn = obsop.ObsData(GC_on_sat,TROPOMI[species],TROPOMI['latitude'],TROPOMI['longitude'],timevals)
			#If saving extra fields, add them here
			if len(extra_obsdata_to_save)>0:
				data_to_add = {}
				for field in extra_obsdata_to_save:
					data_to_add[field] = TROPOMI[field]
				toreturn.addData(**data_to_add)
			if doErrCalc and useObserverError:
				toreturn.addData(err_av=TROPOMI['Error'])
		return toreturn


