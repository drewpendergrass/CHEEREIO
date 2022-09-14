#Some of this code is adapted from the excellent work done by Hannah Nesser, Melissa Sulprizio, Daniel Varon, Lu Shen
#and other contributors to the integrated methane inversion workflow. I am indebted to them and to Elise Penn and Alba Lorente for 
#explaining much of TROPOMI.

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
	elif species=="CH4":
		sl,gp=np.where(qa>0.5)
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
	else:
		raise ValueError('Species not supported')

	if includeObsError:
		if species=='NO2':
			met['Error'] = data['nitrogendioxide_tropospheric_column_precision'].values[0,sl,gp]
		elif species=='CH4':
			met['Error'] = data['methane_mixing_ratio_precision'].values[0,sl,gp]
	
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

	# Store methane prior profile, dry air subcolumns. Not needed for NO2, though surface pressure is
	if species=='CH4':
		data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
		met['methane_profile_apriori']=data['methane_profile_apriori'].values[0,sl,gp,::-1]
		met['dry_air_subcolumns']=data['dry_air_subcolumns'].values[0,sl,gp,::-1]
		met['surface_elevation_sd'] = data['surface_altitude_precision'].values[0,sl,gp]
		pressure_interval = data['pressure_interval'].values[0,sl,gp]/100 #time,scanline,groundpixel
		surface_pressure = data['surface_pressure'].values[0,sl,gp]/100 #time,scanline,groundpixel				# Pa -> hPa
		data.close()
	elif species=='NO2':
		data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
		surface_pressure = data['surface_pressure'].values[0,sl,gp] #time,scanline,groundpixel				# Leave Pa
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
	
	met['pressures'] = pressures
	
	if filterinfo is not None:
		met = obsop.apply_filters(met,filterinfo)

	return met


#This seems to be the memory bottleneck
def GC_to_sat_levels(GC_SPC, GC_edges, sat_edges, GC_M = None,lowmem=False,batchsize=None):
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
	if lowmem:
		satshape = np.shape(sat_edges)
		finalshape = (satshape[0],satshape[1]-1)
		if GC_M is not None:
			GC_M_on_sat = np.zeros(finalshape)
		GC_on_sat = np.zeros(finalshape)
		if (batchsize is not None) and (batchsize>1):
			numiter = int(np.ceil(np.shape(GC_edges)[0]/batchsize))
			for i in range(numiter):
				startvalue = int(i*batchsize)
				endvalue = int(np.min((np.shape(GC_edges)[0],(i+1)*batchsize)))
				# Define vectors that give the "low" and "high" pressure
				# values for each GEOS-Chem and satellite layer.
				GC_lo = GC_edges[startvalue:endvalue, 1:][:, :, None]
				GC_hi = GC_edges[startvalue:endvalue, :-1][:, :, None]
				sat_lo = sat_edges[startvalue:endvalue, 1:][:, None, :]
				sat_hi = sat_edges[startvalue:endvalue, :-1][:, None, :]
				idx = (np.less_equal(sat_lo, GC_hi) & np.greater_equal(sat_hi, GC_lo))
				GC_to_sat = np.minimum(sat_hi, GC_hi) - np.maximum(sat_lo, GC_lo)
				GC_to_sat[~idx] = 0
				# Now map the GC CH4 to the satellite levels
				GC_on_sat_subset = (GC_to_sat*GC_SPC[startvalue:endvalue, :, None]).sum(axis=1)
				GC_on_sat_subset = GC_on_sat_subset/GC_to_sat.sum(axis=1)
				GC_on_sat[startvalue:endvalue,:] = GC_on_sat_subset
				if GC_M is not None:
					GC_M_on_sat_subset = (GC_to_sat*GC_M[startvalue:endvalue, :, None]).sum(axis=1)
					GC_M_on_sat_subset = GC_M_on_sat_subset/GC_to_sat.sum(axis=1)
					GC_M_on_sat[startvalue:endvalue,:] = GC_M_on_sat_subset
		else:
			for i in range(np.shape(GC_edges)[0]):
				# Define vectors that give the "low" and "high" pressure
				# values for each GEOS-Chem and satellite layer.
				GC_lo = GC_edges[i, 1:][:, None]
				GC_hi = GC_edges[i, :-1][:, None]
				sat_lo = sat_edges[i, 1:][None, :]
				sat_hi = sat_edges[i, :-1][None, :]
				idx = (np.less_equal(sat_lo, GC_hi) & np.greater_equal(sat_hi, GC_lo))
				GC_to_sat = np.minimum(sat_hi, GC_hi) - np.maximum(sat_lo, GC_lo)
				GC_to_sat[~idx] = 0
				# Now map the GC CH4 to the satellite levels
				GC_on_sat_subset = (GC_to_sat*GC_SPC[i, :, None]).sum(axis=0)
				GC_on_sat_subset = GC_on_sat_subset/GC_to_sat.sum(axis=0)
				GC_on_sat[i,:] = GC_on_sat_subset
				if GC_M is not None:
					GC_M_on_sat_subset = (GC_to_sat*GC_M[i, :, None]).sum(axis=0)
					GC_M_on_sat_subset = GC_M_on_sat_subset/GC_to_sat.sum(axis=0)
					GC_M_on_sat[i,:] = GC_M_on_sat_subset
	else:
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
		if GC_M is not None:
			GC_M_on_sat = (GC_to_sat*GC_M[:, :, None]).sum(axis=1)
			GC_M_on_sat = GC_M_on_sat/GC_to_sat.sum(axis=1)
	if GC_M is not None:
		return [GC_on_sat,GC_M_on_sat]
	else:
		return GC_on_sat

def apply_avker(sat_avker, sat_pressure_weight, GC_SPC, sat_prior=None,GC_M_on_sat=None,GC_area=None,filt=None):
	'''
	Apply the averaging kernel
	Inputs:
		sat_avker			The averaging kernel for the satellite
		sat_prior			The satellite prior profile in ppb, optional (used for CH4)
		sat_pressure_weight  The relative pressure weights for each level
		GC_SPC			   The GC species on the satellite levels
		filt				 A filter, optional
	'''
	if filt is None:
		filt = np.ones(sat_avker.shape[1])
	else:
		filt = filt.astype(int)
	if GC_M_on_sat is not None: #Take partial columns, which also involves the area and air mass
		GC_SPC = GC_SPC/0.02897 #convert to mol/kg air
		GC_SPC = GC_SPC*GC_M_on_sat #Convert to mol of interest per box
		GC_area = GC_area.reshape((len(GC_area),1)) #reshape to make conformable
		GC_SPC = GC_SPC/GC_area #Convert to mol of interest per box per m2, which is TROPOMI dimensions
	if sat_prior is None:
		GC_col = (filt*sat_pressure_weight*sat_avker*GC_SPC)
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
		if specieskey in list(self.spc_config["filter_obs_poleward_of_n_degrees"].keys()):
			filterinfo['MAIN']=[float(self.spc_config["filter_obs_poleward_of_n_degrees"][specieskey])]
		for obs in obs_list:
			trop_obs.append(read_tropomi(obs,species,filterinfo,includeObsError=includeObsError))
		met = {}
		for key in list(trop_obs[0].keys()):
			met[key] = np.concatenate([metval[key] for metval in trop_obs])
		return met
	def gcCompare(self,specieskey,TROPOMI,GC,GC_area=None,saveAlbedo=False,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		if species=='CH4':
			TROP_PRIOR = 1e9*(TROPOMI['methane_profile_apriori']/TROPOMI['dry_air_subcolumns'])
			synthetic_partial_columns = False
		elif species=='NO2':
			TROP_PRIOR=None
			synthetic_partial_columns = True
		TROP_PW = (-np.diff(TROPOMI['pressures'])/(TROPOMI['pressures'][:, 0] - TROPOMI['pressures'][:, -1])[:, None])
		returnStateMet = self.spc_config['SaveStateMet']=='True'
		if returnStateMet:
			GC_SPC,GC_P,GC_M,GC_area,i,j,t = obsop.getGCCols(GC,TROPOMI,species,returninds=True,returnStateMet=returnStateMet,GC_area=GC_area)
		else:
			GC_SPC,GC_P,GC_area,i,j,t = obsop.getGCCols(GC,TROPOMI,species,returninds=True,returnStateMet=returnStateMet,GC_area=GC_area)
		if species=='CH4':
			GC_SPC*=1e9 #scale to mol/mol
		#TROPOMI_ALL setting extension must be on!!
		memsetting = self.spc_config['LOW_MEMORY_TROPOMI_AVERAGING_KERNEL_CALC'] == 'True'
		if memsetting:
			batchsize = int(self.spc_config['LOW_MEMORY_TROPOMI_AVERAGING_KERNEL_BATCH_SIZE'])
		else:
			batchsize = None
		if returnStateMet:
			GC_on_sat,GC_M_on_sat = GC_to_sat_levels(GC_SPC, GC_P, TROPOMI['pressures'],GC_M = GC_M, lowmem=memsetting,batchsize=batchsize)
		else:
			GC_on_sat = GC_to_sat_levels(GC_SPC, GC_P, TROPOMI['pressures'],lowmem=memsetting,batchsize=batchsize)
			GC_M_on_sat = None
		GC_on_sat = apply_avker(TROPOMI['column_AK'],TROP_PW, GC_on_sat,TROP_PRIOR,GC_M_on_sat,GC_area)
		if self.spc_config['AV_TO_GC_GRID']=="True":
			superObsFunction = self.spc_config['SUPER_OBSERVATION_FUNCTION'][specieskey]
			additional_args_avgGC = {}
			if useObserverError:
				additional_args_avgGC['obsInstrumentError'] = TROPOMI['Error']
				additional_args_avgGC['modelTransportError'] = transportError
			elif prescribed_error is not None:
				additional_args_avgGC['prescribed_error'] = prescribed_error
				additional_args_avgGC['prescribed_error_type'] = prescribed_error_type
			if saveAlbedo:
				additional_args_avgGC['albedo_swir'] = TROPOMI['albedo_swir']
				additional_args_avgGC['albedo_nir'] = TROPOMI['albedo_nir']
				additional_args_avgGC['blended_albedo'] = TROPOMI['blended_albedo']
			if minError is not None:
				additional_args_avgGC['minError'] = minError
			if errorCorr is not None:
				additional_args_avgGC['errorCorr'] = errorCorr
			toreturn = obsop.averageByGC(i,j,t,GC,GC_on_sat,TROPOMI[species],doSuperObs=True,superObsFunction=superObsFunction,**additional_args_avgGC)
		else:
			toreturn = obsop.ObsData(GC_on_sat,TROPOMI[species],TROPOMI['latitude'],TROPOMI['longitude'],TROPOMI['utctime'])
			if saveAlbedo:
				toreturn.addData(swir_av=TROPOMI['albedo_swir'],nir_av=TROPOMI['albedo_nir'],blended_av=TROPOMI['blended_albedo'])
			if useObserverError:
				toreturn.addData(err_av=TROPOMI['Error'])
		return toreturn


