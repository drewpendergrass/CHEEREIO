# TCCON CO operator (v1) prepared by Sina Voshtani -- with pressure adjustment (2024/03/15)
# Update: Include TCCON InSb CO measurements at East Trout Lake. Bug fixed and tested (2024/08/30)

from datetime import datetime
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np
import observation_operators as obsop
from scipy.interpolate import interp1d

def read_tccon(filename, species, filterinfo=None, includeObsError = False):

	met = {}
	
	data = xr.open_dataset(filename)
	idx = data['xco'].values[0,:,:] #time,latitude,longitude
	if species=='CO':
		la,lo=np.where(idx>0) # include all quality controlled XCO
		met[species] = data['xco'].values[0,la,lo] # TCCON column (ppb)
	else:
		raise ValueError('Species not supported')

	if includeObsError:
		if species=='CO':
			met['Error'] = data['xco_error'].values[0,la,lo] #TCCON column error (ppb)
	
	met['longitude'] = data['longitude'].values[:]
	met['latitude'] = data['latitude'].values[:]
	met['utctime'] = data['time_utc'].values[:]

	if species=='CO':
		met['column_AK'] = data['ak_xco'].values[0,la,lo,:] #time,latitude,longitude,layer
		met['pressure_AK'] = data['ak_pressure'].values[0,la,lo,:] # hPa 
		met['pressure_apriori'] = data['prior_pressure'].values[0,la,lo,:] * 101325/100 # hPa  
		met['co_profile_apriori'] = data['prior_co'].values[0,la,lo,:] # ppb
		met['h2o_profile_apriori'] = data['prior_h2o'].values[0,la,lo,:] # parts
		met['altitude_apriori'] = 1000*data['prior_altitude'].values[:] # m
		met['pout'] = data['pout'].values[0,la,lo] # TCCON (surface) pressure hPa
		data.close()
	
	if filterinfo is not None:
		met = obsop.apply_filters(met,filterinfo)

	return met

# map the veritcal levels of the model into the observation levels
def GC_to_sat_levels(GC_SPC, GC_edges, sat_edges):
	# Adjusting bottom and top edges
	idx_bottom = np.less(GC_edges[:, 0], sat_edges[:, 0])
	idx_top = np.greater(GC_edges[:, -1], sat_edges[:, -1])
	GC_edges[idx_bottom, 0] = sat_edges[idx_bottom, 0]
	GC_edges[idx_top, -1] = sat_edges[idx_top, -1]
	
	# Defining vectors that give the "low" and "high" pressure values for each GEOS-Chem and satellite layer.
	GC_lo = GC_edges[:, 1:][:, :, None]
	GC_hi = GC_edges[:, :-1][:, :, None]
	sat_lo = sat_edges[:, 1:][:, None, :]
	sat_hi = sat_edges[:, :-1][:, None, :]
	
	# Get the indices where the GC-to-satellite mapping is non-zero
	idx = (np.less_equal(sat_lo, GC_hi) & np.greater_equal(sat_hi, GC_lo))
	
	# Find the fraction of each GC level that contributes to each TCCON level
	GC_to_sat = np.minimum(sat_hi, GC_hi) - np.maximum(sat_lo, GC_lo)
	GC_to_sat[~idx] = 0
	
	# Now map the GC_SPC to the satellite levels
	GC_on_sat = (GC_to_sat * GC_SPC[:, :, None]).sum(axis=1)
	GC_on_sat = GC_on_sat / GC_to_sat.sum(axis=1)

	if np.any(np.isnan(GC_on_sat)) or np.any(np.isinf(GC_on_sat)):
	    print("Warning: NaN or inf values encountered in GC_on_sat.")

	return GC_on_sat


# compute g at verticl layers for each observation site
def gravity(altitudes, latitudes):
	gm = 3.9862216e+14  # Gravitational constant times Earth's Mass (m3/s2) in GFIT
	omega = 7.292116E-05  # Earth's angular rotational velocity (radians/s) from GFIT
	con = 0.006738  # (a/b)**2-1 where a & b are equatorial & polar radii from GFIT
	shc = 1.6235e-03  # 2nd harmonic coefficient of Earth's gravity field from GFIT
	eqrad = 6378178  # Equatorial Radius (m) from GFIT
	
	# Ensure altitudes and latitudes are numpy arrays
	altitudes = np.asarray(altitudes)
	latitudes = np.asarray(latitudes)
	
	# Initialize an array to store the gravitational accelerations
	g_array = np.zeros((len(altitudes) - 1, len(latitudes)))
	
	for i, lat in enumerate(latitudes):
		lats = np.full_like(altitudes, lat)
		gclat = np.arctan(np.tan(lats * np.pi / 180) / (1 + con))
		radius = altitudes + eqrad / (1 + con * np.sin(gclat) ** 2) ** 0.5
		ff = (radius / eqrad) ** 2
		hh = radius * omega ** 2
		ge = gm / (eqrad ** 2)
		
		g = (ge * (1 - shc * (3 * np.sin(gclat) ** 2 - 1) / ff) / ff - hh * np.cos(gclat) ** 2) * (
		    1 + 0.5 * (np.sin(gclat) * np.cos(gclat) * (hh / ge + 2 * shc / ff ** 2)) ** 2)
		
		# Average every two levels for the middle values
		g_layer = 0.5 * (g[:-1] + g[1:])
		g_array[:, i] = g_layer
		
	return g_array

# this is the main function for TCCON column intergation, using its a-priori and averaging kernels
def integrate_column(gas_profile,h2o_profile,obh2o_profile,obpout,obpressure_profile,altitude_profile,ensemble_profile,oblat,AK):
	# Inputs:
	# gas_profile: the model gas profile of interest in ppm
	# h2o_profile: the model h2o profile in ppm
	# obh2o_profile: the observation h2o profile in ppm
	# pressure_profile: the model pressure profile that corresponds with the
	# gas profile in hPa
	# ensemble_profile: the ensemble profile (xc) from equation 25 in Rodgers
	# and Connor 2000 - this will likely be the a priori profile from GFIT; most often 
	# multiplied by the scaling factors (VSF) from GFIT for the spectra near the aircraft overpass
	# obalt: the altitude of the ground-based site in m - geometric altitude
	# oblat: the latitude (in degrees) of the ground-based site
	# AK: the averaging kernels for all windows of the molecule of interest, in a structure
	# AK_P: the averaging kernel pressure levels	
	# Necessary constants
	Na = 6.0221415e23 # molecules/mol
	m_dry_air = 28.9644*1e-3/Na # in kg/molecule
	m_h2o = 18.01534*1e-3/Na # in kg/molecule
	
	g_array = gravity(altitude_profile,oblat) 
	g_layer = g_array.T # mean gravitational acceleration for each layer	
	
	##### The model integration #####	
	f_h2o = h2o_profile # h2o_profile is already in parts (if in ppm, use: f_h2o = h2o_profile*1e-6)
	f_dry_h2o = f_h2o # h2o is already dry as taken from Met_AVGW (to dry wet mole fraction for h2o, use: f_dry_h2o = f_h2o/(1-f_h2o))
	
	dry_gas_profile = gas_profile # it is dry already (to dry wet gas profile (vmr), use: dry_gas_profile = gas_profile/(1-f_dry_h2o))
	f_dry_gas = dry_gas_profile # gas_profile is already in parts (if in ppb, use: f_dry_gas = dry_gas_profile*1e-9) 

	rr = obpout/obpressure_profile[:, 0] # correct observation pressure profile using surface pressure adjustement
	obpressure_profile_fix = obpressure_profile * rr[:, np.newaxis]
	#obpressure_profile_fix = obpressure_profile # use only if pressure profile is coorect 
	
	dP = 100 * (np.diff(obpressure_profile_fix, axis=1)) # layer thickness for the integration		
	f_dry_gas_layer = f_dry_gas # mean dry mol fraction co for each layer	
	f_dry_h2o_layer = f_dry_h2o # mean dry mol fraction h2o for each layer
	
	VC_gas_integrand = f_dry_gas_layer/(g_layer*m_dry_air*(1+f_dry_h2o_layer*(m_h2o/m_dry_air))) #  integration over layers of gas 
	VC_gas = np.nansum(VC_gas_integrand*abs(dP), axis=1) # vertical column of the gas in molecules/m^2
	
	VC_air_integrand = 1/(g_layer*m_dry_air*(1+f_dry_h2o_layer*(m_h2o/m_dry_air))) # integration over layers of air
	VC_air = np.nansum(VC_air_integrand*abs(dP), axis=1) # vertical column of dry air in molecules/m^2
	
	gas = (VC_gas/VC_air)*10**9 # dry mole fraction of the gas - column   
	
	AK_gas = AK
	AK_gas_layer = 0.5 * (AK_gas[:, :-1] + AK_gas[:, 1:])
	
	##### The ensemble integration #####
	
	f_obh2o = obh2o_profile # no coversion needed as TCCON water profile is already in "parts". vmr for h2o
	f_dry_obh2o = f_obh2o/(1-f_obh2o) # to dry mol fraction for h2o (TCCON prior water profile is wet)
	
	dry_ensemble_profile = ensemble_profile/(1-f_dry_obh2o) # convert a priori vmr into dry vmr
	#dry_ensemble_profile = ensemble_profile # in OSSE it is dry
 
	dry_ensemble_profile_layer = 0.5 * (dry_ensemble_profile[:, :-1] + dry_ensemble_profile[:, 1:])
	f_dry_gas_ensemble = dry_ensemble_profile*1e-9 # ppb to parts for e.g.,  TCCN['co_profile_apriori']
	
	f_dry_gas_ensemble_fix = np.empty_like(f_dry_gas_ensemble) # adjusting gas profile based on the surface pressure correction (obpressure_profile_fix)
	for i in range(f_dry_gas_ensemble.shape[0]):	
		interp = interp1d(obpressure_profile[i, :], f_dry_gas_ensemble[i, :], kind='linear', fill_value="extrapolate", axis=0)
		f_dry_gas_ensemble_fix[i, :] = interp(obpressure_profile_fix[i, :])

	f_dry_gas_layer_ensemble = 0.5 * (f_dry_gas_ensemble_fix[:, :-1] + f_dry_gas_ensemble_fix[:, 1:]) # mean dry mol fraction co for each layer
	
	f_dry_obh2o_fix = np.empty_like(f_dry_obh2o) # adjusting h2o profile based on the surface pressure correction (obpressure_profile_fix)
	for i in range(f_dry_obh2o.shape[0]):	
		interp = interp1d(obpressure_profile[i, :], f_dry_obh2o[i, :], kind='linear', fill_value="extrapolate", axis=0)
		f_dry_obh2o_fix[i, :] = interp(obpressure_profile_fix[i, :])
	f_dry_obh2o_layer = 0.5 * (f_dry_obh2o_fix[:, :-1] + f_dry_obh2o_fix[:, 1:]) # mean dry mol fraction h2o for each layer
	
	VC_gas_integrand_ensemble = f_dry_gas_layer_ensemble/(g_layer*m_dry_air*(1+f_dry_obh2o_layer*(m_h2o/m_dry_air))) # integration over layers of gas
	VC_gas_ensemble = np.nansum(VC_gas_integrand_ensemble*abs(dP), axis=1)  # vertical column of the gas in molecules/m^2
	
	# The ensemble with averaging kernel (for testing purposes only):
	f_dry_gas_ensemble_ak = AK_gas*f_dry_gas_ensemble_fix
	f_dry_gas_layer_ensemble_ak = 0.5 * (f_dry_gas_ensemble_ak[:, :-1] + f_dry_gas_ensemble_ak[:, 1:]) # mean dry mol fraction co for each layer	
	VC_gas_integrand_ensemble_ak = f_dry_gas_layer_ensemble_ak/(g_layer*m_dry_air*(1+f_dry_h2o_layer*(m_h2o/m_dry_air)))
	VC_gas_ensemble_ak = np.nansum(VC_gas_integrand_ensemble_ak*abs(dP), axis=1)
	
	# The true profile with averaging kernel (for testing purposes only):
	f_dry_gas_fix = f_dry_gas ## GC is given in parts
	f_dry_true_gas_ak = AK_gas_layer*f_dry_gas_fix
	f_dry_true_gas_layer_ak = f_dry_true_gas_ak 
	VC_true_gas_integrand_ak = f_dry_true_gas_layer_ak/(g_layer*m_dry_air*(1+f_dry_h2o_layer*(m_h2o/m_dry_air)))
	VC_true_gas_ak = np.nansum(VC_true_gas_integrand_ak*abs(dP), axis=1)
		
	# The difference between the true profile and the ensemble without ak:
	f_dry_gas_diff = f_dry_gas_fix - f_dry_gas_layer_ensemble  
	f_dry_gas_layer_diff = f_dry_gas_diff
	VC_gas_integrand_diff = f_dry_gas_layer_diff/(g_layer*m_dry_air*(1+f_dry_h2o_layer*(m_h2o/m_dry_air)))
	VC_gas_diff = np.nansum(VC_gas_integrand_diff*abs(dP), axis=1)
	
	# The difference between the true profile and the ensemble with ak:
	f_dry_gas_ak_diff = AK_gas_layer*(f_dry_gas_fix - f_dry_gas_layer_ensemble)
	f_dry_gas_layer_ak_diff = f_dry_gas_ak_diff
	VC_gas_integrand_ak_diff = f_dry_gas_layer_ak_diff/(g_layer*m_dry_air*(1+f_dry_h2o_layer*(m_h2o/m_dry_air)))
	VC_gas_ak_diff = np.nansum(VC_gas_integrand_ak_diff*abs(dP), axis=1)

	# column mole fractions in ppb	
	apriori = (VC_gas_ensemble/VC_air)*10**9
	true_gas = (VC_gas_ensemble +  VC_gas_diff)/VC_air*10**9
	true_gas_ak = (VC_gas_ensemble + VC_gas_ak_diff)/VC_air*10**9
	### return apriori, true_gas, true_gas_ak # for verification purpose 
	return true_gas_ak

class TCCON_Translator(obsop.Observation_Translator):
	def __init__(self,verbose=1):
		super().__init__(verbose)
	#Save dictionary of dates for later use
	def initialReadDate(self):
		sourcedirs = self.spc_config['TCCON_dirs']
		TCCON_date_dict = {}
		for key in list(sourcedirs.keys()):
			sourcedir = sourcedirs[key]
			obs_list = glob(f'{sourcedir}/**/tccon_avg_*.nc', recursive=True)
			obs_list.sort()
			TCCON_date_dict[key] = {}
			TCCON_date_dict[key]['start'] = [datetime.strptime(obs.split('_')[-3], "%Y%m%dT%H%M%S") for obs in obs_list]
			TCCON_date_dict[key]['end'] = [datetime.strptime(obs.split('_')[-2], "%Y%m%dT%H%M%S") for obs in obs_list]
		with open(f"{self.scratch}/tccon_dates.pickle", 'wb') as handle:
			pickle.dump(TCCON_date_dict, handle)
		return TCCON_date_dict
	#Timeperiod is two datetime objects
	def globObs(self,species,timeperiod, interval=None):
		sourcedir = self.spc_config['TCCON_dirs'][species]
		if os.path.exists(f"{self.scratch}/tccon_dates.pickle"):
			with open(f"{self.scratch}/tccon_dates.pickle", 'rb') as handle:
				TCCON_date_dict = pickle.load(handle)
		else:
			TCCON_date_dict = self.initialReadDate()
		obs_dates = TCCON_date_dict[species]
		obs_list = glob(f'{sourcedir}/**/tccon_avg_*.nc', recursive=True)
		obs_list.sort()

		if interval:
			obs_list = [obs for obs,t1,t2 in zip(obs_list,obs_dates['start'],obs_dates['end']) if (t1>=timeperiod[0]) and (t2<timeperiod[1]) and ((t1.hour % interval == 0) or (t1.hour % interval == (interval-1)))]
		else:
			obs_list = [obs for obs,t1,t2 in zip(obs_list,obs_dates['start'],obs_dates['end']) if (t1>=timeperiod[0]) and (t2<timeperiod[1])]
			print("obs_list:", obs_list)		
		return obs_list
	def getObservations(self,specieskey,timeperiod, interval=None, includeObsError=True):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		obs_list = self.globObs(species,timeperiod,interval)
		tccon_obs = []
		filterinfo = {}
		if species=='CO':
			if (self.spc_config['Extensions']['TCCON_CO']=="True") and (self.spc_config['TCCON_CO_FILTERS']=="True"):
				filterinfo["TCCON_CO"] = [float(self.spc_config['TCCON_CO_filter_blended_albedo']),float(self.spc_config['TCCON_CO_filter_swir_albedo_low']),float(self.spc_config['TCCON_CO_filter_swir_albedo_high']),float(self.spc_config['TCCON_CO_filter_winter_lat']),float(self.spc_config['TCCON_CO_filter_roughness']),float(self.spc_config['TCCON_CO_filter_swir_aot'])]
		if specieskey in list(self.spc_config["filter_obs_poleward_of_n_degrees"].keys()):
			filterinfo['MAIN']=[float(self.spc_config["filter_obs_poleward_of_n_degrees"][specieskey])]
		
		# read tccon files
		for obs in obs_list:
			if self.spc_config['WHICH_TROPOMI_PRODUCT'] == 'ACMG':
				trop_obs.append(read_tropomi_acmg(obs,species,filterinfo,includeObsError=includeObsError))
			else: 
				tccon_obs.append(read_tccon(obs,species,filterinfo,includeObsError=includeObsError))
		met = {}
		for key in list(tccon_obs[0].keys()):
			met[key] = np.concatenate([metval[key] for metval in tccon_obs])
		return met
	def gcCompare(self,specieskey,TCCON,GC,GC_area=None,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		extra_obsdata_to_save = self.spc_config['EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y'][specieskey]
		returnStateMet = self.spc_config['SaveStateMet']=='True'
		GC_col_data = obsop.getGCCols(GC,TCCON,species,self.spc_config,returninds=True,returnStateMet=returnStateMet,GC_area=GC_area)
		GC_SPC = GC_col_data['GC_SPC']
		GC_P = GC_col_data['GC_P'] # (PEDGE)
		GC_H2O = GC_col_data['GC_H2O']
		i,j,t = GC_col_data['indices']
		if species=='CO':
			synthetic_partial_columns = False
		ff = TCCON['pout']/TCCON['pressure_apriori'][:, 0]
		TCCON_P_fix = TCCON['pressure_apriori'] * ff[:, np.newaxis]	
		GC_on_sat_l_co  = GC_to_sat_levels(GC_SPC, GC_P, TCCON_P_fix) # GC a-priori co on TCCON layer
		GC_on_sat_l_h2o = GC_to_sat_levels(GC_H2O, GC_P, TCCON_P_fix) # GC a-priori h2o on TCCON layer
		GC_on_sat = integrate_column(GC_on_sat_l_co,GC_on_sat_l_h2o,TCCON['h2o_profile_apriori'],TCCON['pout'],TCCON['pressure_apriori'],TCCON['altitude_apriori'][:51],TCCON['co_profile_apriori'],TCCON['latitude'],TCCON['column_AK'])
		#print("GC-TCCON:", GC_on_sat - TCCON[species])
		nan_indices = np.argwhere(np.isnan(GC_on_sat))
		GC_on_sat = np.nan_to_num(GC_on_sat)
		if self.spc_config['AV_TO_GC_GRID'][specieskey]=="True":
			superObsFunction = self.spc_config['SUPER_OBSERVATION_FUNCTION'][specieskey]
			additional_args_avgGC = {}
			if doErrCalc:
				if useObserverError:
					if species=='CO':
						TCCON['Error']=1*TCCON['Error'] # it is in ppb already
					additional_args_avgGC['obsInstrumentError'] = TCCON['Error']
					additional_args_avgGC['modelTransportError'] = transportError
				elif prescribed_error is not None:
					additional_args_avgGC['prescribed_error'] = prescribed_error
					additional_args_avgGC['prescribed_error_type'] = prescribed_error_type
				if minError is not None:
					additional_args_avgGC['minError'] = minError
				if errorCorr is not None:
					additional_args_avgGC['errorCorr'] = errorCorr
			#If saving extra fields, add them here
			if len(extra_obsdata_to_save)>0:
				additional_args_avgGC['other_fields_to_avg'] = {}
				for field in extra_obsdata_to_save:
					additional_args_avgGC['other_fields_to_avg'][field] = TCCON[field]
			toreturn = obsop.averageByGC(i,j,t,GC,GC_on_sat,TCCON[species],doSuperObs=doErrCalc,superObsFunction=superObsFunction,**additional_args_avgGC)
		else:
			timevals = GC.time.values[t]
			toreturn = obsop.ObsData(GC_on_sat,TCCON[species],TCCON['latitude'],TCCON['longitude'],timevals)
			#If saving extra fields, add them here
			if len(extra_obsdata_to_save)>0:
				data_to_add = {}
				for field in extra_obsdata_to_save:
					data_to_add[field] = TCCON[field]
				toreturn.addData(**data_to_add)
			if doErrCalc and useObserverError:
				toreturn.addData(err_av=TCCON['Error'])
		return toreturn

