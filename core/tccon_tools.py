# TCCON CO operator (v1) prepared by Sina Voshtani -- with pressure adjustment (2024/03/15)
# Update: Include TCCON InSb CO measurements at East Trout Lake. Bug fixed and tested (2024/08/30)

## TCCON operator incorporated into main CHEEREIO code branch and extended to N2O by Drew Pendergrass (2025/07/24)
## Update: included TCCON N2O temperature correction from Josh Laughner; will be unnecessary after GGG2020.1 (2025/08/05)

from datetime import datetime
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np
import observation_operators as obsop
from scipy.interpolate import interp1d


def read_tccon(filename, species, filterinfo=None, includeObsError = False,doN2OCorrectionPT700=False):

	met = {}
	
	data = xr.open_dataset(filename)
	if species=='CO':
		met[species] = data['xco'].values # TCCON column (ppb). Dim: site
	elif species=='N2O':
		met[species] = data['xn2o'].values # TCCON column (ppb)
	else:
		raise ValueError('Species not supported')

	if includeObsError:
		if species=='CO':
			met['Error'] = data['xco_error'].values #TCCON column error (ppb)
		elif species=='N2O':
			met['Error'] = data['xn2o_error'].values #TCCON column error (ppb)
	
	met['longitude'] = data['long'].values
	met['latitude'] = data['lat'].values
	met['utctime'] = data['time_utc'].values
	met['station_id'] = data['station_id'].values

	met['pressure_AK'] = data['ak_pressure'].values # hPa 
	met['pressure_apriori'] = data['prior_pressure'].values * 101325/100 # hPa  
	met['h2o_profile_apriori'] = data['prior_h2o'].values # parts
	met['altitude_apriori'] = 1000*data['prior_altitude'].values # m
	met['pout'] = data['pout'].values # TCCON (surface) pressure hPa

	#Do N2O temp correction, if we are using GGG2020 (this should be included in GGG2020.1 when it's released)
	if doN2OCorrectionPT700 and species=='N2O':
		if includeObsError:
			met[species],met['Error']=correct_xn2o_from_pt700(met[species],data['prior_temperature'].values,met['pressure_apriori'],xn2o_error=met['Error'])
		else:
			met[species]=correct_xn2o_from_pt700(met[species],data['prior_temperature'].values,met['pressure_apriori'])

	if species=='CO':
		met['column_AK'] = data['ak_xco'].values #time,latitude,longitude,layer
		met['co_profile_apriori'] = data['prior_co'].values # ppb
	elif species=='N2O':
		met['column_AK'] = data['ak_xn2o'].values #time,latitude,longitude,layer
		met['n2o_profile_apriori'] = data['prior_n2o'].values # ppb
	
	if filterinfo is not None:
		met = obsop.apply_filters(met,filterinfo)
	
	data.close()
	return met


# Based on Josh Laughner's code from GitHub: py_tccon_netcdf/write_tccon_netcdf/bias_corrections.py
def correct_xn2o_from_pt700(xn2o,prior_temperature,prior_pressure,xn2o_error=None,n2o_aicf=0.9821, m=0.000626, b=0.787):
	# We need to remove the AICF because the correction was calculated for pre-AICF XN2O data
	xn2o = xn2o * n2o_aicf
	if xn2o_error is not None:
		xn2o_error = xn2o_error * n2o_aicf
	pt700 = _compute_pt700(prior_temperature,prior_pressure)
	# Apply the temperature correction. Counterintuitively, we do *NOT* need to reapply the AICF.
	# That is only because for N2O the AICF was calculated as the value of this fit at 310 K.
	# Hence, essentially this is applying a *temperature dependent* AICF. 
	xn2o_corr = xn2o / (m * pt700 + b)
	if xn2o_error is not None:
		xn2o_error_corr = xn2o_error / (m * pt700 + b)
		return xn2o_corr, xn2o_error_corr
	else:
		return xn2o_corr

# Based on Josh Laughner's code from GitHub: py_tccon_netcdf/write_tccon_netcdf/bias_corrections.py
def _compute_pt700(prior_temperature,prior_pressure):
	pt700 = np.full(prior_temperature.shape[0], np.nan, dtype=prior_temperature.dtype)
	for (i, (tprof, pprof)) in enumerate(zip(prior_temperature, prior_pressure)):
		# Normally I would interpolate temperature vs. ln(p). But the bias plot does a 
		# straight linear-linear interpolation, so we do the same so that our PT that
		# we use to determine the bias correction is calculated the same way as that
		# used to determine the bias correction slope.
		f = interp1d(pprof, tprof)
		t700_i = f(700.0)
		pt700[i] = t700_i * (1000.0 / 700.0) ** 0.286
	return pt700

# map the vertical levels of the model into the observation levels
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


# compute g at vertical layers for each observation site
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
			TCCON_date_dict[key]['start'] = [datetime.strptime(obs.split('_')[-2], "%Y%m%dT%H%M%S") for obs in obs_list]
			TCCON_date_dict[key]['end'] = [datetime.strptime(obs.split('_')[-1].split('.')[0], "%Y%m%dT%H%M%S") for obs in obs_list]
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
		pt700=False
		if species=='CO':
			if (self.spc_config['Extensions']['TCCON_CO']=="True") and (self.spc_config['TCCON_CO_FILTERS']=="True"):
				pass #no filters implemented
		elif species=='N2O':
			if (self.spc_config['Extensions']['TCCON_N2O']=="True") and (self.spc_config['TCCON_N2O_FILTERS']=="True"):
				pass #no filters implemented
			if (self.spc_config['Extensions']['TCCON_N2O']=="True") and (self.spc_config['PT700_N2O_CORRECTION']=="True"):
				pt700=True #Do GGG2020 N2O correction
		if specieskey in list(self.spc_config["filter_obs_poleward_of_n_degrees"].keys()):
			filterinfo['MAIN']=[float(self.spc_config["filter_obs_poleward_of_n_degrees"][specieskey])]
		# read tccon files
		for obs in obs_list:
			tccon_obs.append(read_tccon(obs,species,filterinfo,includeObsError=includeObsError,doN2OCorrectionPT700=pt700))
		met = {}
		if len(tccon_obs)>0:
			for key in list(tccon_obs[0].keys()):
				met[key] = np.concatenate([metval[key] for metval in tccon_obs])
		else: #If no observations, return something empty
			met[species] = np.array([])
			if includeObsError:
				met['Error'] = np.array([])
			met['longitude'] = np.array([])
			met['latitude'] = np.array([])
			met['station_id'] = np.array([])
			met['utctime'] = np.array([])
			met['pressure_AK'] = np.array([])
			met['pressure_apriori'] = np.array([])
			met['h2o_profile_apriori'] = np.array([])
			met['altitude_apriori'] = np.array([])
			met['pout'] = np.array([])
			met['column_AK'] = np.array([])
			if species=='CO':
				met['co_profile_apriori'] = np.array([])
			elif species=='N2O':
				met['n2o_profile_apriori'] = np.array([])
		return met
	def gcCompare(self,specieskey,TCCON,GC,GC_area=None,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		extra_obsdata_to_save = self.spc_config['EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y'][specieskey]
		returnStateMet = self.spc_config['SaveStateMet']=='True'
		if len(TCCON[species])==0: #return empty arrays
			toreturn = obsop.ObsData(np.array([]),TCCON[species],TCCON['latitude'],TCCON['longitude'],TCCON['utctime'])
			#If saving extra fields, add them here
			if len(extra_obsdata_to_save)>0:
				data_to_add = {}
				for field in extra_obsdata_to_save:
					data_to_add[field] = TCCON[field]
				toreturn.addData(**data_to_add)
			if doErrCalc and useObserverError:
				toreturn.addData(err_av=TCCON['Error'])
		else:
			GC_col_data = obsop.getGCCols(GC,TCCON,species,self.spc_config,returninds=True,returnStateMet=returnStateMet,GC_area=GC_area)
			GC_SPC = GC_col_data['GC_SPC']
			GC_P = GC_col_data['GC_P'] # (PEDGE)
			GC_H2O = GC_col_data['Met_AVGW'] #We require this field from Met History for TCCON observations
			i,j,t = GC_col_data['indices']
			ff = TCCON['pout']/TCCON['pressure_apriori'][:, 0]
			TCCON_P_fix = TCCON['pressure_apriori'] * ff[:, np.newaxis]	
			GC_on_sat_l  = GC_to_sat_levels(GC_SPC, GC_P, TCCON_P_fix) # GC a-priori co on TCCON layer
			GC_on_sat_l_h2o = GC_to_sat_levels(GC_H2O, GC_P, TCCON_P_fix) # GC a-priori h2o on TCCON layer
			if species=="CO":
				GC_on_sat = integrate_column(GC_on_sat_l,GC_on_sat_l_h2o,TCCON['h2o_profile_apriori'],TCCON['pout'],TCCON['pressure_apriori'],TCCON['altitude_apriori'][:51],TCCON['co_profile_apriori'],TCCON['latitude'],TCCON['column_AK'])
			elif species=="N2O":
				GC_on_sat = integrate_column(GC_on_sat_l,GC_on_sat_l_h2o,TCCON['h2o_profile_apriori'],TCCON['pout'],TCCON['pressure_apriori'],TCCON['altitude_apriori'][:51],TCCON['n2o_profile_apriori'],TCCON['latitude'],TCCON['column_AK'])
			else: 
				raise ValueError(f'Species {species} not recognized')
			#print("GC-TCCON:", GC_on_sat - TCCON[species])
			nan_indices = np.argwhere(np.isnan(GC_on_sat))
			GC_on_sat = np.nan_to_num(GC_on_sat)
			if self.spc_config['AV_TO_GC_GRID'][specieskey]=="True":
				superObsFunction = self.spc_config['SUPER_OBSERVATION_FUNCTION'][specieskey]
				additional_args_avgGC = {}
				if doErrCalc:
					if useObserverError:
						TCCON['Error']=1*TCCON['Error'] # it is in ppb already
						additional_args_avgGC['obsInstrumentError'] = TCCON['Error']
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

