#Tools for reading IASI NH3

from datetime import datetime
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np
import observation_operators as obsop

def read_iasi(filename, species, filterinfo=None, includeObsError = False):
	"""
	Read IASI data and save important variables to dictionary.

	Arguments
		filename [str]  : IASI netcdf data file to read
		species [str]   : Species string (currently only NH3 supplorted)
		filterinfo [list] : NOT CURRENTLY SUPPORTED; leave as None
		includeObsError [bool] : if True, combine systematic and random uncertainty in one number

	Returns
		met	  [dict] : Dictionary of important variables from TROPOMI:
							- species
 							- Latitude
							- Longitude
	   						- QA value
 							- UTC time
							- Averaging kernel
	"""

	if (filterinfo is not None) or (species != 'NH3'):
		raise ValueError('Only supports NH3 with no filter info at this time.')

	# Initialize list for TROPOMI data
	met = {}
	
	# Store species, QA, lat, lon, time, averaging kernel
	data = xr.open_dataset(filename)
	qa_init = data['prefilter'].values
	ind = np.where(qa==1)[0] #Keep only those passing prefilter; will apply postfilter later

	met['qa_value'] = qa[ind]
	met['utc_time'] = data['AERIStime'].values[ind]
	met['lat'] = data['latitude'].values[ind]
	met['lon'] = data['longitude'].values[ind]
	met['level_edge'] = data['levels'].values #in km, same for all obs (length 15)
	met['level_middle'] = data['midlevels'].values #in km, same for all obs (length 14)
	met[species] = data['nh3_total_column'].values[ind]

	#We are applying Method 2 from the avkReadMe (accompanying Clarisse et al 2023). 
	met['column_AK'] = (data['nh3_AvKnorm']**-1)*(met[species]/data['nh3_Zcolumn'].values[ind,:]) #Eqn 10, B is zero for NH3.
	met['HRI'] = data['HRI'] #For new postfilter
	
	return met


#We only have height from surface, so we will use GC Boxheight diagnostic to do a regrid approximation.
def GC_to_sat_levels(GC_SPC, GC_bxheight, sat_edges):
	'''
	The provided edges for GEOS-Chem and the satellite should
	have dimension number of observations x number of edges
	'''
	GC_edges = np.concatenate([np.array([0]),np.cumsum(GC_bxheight,axis=1)]) #Try to put a 0 at the bottom of the cumsum
	# We want to account for the case when the GEOS-Chem surface
	# is above the satellite surface (altitude wise) or the GEOS-Chem
	# top is below the satellite top. We do this by adjusting the
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
	# Now map the GC NH3 to the satellite levels
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
		GC_col = (filt*sat_pressure_weight*sat_avker*GC_SPC)
	else:
		GC_col = (filt*sat_pressure_weight
				  *(sat_prior + sat_avker*(GC_SPC - sat_prior)))
	GC_col = GC_col.sum(axis=1)
	return GC_col 

class IASI_Translator(obsop.Observation_Translator):
	def __init__(self,verbose=1):
		super().__init__(verbose)
	#Save dictionary of dates for later use
	#Downloaded all of July 2019 with this code: 
	#wget https://cds-espri.ipsl.upmc.fr/iasibl2/iasi_nh3/V4.0.0/2019/07/IASI_METOPB_L2_NH3_201907{01..31}_ULB-LATMOS_V4.0.0.nc
	def initialReadDate(self):
		sourcedirs = self.spc_config['IASI_dirs']
		IASI_date_dict = {}
		for key in list(sourcedirs.keys()):
			sourcedir = sourcedirs[key]
			obs_list = glob(f'{sourcedir}/**/IASI_*.nc', recursive=True)
			obs_list.sort()
			IASI_date_dict[key] = [datetime.strptime(obs.split('_')[-3], "%Y%m%d") for obs in obs_list]
		with open(f"{self.scratch}/iasi_dates.pickle", 'wb') as handle:
			pickle.dump(IASI_date_dict, handle)
		return IASI_date_dict
	#Timeperiod is two datetime objects
	def globObs(self,species,timeperiod):
		sourcedir = self.spc_config['IASI_dirs'][species]
		if os.path.exists(f"{self.scratch}/iasi_dates.pickle"):
			with open(f"{self.scratch}/iasi_dates.pickle", 'rb') as handle:
				IASI_date_dict = pickle.load(handle)
		else:
			IASI_date_dict = self.initialReadDate()
		obs_dates = IASI_date_dict[species]
		obs_list = glob(f'{sourcedir}/**/IASI_*.nc', recursive=True)
		obs_list.sort()
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
			if self.spc_config['WHICH_TROPOMI_PRODUCT'] == 'ACMG':
				trop_obs.append(read_tropomi_acmg(obs,species,filterinfo,includeObsError=includeObsError))
			elif self.spc_config['WHICH_TROPOMI_PRODUCT'] == 'BLENDED':
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
		GC_col_data = obsop.getGCCols(GC,TROPOMI,species,self.spc_config,returninds=True,returnStateMet=True,GC_area=GC_area)
		GC_SPC = GC_col_data['GC_SPC']
		if 'Met_BXHEIGHT' not in GC_col_data:
			raise ValueError('ERROR: Met_BXHEIGHT not found. Ensure Met_BXHEIGHT is saved out in HistoryRC and listed in ens_config.')
		GC_bxheight = GC_col_data['Met_BXHEIGHT']
		i,j,t = GC_col_data['indices']
		if species=='CH4':
			GC_SPC*=1e9 #scale to mol/mol
			TROP_PRIOR = 1e9*(TROPOMI['methane_profile_apriori']/TROPOMI['dry_air_subcolumns'])
			synthetic_partial_columns = False
		elif species=='NO2':
			TROP_PRIOR=None
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
					additional_args_avgGC['other_fields_to_avg'][field] = TROPOMI[field]
		else:
			superobs=False
		#If user is taking super observations and wants to average before applying averaging kernel, do so now
		if superobs and (self.spc_config['SUPER_OBS_BEFORE_Hx'][specieskey]=="True"):
			presuperobs=True
			trop_agg_args = {'GC_SPC':GC_SPC,'GC_P':GC_P, 'T_press':TROPOMI['pressures'], 'ak':TROPOMI['column_AK'],'TROP_PW':TROP_PW,'TROP_PRIOR':TROP_PRIOR}
			agg_data = obsop.averageFieldsToGC(i,j,t,GC,TROPOMI[species],doSuperObs=doErrCalc,superObsFunction=superObsFunction,**additional_args_avgGC,**trop_agg_args)
			GC_on_sat_l = GC_to_sat_levels(agg_data['GC_SPC'], agg_data['GC_P'], agg_data['T_press'],species)
			GC_on_sat = apply_avker(agg_data['ak'],agg_data['TROP_PW'], GC_on_sat_l,agg_data['TROP_PRIOR'])
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


