#Tools for reading IASI NH3. Lieven Clarisse and Shixian Zhai helped me out on this (DCP 2024 Dec 18)

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
		filterinfo [list] : Currently only support the MAIN filters.
		includeObsError [bool] : if True, combine systematic and random uncertainty in one number
	Returns
		met	  [dict] : Dictionary of important variables from IASI:
							- species column
 							- Latitude
							- Longitude
	   						- QA value
 							- UTC time
							- Averaging kernel
							- HRI
	"""
	if (species != 'NH3'):
		raise ValueError('Only supports NH3 at this time.')
	# Initialize list for IASI data
	met = {}
	# Store species, QA, lat, lon, time, averaging kernel
	if filename=='test': #for developer only.
		filename='/n/holylfs05/LABS/jacob_lab/Users/drewpendergrass/IASI_NH3/v4/2019_07/IASI_METOPB_L2_NH3_20190701_ULB-LATMOS_V4.0.0.nc'
	data = xr.open_dataset(filename)
	qa = data['prefilter'].values
	ind = np.where(qa==1)[0] #Keep only those passing prefilter; will apply postfilter later
	#Store data with prefilter applied
	met['qa_prefilter'] = qa[ind]
	met['utctime'] = data['AERIStime'].values[ind]
	met['latitude'] = data['latitude'].values[ind]
	met['longitude'] = data['longitude'].values[ind]
	met['level_edge'] = data['levels'].values #in km, same for all obs (length 15)
	met['level_middle'] = data['midlevels'].values #in km, same for all obs (length 14)
	met[species] = data['nh3_total_column'].values[ind]
	#We are applying Method 2 from the avkReadMe (accompanying Clarisse et al 2023). 
	numerator = met[species]/data['nh3_AvKnorm'].values[ind]
	met['column_AK'] = numerator[:, np.newaxis]*(data['nh3_Zcolumn'].values[ind,:]**-1) #Eqn 10, B is zero for NH3.
	met['HRI'] = data['HRI'].values #For new postfilter
	data.close()
	#Apply main filter
	if filterinfo is not None:
		filterinfo["TO_SKIP"] = ['level_edge','level_middle'] #exclude these; no need to filter.
		met = obsop.apply_filters(met,filterinfo)
	return met

#We only have height from surface, so we will use GC Boxheight diagnostic to do a regrid approximation.
#Sat edges from iasi are 1d vector (same for all points) so we will convert
def GC_to_sat_levels(GC_SPC, GC_bxheight, sat_edges):
	#Make conversions into compatible matrices.
	GC_edges = np.column_stack([np.zeros(GC_bxheight.shape[0]),np.cumsum(GC_bxheight,axis=1)])*1e-3 #Put a 0 at the bottom of the cumsum, convert to km
	sat_edges = np.tile(sat_edges,(GC_bxheight.shape[0],1))
	'''
	The provided edges for GEOS-Chem and the satellite should
	have dimension number of observations x number of edges
	'''
	# We want to account for the case when the GEOS-Chem surface
	# is above the satellite surface (altitude wise) or the GEOS-Chem
	# top is below the satellite top. We do this by adjusting the
	# GEOS-Chem surface pressure up to the IASI surface pressure
	idx_bottom = np.less(GC_edges[:, 0], sat_edges[:, 0])
	idx_top = np.greater(GC_edges[:, -1], sat_edges[:, -1])
	GC_edges[idx_bottom, 0] = sat_edges[idx_bottom, 0]
	GC_edges[idx_top, -1] = sat_edges[idx_top, -1]
	# Define vectors that give the "low" and "high" altitude
	# values for each GEOS-Chem and satellite layer.
	GC_lo = GC_edges[:, 1:][:, :, None]
	GC_hi = GC_edges[:, :-1][:, :, None]
	sat_lo = sat_edges[:, 1:][:, None, :]
	sat_hi = sat_edges[:, :-1][:, None, :]
	# Get the indices where the GC-to-satellite mapping, which is
	# a nobs x ngc x nsat matrix, is non-zero
	idx = (np.less_equal(sat_lo, GC_hi) & np.greater_equal(sat_hi, GC_lo))
	# Find the fraction of each GC level that contributes to each
	# IASI level. We should first divide (to normalize) and then
	# multiply (to apply the map to the column) by the GC pressure
	# difference, but we exclude this (since it's the same as x1).
	GC_to_sat = np.minimum(sat_hi, GC_hi) - np.maximum(sat_lo, GC_lo)
	GC_to_sat[~idx] = 0
	# Now map the GC NH3 to the satellite levels
	GC_on_sat = (GC_to_sat*GC_SPC[:, :, None]).sum(axis=1)
	GC_on_sat = GC_on_sat/GC_to_sat.sum(axis=1)
	return [GC_on_sat,AD_on_sat]

#Map GC data (in molec/cm3), already on sat levels, to be equivalent to IASI. Both are changed. 
#Returns GC data in mol/m2 (Mm) and IASI column retrieved with GC prior (Xm), which are comparable.
def apply_avker(GC_on_sat_l,IASI_BXHEIGHT,IASI_AK,IASI_col):
	#Convert GC from mol/m3 to mol/m2 with IASI box height
	Mzm=GC_on_sat_l*IASI_BXHEIGHT # convert to column mol/m2 of dry air but keep binned. M_z^m in eq 11
	Mm=np.copy(Mzm).sum(axis=1) #Sum for total column. M^m in eq 11
	#Normalized modelled profile m_z defined in eq 11: mz = (M_z^m - Bz)/(M^m - B). No B for ammonia
	mz = Mzm / Mm
	#Calculate Xm (IASI column retrieved with GC normalized profile), defined in eq 13.
	Xm = IASI_col/(mz*IASI_AK).sum(axis=1)
	return [Mm,Xm]

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
		obs_list = [obs for obs,t in zip(obs_list,obs_dates) if (t>=timeperiod[0]) and (t<=timeperiod[1])]
		return obs_list
	def getObservations(self,specieskey,timeperiod, interval=None, includeObsError=False):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		obs_list = self.globObs(species,timeperiod)
		iasi_obs = []
		filterinfo = {}
		if specieskey in list(self.spc_config["filter_obs_poleward_of_n_degrees"].keys()):
			filterinfo['MAIN']=[float(self.spc_config["filter_obs_poleward_of_n_degrees"][specieskey])]
		for obs in obs_list:
			iasi_obs.append(read_iasi(obs,species,filterinfo,includeObsError=includeObsError))
		met = {}
		for key in list(iasi_obs[0].keys()):
			met[key] = np.concatenate([metval[key] for metval in iasi_obs])
		return met
	def gcCompare(self,specieskey,IASI,GC,GC_area=None,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):
		species = self.spc_config['OBSERVED_SPECIES'][specieskey]
		extra_obsdata_to_save = self.spc_config['EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y'][specieskey]
		GC_col_data = obsop.getGCCols(GC,IASI,species,self.spc_config,returninds=True,returnLevelEdge=False,returnStateMet=True,GC_area=GC_area)
		GC_SPC = GC_col_data['GC_SPC']
		if ('Met_BXHEIGHT' not in GC_col_data) or ('Met_AIRDEN' not in GC_col_data):
			raise ValueError('ERROR: missing met fields. Ensure Met_BXHEIGHT and Met_AIRDEN are saved out in HistoryRC and listed in ens_config.')
		GC_bxheight = GC_col_data['Met_BXHEIGHT']
		GC_AIRDEN = GC_col_data['Met_AIRDEN']
		IASI_EDGES = IASI['level_edge']*1e-3 #km to m
		IASI_BXHEIGHT = IASI_EDGES[1::]-IASI_EDGES[0:-1] #calculate box height
		#Convert GC_SPC from mol/mol to mol/m3 (still fine for regridding)
		GC_SPC = GC_SPC*(GC_AIRDEN / 0.028964) # apply model layer dry air density in mol/m3 (air molar mass: 0.028964 kg/mol), Met_AIRDEN(kg/m3)
		i,j,t = GC_col_data['indices']
		#If super observations, prep the error calculation
		if self.spc_config['AV_TO_GC_GRID'][specieskey]=="True":
			superobs=True
			superObsFunction = self.spc_config['SUPER_OBSERVATION_FUNCTION'][specieskey]
			additional_args_avgGC = {}
			if doErrCalc:
				if useObserverError:
					additional_args_avgGC['obsInstrumentError'] = IASI['Error']
					additional_args_avgGC['modelTransportError'] = transportError
				elif prescribed_error is not None:
					additional_args_avgGC['prescribed_error'] = prescribed_error
					additional_args_avgGC['prescribed_error_type'] = prescribed_error_type
				if minError is not None:
					additional_args_avgGC['minError'] = minError
				if errorCorr is not None:
					additional_args_avgGC['errorCorr'] = errorCorr
			additional_args_avgGC['other_fields_to_avg'] = {}
			#save HRI for postfilter
			additional_args_avgGC['other_fields_to_avg']['HRI'] = IASI['HRI']
			for field in extra_obsdata_to_save:
				additional_args_avgGC['other_fields_to_avg'][field] = IASI[field]
		else:
			superobs=False
		#If user is taking super observations and wants to average before applying averaging kernel, do so now
		if superobs and (self.spc_config['SUPER_OBS_BEFORE_Hx'][specieskey]=="True"):
			presuperobs=True
			iasi_agg_args = {'GC_SPC':GC_SPC, 'ak':IASI['column_AK'],'GC_bxheight':GC_bxheight}
			agg_data = obsop.averageFieldsToGC(i,j,t,GC,IASI[species],doSuperObs=doErrCalc,superObsFunction=superObsFunction,**additional_args_avgGC,**iasi_agg_args)
			GC_on_sat_l = GC_to_sat_levels(agg_data['GC_SPC'],agg_data['GC_bxheight'], IASI_EDGES)
			GC_SPC_final,IASI_SPC_final = apply_avker(GC_on_sat_l,IASI_BXHEIGHT,agg_data['ak'],agg_data['obs'])
		else:
			presuperobs=False
			GC_on_sat_l = GC_to_sat_levels(GC_SPC, GC_bxheight, IASI_EDGES)
			GC_SPC_final,IASI_SPC_final = apply_avker(GC_on_sat_l,IASI_BXHEIGHT,IASI['column_AK'],IASI[species])
		nan_indices = np.argwhere(np.isnan(GC_SPC_final))
		GC_SPC_final = np.nan_to_num(GC_SPC_final)
		if superobs:
			if presuperobs:
				toreturn = obsop.ObsData(GC_SPC_final,IASI_SPC_final,agg_data['lat'],agg_data['lon'], agg_data['time'],num_av=agg_data['num'])
				if doErrCalc:
					toreturn.addData(err_av=agg_data['err'])
				if 'additional_fields' in agg_data:
					toreturn.addData(**agg_data['additional_fields'])
			else:
				toreturn = obsop.averageByGC(i,j,t,GC,GC_SPC_final,IASI_SPC_final,doSuperObs=doErrCalc,superObsFunction=superObsFunction,**additional_args_avgGC)
		else:
			timevals = GC.time.values[t]
			toreturn = obsop.ObsData(GC_SPC_final,IASI_SPC_final,IASI['latitude'],IASI['longitude'],timevals)
			data_to_add = {}
			data_to_add['HRI'] = IASI['HRI']
			#If saving extra fields, add them here
			for field in extra_obsdata_to_save:
				data_to_add[field] = IASI[field]
			toreturn.addData(**data_to_add)
			if doErrCalc and useObserverError:
				toreturn.addData(err_av=IASI['Error'])
		#Apply postfilter.
		SFm_inv_abs = np.abs(toreturn.getDataByKey('HRI')/(toreturn.getObsCol()*6.02214076e19))**-1 #convert to molec/cm2 for threshholding.
		valid_after_postfilter = (SFm_inv_abs<1.5e16) & np.logical_not((np.abs(toreturn.getDataByKey('HRI'))>1.5) & (toreturn.getObsCol() < 0)) #These points are all valid after threshholding.
		#go through elements of obsdata and apply filter.
		toreturn.gccol = toreturn.gccol[valid_after_postfilter]
		toreturn.obscol = toreturn.obscol[valid_after_postfilter]
		toreturn.obslat = toreturn.obslat[valid_after_postfilter]
		toreturn.obslon = toreturn.obslon[valid_after_postfilter]
		toreturn.obstime = toreturn.obstime[valid_after_postfilter]
		for field in toreturn.additional_data: #Apply filter to every other item with right dimensionality
			to_edit = to_return.getDataByKey(field)
			if len(to_edit)==len(valid_after_postfilter):
				to_return.addData(**{field:to_edit[valid_after_postfilter]})
		return toreturn

#Testing zone.
# import testing_tools as tt
# a = tt.makeAssimilator(date='20190704_0000')

