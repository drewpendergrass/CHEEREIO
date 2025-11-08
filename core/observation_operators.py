from datetime import datetime
import settings_interface as si 
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np
import functools

def produceSuperObservationFunction(fname):
	if (fname is None) or (fname == "default"):
		def super_obs(mean_error,num_obs,errorCorr=0,min_error=0,transportError=0):
			return np.max([np.sqrt((mean_error**2 * (((1-errorCorr)/num_obs) + errorCorr) )+transportError**2),min_error])
	elif fname == "sqrt":
		def super_obs(mean_error,num_obs,errorCorr=0,min_error=0):
			return np.max([(mean_error * np.sqrt(((1-errorCorr)/num_obs) + errorCorr) ),min_error])
	elif fname == "constant":
		def super_obs(mean_error,num_obs,errorCorr=0,min_error=0):
			return mean_error
	else:
		raise ValueError('Superobservation error reduction function not recognized.')
	return super_obs

#All observations are passed as dictionaries. Filter info is a dictionary. Filter info keys indicate filter family.
#The OBSDATA dictionary must have a numpy array labeled "longitude", one labeled "latitude"
#and one labeled "utctime". The "utctime" entry must be formatted by ISO 8601 data time format
# ISO 8601 represents date and time by starting with the year, followed by the month, the day, 
#the hour, the minutes, seconds and milliseconds. For example, 2020-07-10 15:00:00.000, 
#represents the 10th of July 2020 at 3 p.m. Timezone is assumed to be UTC
def apply_filters(OBSDATA,filterinfo):
	to_keep = []
	filter_families = list(filterinfo.keys())
	if "MAIN" in filter_families:
		filter_data = filterinfo["MAIN"]
		filter_latitude = filter_data[0]
		if ~np.isnan(filter_latitude):
			to_keep.append(np.where(np.abs(OBSDATA['latitude'])<=filter_latitude)[0])
	if "OMI_NO2" in filter_families:
		filter_data = filterinfo["OMI_NO2"]
		filter_sza = filter_data[0]
		filter_cloud_radiance_frac = filter_data[1]
		filter_surface_albedo = filter_data[2]
		if ~np.isnan(filter_sza):
			to_keep.append(np.where(OBSDATA['SolarZenithAngle']<filter_sza)[0])
		if ~np.isnan(filter_cloud_radiance_frac):
			to_keep.append(np.where(OBSDATA['CloudRadianceFraction']<filter_cloud_radiance_frac)[0])
		if ~np.isnan(filter_surface_albedo):
			to_keep.append(np.where(OBSDATA['TerrainReflectivity']<filter_surface_albedo)[0])
	if "TROPOMI_CH4" in filter_families:
		filter_data = filterinfo["TROPOMI_CH4"]
		filter_blended_albedo = filter_data[0]
		filter_swir_albedo_low = filter_data[1]
		filter_swir_albedo_high = filter_data[2]
		filter_winter_lat = filter_data[3]
		filter_roughness = filter_data[4]
		filter_swir_aot = filter_data[5]
		if ~np.isnan(filter_blended_albedo):
			to_keep.append(np.where(OBSDATA['blended_albedo']<filter_blended_albedo)[0])
		if ~np.isnan(filter_swir_albedo_low):
			to_keep.append(np.where(OBSDATA['albedo_swir']>filter_swir_albedo_low)[0])
		if ~np.isnan(filter_swir_albedo_high):
			to_keep.append(np.where(OBSDATA['albedo_swir']<filter_swir_albedo_high)[0])
		if ~np.isnan(filter_winter_lat):
			months = OBSDATA['utctime'].astype('datetime64[M]').astype(int) % 12 + 1
			nh_winter = np.array([1,2,3,11,12])
			sh_winter = np.array([5,6,7,8,9])
			to_keep.append(np.where( ~( ( (OBSDATA['latitude']>filter_winter_lat)&(np.isin(months,nh_winter)) )| ( (OBSDATA['latitude']<(-1*filter_winter_lat))&(np.isin(months,sh_winter)) ) ) )[0])
		if ~np.isnan(filter_roughness):
			to_keep.append(np.where(OBSDATA['surface_elevation_sd']<filter_roughness)[0])
		if ~np.isnan(filter_swir_aot):
			to_keep.append(np.where(OBSDATA['swir_aot']<filter_swir_aot)[0])
	if len(to_keep)==0:
		return OBSDATA
	else:
		if len(to_keep)==1:
			to_keep = to_keep[0]
		else:
			to_keep = functools.reduce(np.intersect1d, to_keep)
		keys = list(OBSDATA.keys())
		for key in keys:
			if "TO_SKIP" in filter_families:
				if key in filterinfo["TO_SKIP"]:
					continue
			if len(np.shape(OBSDATA[key])) == 1:
				OBSDATA[key] = OBSDATA[key][to_keep]
			else:
				OBSDATA[key] = OBSDATA[key][to_keep,:]
		return OBSDATA

#GC is an xarray dataset representing GEOS-Chem results for time period of interest
def nearest_loc(GC,OBSDATA):
	# Find the grid box and time indices corresponding to OBSDATA obs
	# i index
	iGC = np.abs(GC.lon.values.reshape((-1, 1))
				 - OBSDATA['longitude'].reshape((1, -1)))
	iGC = iGC.argmin(axis=0)
	# j index
	jGC = np.abs(GC.lat.values.reshape((-1, 1))
				 - OBSDATA['latitude'].reshape((1, -1)))
	jGC = jGC.argmin(axis=0)
	# Time index
	tGC = np.abs(GC.time.values.reshape((-1, 1))
				 - OBSDATA['utctime'].astype('datetime64').reshape((1, -1)))
	tGC = tGC.argmin(axis=0)
	return iGC, jGC, tGC

def getGC_SPC(GC,species,spc_config):
	gc_version = float(spc_config['GC_VERSION'][0:-2]) #major plus minor version
	if gc_version>=14.1:
		spcconc_name = "SpeciesConcVV"
	else:
		spcconc_name = "SpeciesConc" #Starting in 14.1 we have to specify VV
	gc_overrides = si.checkGCSpeciesOverride(spc_config) #Here we check to see if GC stores the species under a different name. Only applies for N2O.
	if (len(gc_overrides)>0)&(species in gc_overrides):
		GC_SPC = GC[f'{spcconc_name}_{gc_overrides[species]}']
	else:
		GC_SPC = GC[f'{spcconc_name}_{species}']
	return GC_SPC


def getGCCols(GC,OBSDATA,species,spc_config,returninds=False,returnLevelEdge=True,returnStateMet=False,GC_area=None):
	i,j,t = nearest_loc(GC,OBSDATA)
	to_return = {}
	to_return['GC_SPC'] = getGC_SPC(GC,species,spc_config).values[t,:,j,i]
	if returnLevelEdge:
		to_return['GC_P'] = GC[f'Met_PEDGE'].values[t,:,j,i]
	if returnStateMet:
		for metcoll in spc_config['HistoryStateMetToSave']:
			if len(np.shape(GC[metcoll].values))==3:
				to_return[metcoll] = GC[metcoll].values[t,j,i]
			else:
				to_return[metcoll] = GC[metcoll].values[t,:,j,i]
	if GC_area is not None:
		to_return['GC_area']=GC_area.values[j,i]
	if returninds:
		to_return['indices'] = [i,j,t]
	return to_return

#No index, puts loc at GC grid values
#other_fields_to_avg is a dictionary with keys "variable name" and values of arrays.
def averageByGC(iGC, jGC, tGC, GC,GCmappedtoobs,obsvals,doSuperObs,superObsFunction=None,other_fields_to_avg=None, prescribed_error=None,prescribed_error_type=None, obsInstrumentError = None, modelTransportError = None, errorCorr = None,minError=None):
	#Check for nans and remove consistently. This really only comes up with IASI NH3 method 2
	valid_inds = np.logical_not(np.isnan(GCmappedtoobs))&np.logical_not(np.isnan(obsvals))
	iGC = iGC[valid_inds]
	jGC = jGC[valid_inds]
	tGC = tGC[valid_inds]
	GCmappedtoobs = GCmappedtoobs[valid_inds]
	obsvals = obsvals[valid_inds]
	#Now go about averaging.
	index = ((iGC+1)*100000000)+((jGC+1)*10000)+(tGC+1)
	unique_inds = np.unique(index)
	i_unique = np.floor(unique_inds/100000000).astype(int)-1
	j_unique = (np.floor(unique_inds/10000).astype(int) % 10000)-1
	t_unique = (unique_inds % 10000).astype(int)-1
	lonvals = GC.lon.values[i_unique]
	latvals = GC.lat.values[j_unique]
	timevals = GC.time.values[t_unique]
	av_len = len(unique_inds)
	gc_av = np.zeros(av_len)
	obs_av = np.zeros(av_len)
	obslat_av = np.zeros(av_len)
	obslon_av = np.zeros(av_len)
	obstime_av = np.zeros(av_len)
	num_av = np.zeros(av_len)
	if other_fields_to_avg is not None:
		#Initialize additional fields
		additional_fields = {}
		for field in other_fields_to_avg:
			additional_fields[field] = np.zeros(av_len)
	if doSuperObs:
		err_av = np.zeros(av_len)
	for count,ind in enumerate(unique_inds):
		indmatch = np.where(index==ind)[0]
		gc_av[count] = np.nanmean(GCmappedtoobs[indmatch])
		obs_av[count] = np.nanmean(obsvals[indmatch])
		obslat_av[count] = latvals[count]
		obslon_av[count] = lonvals[count]
		obstime_av[count] = timevals[count]
		num_av[count] = np.sum(~np.isnan(GCmappedtoobs[indmatch]))
		#Average the additional fields.
		if other_fields_to_avg is not None:
			for field in other_fields_to_avg:
				additional_fields[field][count] = np.nanmean(other_fields_to_avg[field][indmatch])
		if doSuperObs:
			#SuperObservation function selected by user
			obs_f = produceSuperObservationFunction(superObsFunction)
			obs_args = {}
			#Add arguments that are required by the function:
			if modelTransportError is not None:
				obs_args['transportError'] = modelTransportError
			if errorCorr is not None:
				obs_args['errorCorr'] = errorCorr
			if minError is not None: 
				obs_args['min_error'] = minError
			#Calculate mean error to be reduced by superobservation function. 
			if obsInstrumentError is not None:
				mean_err = np.nanmean(obsInstrumentError[indmatch])
			elif prescribed_error is not None:
				if prescribed_error_type == "relative":
					mean_err = obs_av[count]*prescribed_error
				elif prescribed_error_type == "absolute":
					mean_err = prescribed_error
				else:
					raise ValueError("Errors must be prescribed or included with observations; missing needed information")
			else:
				raise ValueError("Errors must be prescribed or included with observations; missing needed information")
			#Baseline model transport error doesn't average out; this is Zhen Qu's formulation; error correlation accounted for following Miyazaki et al 2012 and Eskes et al., 2003
			err_av[count] = obs_f(mean_error=mean_err,num_obs=num_av[count],**obs_args)
	to_return = ObsData(gc_av,obs_av,obslat_av,obslon_av,obstime_av,num_av=num_av)
	if other_fields_to_avg is not None:
		to_return.addData(**additional_fields)
	if doSuperObs:
		to_return.addData(err_av=err_av)
	return to_return

#Utility function for averageFieldsToGC (below)
def makeEmptyArrayDims(index_len, unique_len, field, returnAxis=False):
	to_return = list(field.shape)
	averaging_axis = np.where(np.array(to_return)==index_len)[0][0]
	if returnAxis:
		return averaging_axis
	else:
		to_return[averaging_axis] = unique_len
		return to_return

#Utility function for averageFieldsToGC (below)
def simple_slice(arr, inds, axis):
	sl = [slice(None)] * arr.ndim
	sl[axis] = inds
	return tuple(sl)

#In case where you want to average to GC grid before mapping to observations, use this function
#Can handle multidimension inputs (e.g. averaging kernels, prior profiles, columns) via kwargs. Autodetects axis with number of observations (defaulting to first axis)
def averageFieldsToGC(iGC, jGC, tGC, GC, obsvals, doSuperObs=False,superObsFunction=None,other_fields_to_avg=None, prescribed_error=None,prescribed_error_type=None, obsInstrumentError = None, modelTransportError = None, errorCorr = None,minError=None, **kwargs):
	index = ((iGC+1)*100000000)+((jGC+1)*10000)+(tGC+1)
	unique_inds = np.unique(index)
	i_unique = np.floor(unique_inds/100000000).astype(int)-1
	j_unique = (np.floor(unique_inds/10000).astype(int) % 10000)-1
	t_unique = (unique_inds % 10000).astype(int)-1
	lonvals = GC.lon.values[i_unique]
	latvals = GC.lat.values[j_unique]
	timevals = GC.time.values[t_unique]
	av_len = len(unique_inds)
	obsaxis = {k:makeEmptyArrayDims(len(iGC),av_len,kwargs[k],returnAxis=True) for k in kwargs}
	to_return = {k:np.zeros(makeEmptyArrayDims(len(iGC),av_len,kwargs[k])) for k in kwargs}
	to_return['lat'] = np.zeros(av_len)
	to_return['lon'] = np.zeros(av_len)
	to_return['time'] = np.zeros(av_len)
	to_return['num'] = np.zeros(av_len)
	to_return['obs'] = np.zeros(av_len)
	if other_fields_to_avg is not None:
		#Initialize additional fields
		to_return['additional_fields'] = {}
		for field in other_fields_to_avg:
			to_return['additional_fields'][field] = np.zeros(av_len)
	if doSuperObs:
		to_return['err'] = np.zeros(av_len)
	for count,ind in enumerate(unique_inds):
		indmatch = np.where(index==ind)[0]
		to_return['lat'][count] = latvals[count]
		to_return['lon'][count] = lonvals[count]
		to_return['time'][count] = timevals[count]
		to_return['num'][count] = len(indmatch)
		to_return['obs'][count] = np.mean(obsvals[indmatch])
		#Average the additional fields.
		if other_fields_to_avg is not None:
			for field in other_fields_to_avg:
				to_return['additional_fields'][field][count] = np.mean(other_fields_to_avg[field][indmatch])
		for k in kwargs:
			to_return[k][simple_slice(to_return[k],count,obsaxis[k])] = np.mean(kwargs[k][simple_slice(kwargs[k],indmatch,obsaxis[k])],axis=obsaxis[k])
		if doSuperObs:
			#SuperObservation function selected by user
			obs_f = produceSuperObservationFunction(superObsFunction)
			obs_args = {}
			#Add arguments that are required by the function:
			if modelTransportError is not None:
				obs_args['transportError'] = modelTransportError
			if errorCorr is not None:
				obs_args['errorCorr'] = errorCorr
			if minError is not None: 
				obs_args['min_error'] = minError
			#Calculate mean error to be reduced by superobservation function. 
			if obsInstrumentError is not None:
				mean_err = np.mean(obsInstrumentError[indmatch])
			elif prescribed_error is not None:
				if prescribed_error_type == "relative":
					mean_err = to_return['obs'][count]*prescribed_error
				elif prescribed_error_type == "absolute":
					mean_err = prescribed_error
				else:
					raise ValueError("Errors must be prescribed or included with observations; missing needed information")
			else:
				raise ValueError("Errors must be prescribed or included with observations; missing needed information")
			#Baseline model transport error doesn't average out; this is Zhen Qu's formulation; error correlation accounted for following Miyazaki et al 2012 and Eskes et al., 2003
			to_return['err'][count] = obs_f(mean_error=mean_err,num_obs=to_return['num'][count],**obs_args)
	return to_return

class Observation_Translator(object):
	def __init__(self,verbose=1):
		self.verbose = verbose
		self.spc_config = si.getSpeciesConfig()
		self.scratch = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/scratch"
	#The globObs function returns a dictionary of observations that are within a user-specified timeperiod, and, if specified
	#that take place at a user-specified interval. The inherited function must have this signature.
	#Please note that the "specieskey" variable MUST be the key in the dictionary OBSERVED_SPECIES in ens_config.
	#The returned dictionary must have keys for "latitude", "longitude", and "utctime",
	#where UTC time is an ISO 8601 date time string
	def getObservations(self,specieskey,timeperiod, interval=None, includeObsError=False):
		#Returns a specifically formatted dictionary (see above for instructions)
		raise NotImplementedError
	#The function that gets the comparison between GEOS-Chem and the observations (OBSDATA, formatted in a dictionary as above).
	#Please note that the "specieskey" variable MUST be the key in the dictionary OBSERVED_SPECIES in ens_config.
	#). Inherited function must have this signature and return an ObsData object.
	def gcCompare(self,specieskey,OBSDATA,GC,GC_area=None,doErrCalc=True,useObserverError=False, prescribed_error=None,prescribed_error_type=None,transportError = None, errorCorr = None,minError=None):
		#Returns an ObsData object
		raise NotImplementedError

#Class containing simulated and actual observations, feeds directly into HIST_Ens and then into the LETKF routine in Assimilator.
#If you are adding a piece of additional data not already supported (e.g. you want to pass data along and save to plot in postprocessing, like TROPOMI albedo),
#you will have to implement this in HIST_Ens and relevant postprocessing scripts
class ObsData(object):
	def __init__(self,gccol,obscol,obslat,obslon,obstime,**additional_data):
		self.gccol = gccol #Note: this is sometimes a set of columns, as when stored in bigYDict in HIST_Ens
		self.obscol = obscol
		self.obslat = obslat
		self.obslon = obslon
		self.obstime = obstime
		self.additional_data = additional_data #dictionary of additional data for various observations operators (e.g. albedo, averaging information)
	def getGCCol(self):
		return self.gccol
	def setGCCol(self,gccol):
		self.gccol = gccol
	def getObsCol(self):
		return self.obscol
	def getCols(self):
		return [self.gccol,self.obscol]
	def getLatLon(self):
		return [self.obslat,self.obslon]
	def getTime(self):
		return self.obstime
	def addData(self,**data_to_add):
		for key in data_to_add.keys():
			self.additional_data[key] = data_to_add[key]
	def getDataByKey(self,key):
		if type(key) is list:
			to_return = [self.additional_data[keyval] for keyval in key]
		else:
			to_return = self.additional_data[key]
		return to_return

