from datetime import datetime
import settings_interface as si 
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np
import functools


#All observations are passed as dictionaries. Filter info is a list
#The OBSDATA dictionary must have a numpy array labeled "longitude", one labeled "latitude"
#and one labeled "utctime". The "utctime" entry must be formatted by ISO 8601 data time format
# ISO 8601 represents date and time by starting with the year, followed by the month, the day, 
#the hour, the minutes, seconds and milliseconds. For example, 2020-07-10 15:00:00.000, 
#represents the 10th of July 2020 at 3 p.m. Timezone is assumed to be UTC
def apply_filters(OBSDATA,filterinfo):
	to_keep = []
	if filterinfo[0] == "TROPOMI_CH4":
		filter_blended_albedo = filterinfo[1]
		filter_swir_albedo_low = filterinfo[2]
		filter_swir_albedo_high = filterinfo[3]
		filter_winter_lat = filterinfo[4]
		filter_roughness = filterinfo[5]
		filter_swir_aot = filterinfo[6]
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

def getGCCols(GC,OBSDATA,species,returninds=False,returnStateMet=False,GC_area=None):
	i,j,t = nearest_loc(GC,OBSDATA)
	to_return = [GC[f'SpeciesConc_{species}'].values[t,:,j,i],GC[f'Met_PEDGE'].values[t,:,j,i]]
	if returnStateMet:
		to_return.append(GC[f'Met_AD'].values[t,:,j,i])
	if GC_area is not None:
		to_return.append(GC_area.values[j,i])
	else:
		to_return.append(None)
	if returninds:
		to_return = to_return + [i,j,t]
	return to_return

#No index, puts loc at GC grid values
def averageByGC(iGC, jGC, tGC, GC,GCmappedtoobs,obsvals,albedo_swir=None,albedo_nir=None,blended_albedo=None, satError = None, modelTransportError = None, errorCorr = None):
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
	if albedo_swir is not None:
		swir_av = np.zeros(av_len)
		nir_av = np.zeros(av_len)
		blended_av = np.zeros(av_len)
	if satError is not None:
		err_av = np.zeros(av_len)
	for count,ind in enumerate(unique_inds):
		indmatch = np.where(index==ind)[0]
		gc_av[count] = np.mean(GCmappedtoobs[indmatch])
		obs_av[count] = np.mean(obsvals[indmatch])
		obslat_av[count] = latvals[count]
		obslon_av[count] = lonvals[count]
		obstime_av[count] = timevals[count]
		num_av[count] = len(indmatch)
		if albedo_swir is not None:
			swir_av[count] = np.mean(albedo_swir[indmatch])
			nir_av[count] = np.mean(albedo_nir[indmatch])
			blended_av[count] = np.mean(blended_albedo[indmatch])
		if satError is not None:
			#Baseline model transport error doesn't average out; this is Zhen Qu's formulation; error correlation accounted for following Miyazaki et al 2012 and Eskes et al., 2003
			err_av[count] = (np.mean(satError[indmatch]) * np.sqrt(((1-errorCorr)/num_av[count]) + errorCorr) )+modelTransportError
	to_return = ObsData(gc_av,obs_av,obslat_av,obslon_av,obstime_av,num_av=num_av)
	if albedo_swir is not None:
		to_return.addData(swir_av=swir_av,nir_av=nir_av,blended_av=blended_av)
	if satError is not None:
		to_return.addData(err_av=err_av)
	return to_return

class Observation_Translator(object):
	def __init__(self,verbose=1):
		self.verbose = verbose
		self.spc_config = si.getSpeciesConfig()
		self.scratch = f"{self.spc_config['MY_PATH']}/{self.spc_config['RUN_NAME']}/scratch"
	#The globObs function returns a dictionary of observations that are within a user-specified timeperiod, and, if specified
	#that take place at a user-specified interval. The inherited function must have this signature.
	#The returned dictionary must have keys for "latitude", "longitude", and "utctime",
	#where UTC time is an ISO 8601 date time string
	def getObservations(self,species,timeperiod, interval=None, calcError=False):
		#Returns a specifically formatted dictionary (see above for instructions)
		raise NotImplementedError
	#The function that gets the comparison between GEOS-Chem and the observations (OBSDATA, formatted in a dictionary as above).
	#). Inherited function must have this signature and return an ObsData object.
	def gcCompare(self,species,OBSDATA,GC,GC_area=None,saveAlbedo=False,saveError=False, transportError = 0):
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

