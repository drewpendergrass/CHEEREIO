#Some of this code is adapted from the excellent work done by Hannah Nesser, Melissa Sulprizio, Daniel Varon, Lu Shen
#and other contributors to the integrated methane inversion workflow. I am indebted to them and to Elise Penn and Alba Lorente for 
#explaining much of TROPOMI.

from datetime import datetime
import toolbox as tx
from glob import glob
import pickle
import os.path
import xarray as xr
import numpy as np

def read_tropomi(filename, species, filterinfo=None):
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
	
	if filterinfo:
		met = apply_filters(met,filterinfo)

	return met

def apply_filters(TROPOMI,filterinfo):
	to_keep = []
	if filterinfo[0] == "TROPOMI_CH4":
		filter_blended_albedo = filterinfo[1]
		filter_swir_albedo_low = filterinfo[2]
		filter_swir_albedo_high = filterinfo[3]
		filter_winter_lat = filterinfo[4]
		filter_roughness = filterinfo[5]
		filter_swir_aot = filterinfo[6]
		if ~np.isnan(filter_blended_albedo):
			to_keep.append(np.where(TROPOMI['blended_albedo']<filter_blended_albedo)[0])
		if ~np.isnan(filter_swir_albedo_low):
			to_keep.append(np.where(TROPOMI['albedo_swir']>filter_swir_albedo_low)[0])
		if ~np.isnan(filter_swir_albedo_high):
			to_keep.append(np.where(TROPOMI['albedo_swir']<filter_swir_albedo_high)[0])
		if ~np.isnan(filter_winter_lat):
			months = TROPOMI['utctime'].astype('datetime64[M]').astype(int) % 12 + 1
			nh_winter = np.array([1,2,3,11,12])
			sh_winter = np.array([5,6,7,8,9])
			to_keep.append(np.where( ~( ( (TROPOMI['latitude']>filter_winter_lat)&(np.isin(months,nh_winter)) )| ( (TROPOMI['latitude']<(-1*filter_winter_lat))&(np.isin(months,sh_winter)) ) ) )[0])
		if ~np.isnan(filter_roughness):
			to_keep.append(np.where(TROPOMI['surface_elevation_sd']<filter_roughness)[0])
		if ~np.isnan(filter_swir_aot):
			to_keep.append(np.where(TROPOMI['swir_aot']<filter_swir_aot)[0])
	if len(to_keep)==0:
		return TROPOMI
	else:
		if len(to_keep)==1:
			to_keep = to_keep[0]
		elif len(to_keep)==2:
			to_keep = np.intersect1d(to_keep[0],to_keep[1])
		elif len(to_keep)==3:
			to_keep = np.intersect1d(np.intersect1d(to_keep[0],to_keep[1]),to_keep[2])
		elif len(to_keep)==4:
			to_keep = np.intersect1d(np.intersect1d(to_keep[0],to_keep[1]),np.intersect1d(to_keep[2],to_keep[3]))
		elif len(to_keep)==5:
			to_keep = np.intersect1d(np.intersect1d(np.intersect1d(to_keep[0],to_keep[1]),np.intersect1d(to_keep[2],to_keep[3])),to_keep[4])
		elif len(to_keep)==6:
			to_keep = np.intersect1d(np.intersect1d(np.intersect1d(to_keep[0],to_keep[1]),np.intersect1d(to_keep[2],to_keep[3])),np.intersect1d(to_keep[4],to_keep[5]))
		keys = list(TROPOMI.keys())
		for key in keys:
			if len(np.shape(TROPOMI[key])) == 1:
				TROPOMI[key] = TROPOMI[key][to_keep]
			else:
				TROPOMI[key] = TROPOMI[key][to_keep,:]
		return TROPOMI

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

def getGCCols(GC,TROPOMI,species,returninds=False,returnStateMet=False,GC_area=None):
	i,j,t = nearest_loc(GC,TROPOMI)
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

#Takes index in form index = (iGC*100000000)+(jGC*10000)+tGC
def averageByGCTROPOMIlocs(index,GConsat,satvals,satlat,satlon,sattime):
	unique_inds = np.unique(index)
	av_len = len(unique_inds)
	gc_av = np.zeros(av_len)
	sat_av = np.zeros(av_len)
	satlat_av = np.zeros(av_len)
	satlon_av = np.zeros(av_len)
	sattime_av = np.zeros(av_len)
	num_av = np.zeros(av_len)
	delta_sec = np.timedelta64(1, 's')
	epoch = np.datetime64('1970-01-01T00:00:00')
	for count,ind in enumerate(unique_inds):
		indmatch = np.where(index==ind)[0]
		gc_av[count] = np.mean(GConsat[indmatch])
		sat_av[count] = np.mean(satvals[indmatch])
		satlat_av[count] = np.mean(satlat[indmatch])
		satlon_av[count] = np.mean(satlon[indmatch])
		dt = sattime[indmatch].astype('datetime64')
		epoch_sec = (dt - epoch) / delta_sec
		epoch_sec_mean = np.mean(epoch_sec)
		sattime_av[count] = epoch + np.timedelta64(int(epoch_sec_mean), 's')
		num_av[count] = len(indmatch)
	return [gc_av,sat_av,satlat_av,satlon_av,sattime_av,num_av]

#No index, puts loc at GC grid values
def averageByGC(iGC, jGC, tGC, GC,GConsat,satvals,satlat,satlon,sattime,albedo_swir=None,albedo_nir=None,blended_albedo=None):
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
	sat_av = np.zeros(av_len)
	satlat_av = np.zeros(av_len)
	satlon_av = np.zeros(av_len)
	sattime_av = np.zeros(av_len)
	num_av = np.zeros(av_len)
	if albedo_swir is not None:
		swir_av = np.zeros(av_len)
		nir_av = np.zeros(av_len)
		blended_av = np.zeros(av_len)
	for count,ind in enumerate(unique_inds):
		indmatch = np.where(index==ind)[0]
		gc_av[count] = np.mean(GConsat[indmatch])
		sat_av[count] = np.mean(satvals[indmatch])
		satlat_av[count] = latvals[count]
		satlon_av[count] = lonvals[count]
		sattime_av[count] = timevals[count]
		num_av[count] = len(indmatch)
		if albedo_swir is not None:
			swir_av[count] = np.mean(albedo_swir[indmatch])
			nir_av[count] = np.mean(albedo_nir[indmatch])
			blended_av[count] = np.mean(blended_albedo[indmatch])
	if albedo_swir is not None:
		return [gc_av,sat_av,satlat_av,satlon_av,sattime_av,num_av,swir_av,nir_av,blended_av]
	else:
		return [gc_av,sat_av,satlat_av,satlon_av,sattime_av,num_av]

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
	def getSatellite(self,species,timeperiod, interval=None):
		obs_list = self.globObs(species,timeperiod,interval)
		trop_obs = []
		if species=='CH4':
			if self.spc_config['TROPOMI_CH4_FILTERS']=="True":
				filterinfo = ["TROPOMI_CH4",float(self.spc_config['TROPOMI_CH4_filter_blended_albedo']),float(self.spc_config['TROPOMI_CH4_filter_swir_albedo_low']),float(self.spc_config['TROPOMI_CH4_filter_swir_albedo_high']),float(self.spc_config['TROPOMI_CH4_filter_winter_lat']),float(self.spc_config['TROPOMI_CH4_filter_roughness']),float(self.spc_config['TROPOMI_CH4_filter_swir_aot'])]
			else:
				filterinfo = None
		else:
			filterinfo = None
		for obs in obs_list:
			trop_obs.append(read_tropomi(obs,species,filterinfo))
		met = {}
		for key in list(trop_obs[0].keys()):
			met[key] = np.concatenate([metval[key] for metval in trop_obs])
		return met
	def gcCompare(self,species,timeperiod,TROPOMI,GC,GC_area=None,saveAlbedo=False):
		if species=='CH4':
			TROP_PRIOR = 1e9*(TROPOMI['methane_profile_apriori']/TROPOMI['dry_air_subcolumns'])
			synthetic_partial_columns = False
		elif species=='NO2':
			TROP_PRIOR=None
			synthetic_partial_columns = True
		TROP_PW = (-np.diff(TROPOMI['pressures'])/(TROPOMI['pressures'][:, 0] - TROPOMI['pressures'][:, -1])[:, None])
		returnStateMet = self.spc_config['SaveStateMet']=='True'
		if returnStateMet:
			GC_SPC,GC_P,GC_M,GC_area,i,j,t = getGCCols(GC,TROPOMI,species,returninds=True,returnStateMet=returnStateMet,GC_area=GC_area)
		else:
			GC_SPC,GC_P,GC_area,i,j,t = getGCCols(GC,TROPOMI,species,returninds=True,returnStateMet=returnStateMet,GC_area=GC_area)
		if species=='CH4':
			GC_SPC*=1e9 #scale to mol/mol
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
			if saveAlbedo:
				gc_av,sat_av,satlat_av,satlon_av,sattime_av,num_av,swir_av,nir_av,blended_av = averageByGC(i,j,t,GC,GC_on_sat,TROPOMI[species],TROPOMI['latitude'],TROPOMI['longitude'],TROPOMI['utctime'],TROPOMI['albedo_swir'],TROPOMI['albedo_nir'],TROPOMI['blended_albedo'])
				toreturn = [gc_av,sat_av,satlat_av,satlon_av,sattime_av,num_av,swir_av,nir_av,blended_av]
			else:
				gc_av,sat_av,satlat_av,satlon_av,sattime_av,num_av = averageByGC(i,j,t,GC,GC_on_sat,TROPOMI[species],TROPOMI['latitude'],TROPOMI['longitude'],TROPOMI['utctime'])
				toreturn = [gc_av,sat_av,satlat_av,satlon_av,sattime_av,num_av]
		else:
			if saveAlbedo:
				toreturn = [GC_on_sat,TROPOMI[species],TROPOMI['latitude'],TROPOMI['longitude'],TROPOMI['utctime'],TROPOMI['albedo_swir'],TROPOMI['albedo_nir'],TROPOMI['blended_albedo']]
			else:
				toreturn = [GC_on_sat,TROPOMI[species],TROPOMI['latitude'],TROPOMI['longitude'],TROPOMI['utctime']]
		return toreturn