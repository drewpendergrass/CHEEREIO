import numpy as np 
import xarray as xr
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter
from datetime import datetime,timedelta
import sys
import pickle
import pandas as pd
sys.path.append('../core')
import settings_interface as si 

def getSpcConcName():
	spc_config = si.getSpeciesConfig()
	gc_version = float(spc_config['GC_VERSION'][0:-2]) #major plus minor version
	if gc_version>=14.1:
		spcconc_name = "SpeciesConcVV"
	else:
		spcconc_name = "SpeciesConc" #Starting in 14.1 we have to specify VV
	return spcconc_name


def globDirs(ensemble_dir,removeNature=False,includeOutputDir=False):
	subdirs = glob(f"{ensemble_dir}/*/")
	subdirs.remove(f"{ensemble_dir}/logs/")
	subdirs.sort()
	dirnames = [d.split('/')[-2] for d in subdirs]
	subdir_numbers = [int(n.split('_')[-1]) for n in dirnames]
	if removeNature:
		try:
			ind = subdir_numbers.index(0)
			del subdirs[ind]
			del dirnames[ind]
			del subdir_numbers[ind]
		except ValueError:
			pass
	if includeOutputDir:
		subdirs = [d+'OutputDir' for d in subdirs]
	return([subdirs,dirnames,subdir_numbers])

def globSubDir(hist_dir,timeperiod=None,hourlysub = 6):
	specconc_list = glob(f'{hist_dir}/GEOSChem.SpeciesConc*.nc4')
	specconc_list.sort()
	ts = [datetime.strptime(spc.split('.')[-2][0:13], "%Y%m%d_%H%M") for spc in specconc_list]
	if timeperiod:
		specconc_list = [spc for spc,t in zip(specconc_list,ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
	specconc_list = [spc for spc,t in zip(specconc_list,ts) if t.hour % hourlysub == 0]
	return specconc_list

def globSubDirLevelEdge(hist_dir,timeperiod=None,hourlysub = 6):
	edgeconc_list = glob(f'{hist_dir}/GEOSChem.LevelEdgeDiags*.nc4')
	edgeconc_list.sort()
	ts = [datetime.strptime(edge.split('.')[-2][0:13], "%Y%m%d_%H%M") for edge in edgeconc_list]
	if timeperiod:
		edgeconc_list = [edge for edge,t in zip(edgeconc_list,ts) if (t>=timeperiod[0]) and (t<timeperiod[1])]
	edgeconc_list = [edge for edge,t in zip(edgeconc_list,ts) if t.hour % hourlysub == 0]
	return edgeconc_list

def combineScaleFactors(ensemble_dir,output_dir,timeperiod=None,flag_snapshot=False,return_not_write=False):
	subdirs,dirnames,subdir_numbers = globDirs(ensemble_dir,removeNature=True)
	path_to_sfs = glob(f'{subdirs[0]}*_SCALEFACTOR.nc')
	path_to_sfs.sort()
	sf_names = [pts.split('/')[-1] for pts in path_to_sfs]
	if return_not_write:
		to_return = {}
	for name in sf_names:
		pathstocombine = [path+name for path in subdirs]
		ds_files = []
		for path in pathstocombine:
			ds_files.append(xr.open_dataset(path))
		ds = xr.concat(ds_files,'Ensemble')
		ds.assign_coords({'Ensemble':np.array(subdir_numbers)})
		if timeperiod is not None:
			ds = ds.sel(time=slice(timeperiod[0], timeperiod[1]))
		if return_not_write:
			to_return[name] = ds
		else:
			if flag_snapshot:
				ds.to_netcdf(output_dir+'/SNAPSHOT_'+name)
			else:
				ds.to_netcdf(output_dir+'/'+name)
	if return_not_write:
		return to_return

def combineHemcoDiag(ensemble_dir,output_dir,timeperiod=None,prefix=''):
	subdirs,dirnames,subdir_numbers = globDirs(ensemble_dir,removeNature=True,includeOutputDir=True)
	combined_ds = []
	for subdir in subdirs:
		paths = glob(f'{subdir}/HEMCO_diagnostics.*.nc')
		paths.sort()
		ds_files = []
		for path in paths:
			ds_files.append(xr.open_dataset(path))
		combined_ds.append(xr.concat(ds_files,'time'))
	ds = xr.concat(combined_ds,'Ensemble')
	ds.assign_coords({'Ensemble':np.array(subdir_numbers)})
	if timeperiod is not None:
		ds = ds.sel(time=slice(timeperiod[0], timeperiod[1]))
	ds.to_netcdf(f'{output_dir}/{prefix}combined_HEMCO_diagnostics.nc')

def combineHemcoDiagControl(control_dir,output_dir,timeperiod=None,prefix=''):
	paths = glob(f'{control_dir}/OutputDir/HEMCO_diagnostics.*.nc')
	paths.sort()
	ds_files = []
	for path in paths:
		ds_files.append(xr.open_dataset(path))
	ds = xr.concat(ds_files,'time')
	if timeperiod is not None:
		ds = ds.sel(time=slice(timeperiod[0], timeperiod[1]))
	ds.to_netcdf(f'{output_dir}/{prefix}control_HEMCO_diagnostics.nc')


def makeDatasetForDirectory(hist_dir,species_names,timeperiod=None,hourlysub = 6,subset_rule = 'SURFACE', fullpath_output_name = None):
	specconc_list = globSubDir(hist_dir,timeperiod,hourlysub)
	concstrings = [f'{getSpcConcName()}_{name}' for name in species_names]
	ds = xr.open_mfdataset(specconc_list,concat_dim='time',combine="nested",data_vars='minimal', coords='minimal', compat='override')
	ds = ds[concstrings]
	if subset_rule=='SURFACE':
		ds = ds.isel(lev=0)
	elif subset_rule=='850':
		ds = ds.isel(lev=10) #850 hPa pressure level
	elif subset_rule=='ALL':
		pass
	if fullpath_output_name:
		ds.to_netcdf(fullpath_output_name)
	return ds

def makeDatasetForEnsemble(ensemble_dir,species_names,timeperiod=None,hourlysub = 6,subset_rule = 'SURFACE',fullpath_output_name = None):
	subdirs,dirnames,subdir_numbers = globDirs(ensemble_dir,includeOutputDir=True,removeNature=True)
	array_list = []
	for subdir in subdirs:
		print(f'Processing {subdir}')
		array_list.append(makeDatasetForDirectory(subdir,species_names,timeperiod,hourlysub,subset_rule))
	ds = xr.concat(array_list,'Ensemble')
	ds.assign_coords({'Ensemble':np.array(subdir_numbers)})
	if fullpath_output_name:
		ds.to_netcdf(fullpath_output_name)
	return ds

def getArea(ensemble_dir,pp_dir):
	subdirs,_,_ = globDirs(ensemble_dir,includeOutputDir=True,removeNature=True)
	specconc_list = globSubDir(subdirs[0],None,1)
	file = specconc_list[0]
	ds = xr.open_dataset(file)
	area = ds['AREA'].values
	np.save(f'{pp_dir}/area.npy',area)

#If you want to ensure SWOOSH data handled correctly, make sure to pass 
def makeYEachAssimPeriod(path_to_bigy_subsets,assim_time,OBS_TYPE,startdate=None,enddate=None,fullpath_output_name = None):
	masterY = {}
	bigy_list = glob(f'{path_to_bigy_subsets}/*.pkl')
	bigy_list.sort()
	bigy_list = [by for by in bigy_list if 'scaling' not in by.split('/')[-1].split('.')[0]]
	timestamps = [by.split('/')[-1].split('.')[0] for by in bigy_list]
	timestamps_datetime = [datetime.strptime(timestamp, "%Y%m%d_%H%M") for timestamp in timestamps]
	timebounds = getDatesToSubsetEachYAssimPeriod(timestamps_datetime, assim_time)
	for bigy_file,timestamp,timedate in zip(bigy_list,timestamps,timestamps_datetime):
		#Don't process values during burn in period
		if startdate is not None:
			if (timedate<startdate):
				continue
		if enddate is not None:
			if (timedate>enddate): 
				continue
		print(f'Processing the Y dictionary for time {timestamp}')
		with open(bigy_file,'rb') as f:
			bigy=pickle.load(f)
		#subset bigy to right times; needed for RIP (but not SWOOSH)
		bigy_start,bigy_end = timebounds[timedate]
		for spec in list(bigy.keys()):
			if OBS_TYPE[spec]!='SWOOSH':
				df = bigy[spec]
				df = df[(df['time']>=startdate) & (df['time']<enddate)]
				bigy[spec] = df
		masterY[timestamp] = bigy
	if fullpath_output_name:
		f = open(fullpath_output_name,"wb")
		pickle.dump(masterY,f)
		f.close()
	return masterY

def getDatesToSubsetEachYAssimPeriod(timestamps, assim_time):
	assim_time = timedelta(hours=assim_time)
	timebounds = {} #dictionary, where each timestamp has the start and end datetime for processing
	for ind,timestamp in enumerate(timestamps):
		if ind == 0:
			start = timestamp-assim_time
			timebounds[timestamp] = [start,timestamp]
		else:
			prev_time = timestamps[ind-1]
			timebounds[timestamp] = [prev_time,timestamp]
	return timebounds


def plotSurfaceCell(ds,species_name,latind,lonind,outfile=None,unit='ppt',includesNature=False):
	if unit=='ppm':
		multiplier = 1e6
	elif unit=='ppb':
		multiplier = 1e9
	elif unit=='ppt':
		multiplier = 1e12
	else:
		raise ValueError('Unit not recognized.')
	da = ds[f'{getSpcConcName()}_{species_name}'].isel(lat=latind,lon=lonind)
	time = np.array(ds['time'])
	if includesNature:
		ens = da[1::,:]*multiplier
		nature=da[0,:]*multiplier
	else:
		ens = da*multiplier
		nature=None
	ensmean = ens.mean(axis=0)
	enssd = ens.std(axis=0)
	tsPlot(time,ensmean,enssd,species_name,unit,nature,outfile=outfile)

def plotSurfaceMean(ds,species_name,outfile=None,unit='ppt',includesNature=False):
	if unit=='ppm':
		multiplier = 1e6
	elif unit=='ppb':
		multiplier = 1e9
	elif unit=='ppt':
		multiplier = 1e12
	else:
		raise ValueError('Unit not recognized.')
	da = ds[f'{getSpcConcName()}_{species_name}'].mean(axis=(2,3))
	time = np.array(ds['time'])
	if includesNature:
		ens = da[1::,:]*multiplier
		nature=da[0,:]*multiplier
	else:
		ens = da*multiplier
		nature=None
	ensmean = ens.mean(axis=0)
	enssd = ens.std(axis=0)
	tsPlot(time,ensmean,enssd,species_name,unit,nature,outfile=outfile)

def tsPlotTotalEmissions(ds_ensemble,ds_prior,area,collectionName,useLognormal = False, timeslice=None,outfile=None,conversion_factor=None):
	if timeslice is not None:
		ds_ensemble = ds_ensemble.sel(time=slice(timeslice[0],timeslice[1]))
		ds_prior = ds_prior.sel(time=slice(timeslice[0],timeslice[1]))
	if len(np.shape(ds_ensemble[collectionName])) == 5:
		axis_to_average = (2,3,4) #Emit into higher levels (e.g. aircraft)
		prior_axis_to_average = (1,2,3)
	else:
		axis_to_average = (2,3) #Surface emissions only
		prior_axis_to_average = (1,2)
	#Remove the area unit (e.g. kg/m2/s) by multiplying by area (m2)
	da = ds_ensemble[collectionName]*area
	da = da.sum(axis=axis_to_average) #sum up all emissions
	da_prior = ds_prior[collectionName]*area 
	da_prior = da_prior.sum(axis=prior_axis_to_average)
	enstime = np.array(ds_ensemble['time'])
	priortime = np.array(ds_prior['time'])
	if useLognormal:
		ensmean = np.exp(np.mean(np.log(da),axis=0))
		enssd = np.exp(np.std(np.log(da),axis=0))
	else:
		ensmean = da.mean(axis=0)
		enssd = da.std(axis=0)
	if conversion_factor is not None:
		ensmean*=conversion_factor #kg/s to Tg/d (8.64e-05 in this case) or whatever you'd like 
		enssd*=conversion_factor
		da_prior*=conversion_factor
	tsPlot(enstime,ensmean,enssd,collectionName,'kg/s',priortime=priortime,prior=da_prior,priorcolor='r',outfile=outfile)

def tsPlotSatCompare(bigY,species,numens,unit='ppb',observer_name='Observations',useControl=False,outfile=None):
	ensmeans = []
	ensstds = []
	obsmeans = []
	if useControl:
		ctrlmeans = []
	datestrs = list(bigY.keys())
	datevals = [datetime.strptime(dateval,'%Y%m%d_%H%M') for dateval in datestrs]
	for date in datestrs:
		conc2D=np.array(bigY[date][species].iloc[:,0:numens])
		if len(conc2D)==0: #In case we have no observations for a given date.
			ensmean=np.nan
			enssd=np.nan
			obsmean=np.nan
			if useControl:
				ctrlmeans.append(np.nan)
		else:
			assimperiodensmean = np.mean(conc2D,axis=0) #One average for each ensemble member
			ensmean = np.mean(assimperiodensmean) #Ensemble mean for total average
			enssd = np.std(assimperiodensmean)
			obscol=np.array(bigY[date][species]['Observations'])
			obsmean = np.mean(obscol)
			if useControl:
				ctrlcol=np.array(bigY[date][species]['Control'])
				ctrlmean = np.mean(ctrlcol)
				ctrlmeans.append(ctrlmean)
		ensmeans.append(ensmean)
		ensstds.append(enssd)
		obsmeans.append(obsmean)
	ensmeans = np.array(ensmeans)
	ensstds = np.array(ensstds)
	obsmeans = np.array(obsmeans)
	plt.rcParams.update({'font.size': 16})
	plt.figure(figsize=(6,4))
	plt.plot(datevals,ensmeans,color='b',label='Ensemble mean')
	plt.plot(datevals,ensmeans+ensstds,':',color='b')
	plt.plot(datevals,ensmeans-ensstds,':',color='b')
	plt.plot(datevals,obsmeans,color='g',label=observer_name)
	if useControl:
		ctrlmeans = np.array(ctrlmeans)
		plt.plot(datevals,ctrlmeans,color='r',label='Control')
	plt.legend()
	plt.xlabel('Time')
	plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
	plt.gca().xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=15))
	plt.gca().xaxis.set_major_formatter(NullFormatter())
	plt.gca().xaxis.set_minor_formatter(mdates.DateFormatter('%b'))
	plt.ylabel(f'{species} ({unit})')
	plt.gcf().autofmt_xdate()
	plt.gcf().tight_layout()
	if outfile:
		plt.savefig(outfile)
	else:
		plt.show()


def tsPlot(time,ensmean,enssd,species_name,unit,nature=None,priortime=None,prior=None,priorcolor='g',outfile=None):
	plt.rcParams.update({'font.size': 16})
	plt.figure(figsize=(6,4))
	plt.plot(time,ensmean,color='b',label='Ensemble')
	plt.plot(time,ensmean+enssd,':',color='b')
	plt.plot(time,ensmean-enssd,':',color='b')
	if nature is not None:
		plt.plot(time,nature,color='g',label='Nature')
		plt.legend()
	if prior is not None:
		plt.plot(priortime,prior,color=priorcolor,label='Prior')
		plt.legend()
	plt.xlabel('Time')
	plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
	plt.gca().xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=15))
	plt.gca().xaxis.set_major_formatter(NullFormatter())
	plt.gca().xaxis.set_minor_formatter(mdates.DateFormatter('%b'))
	plt.ylabel(f'{species_name} ({unit})')
	plt.gcf().autofmt_xdate()
	plt.gcf().tight_layout()
	if outfile:
		plt.savefig(outfile)
	else:
		plt.show()

#Process BigY as output by the Assimilator so that it is either (1) gridded and ready for plotting/analysis, or (2) combined by sensor.
def makeBigYArrays(bigy,gclat,gclon,nEnsemble,av_to_grid,observers_to_plot_as_point,OBS_TYPE,extra_obsdata_fields=None,useControl=False):
	dates = list(bigy.keys())
	specieslist = list(bigy[dates[0]].keys())
	to_return = {}
	#Make dictionary to return; if averaged to GC grid, we are gridding for maps.
	#If we aren't averaging to GC grid for DA, we do that here for plotting.
	#However, for surface observations, we might treat these differently, and instead aggregate by id.
	#The id used to aggregate is supplied as the value in the key/value pair OBSERVERS_TO_PLOT_AS_POINTS.
	#Finally, for SWOOSH data, we average to vertical profiles
	for species in specieslist:
		to_return[species] = {}
		is_Av = av_to_grid[species] == 'True'
		if extra_obsdata_fields is not None:
			bonus_fields = extra_obsdata_fields[species]
		else:
			bonus_fields = []
		#If plotting as points, dictionary will have keys for location codes (given as values in observers_to_plot_as_points dictionary)
		#Otherwise, we will be generating maps, and will prep them with empty arrys
		#Unless we have SWOOSH data, in which case we will be doing vertical profiles
		if OBS_TYPE[species] == "SWOOSH":
			to_return[species]['interpret_as'] = 'profile'
		if species in observers_to_plot_as_points:
			to_return[species]['interpret_as'] = 'points'
		else:
			to_return[species]['interpret_as'] = 'map'
			to_return[species]['sim_obs'] = np.zeros([len(dates),len(gclat),len(gclon)])*np.nan
			to_return[species]['obs'] = np.zeros([len(dates),len(gclat),len(gclon)])*np.nan
			to_return[species]['dates'] = dates
			to_return[species]['obscount'] = np.zeros([len(dates),len(gclat),len(gclon)])
			if is_Av:
				to_return[species]['obscount_avg'] = np.zeros([len(dates),len(gclat),len(gclon)])
			for field in bonus_fields:
				to_return[species][field] = np.zeros([len(dates),len(gclat),len(gclon)])*np.nan
			if useControl:
				to_return[species]['control'] = np.zeros([len(dates),len(gclat),len(gclon)])*np.nan
		#Go ahead and do the regridding or aggregation to point data.
		for date_index, date in enumerate(dates): 
			ydict = bigy[date][species]
			#Regrid for maps
			if to_return[species]['interpret_as'] == 'map':
				latdict = np.array(ydict['Latitude'])
				londict = np.array(ydict['Longitude'])
				trueobsdict = np.array(ydict['Observations'])
				simobsdict = np.mean(np.array(ydict.iloc[:,0:nEnsemble]),axis=1) #Get the ensemble sim obs values, average to get sim obs dict
				if is_Av:
					countdict = np.array(ydict['Num_Averaged'])
				if useControl:
					controldict = np.array(ydict['Control'])
				#Aggregate to map
				latind = np.abs(latdict.reshape((1, -1)) - gclat.reshape((-1, 1)))
				latind = latind.argmin(axis=0)
				lonind = np.abs(londict.reshape((1, -1)) - gclon.reshape((-1, 1)))
				lonind = lonind.argmin(axis=0)
				pairedind = ((latind+1)*10000)+(lonind+1)
				uniqueind,countind = np.unique(pairedind,return_counts=True)
				uniquelonind = (uniqueind % 10000)-1
				uniquelatind = np.floor(uniqueind / 10000).astype(int)-1
				#Number of matches is interpreted differently if we have already averaged to the GC grid in assimilation
				if is_Av: 
					to_return[species]['obscount_avg'][date_index,uniquelatind,uniquelonind]=countind
				else:
					to_return[species]['obscount'][date_index,uniquelatind,uniquelonind]=countind
				#Aggregate everything at matching latvals and lonvals
				for lonindval,latindval in zip(uniquelonind,uniquelatind):
					dictind = np.where((latdict==gclat[latindval])&(londict==gclon[lonindval]))[0]
					if is_Av:
						totalcount = np.sum(countdict[dictind])
						to_return[species]['obscount'][date_index,latindval,lonindval]=totalcount
					to_return[species]['obs'][date_index,latindval,lonindval]=np.mean(trueobsdict[dictind])
					to_return[species]['sim_obs'][date_index,latindval,lonindval]=np.mean(simobsdict[dictind])
					if useControl:
						to_return[species]['control'][date_index,latindval,lonindval]=np.mean(controldict[dictind])
					for field in bonus_fields:
						to_return[species][field][date_index,latindval,lonindval] = np.mean(np.array(ydict[field])[dictind])
			#Aggregate for points
			elif to_return[species]['interpret_as'] == 'points':
				agg_dict = np.array(ydict[observers_to_plot_as_points[species]]) #Get the point codes for aggregation
				uniqueind,countind = np.unique(agg_dict,return_counts=True)
				for station in uniqueind:
					dictind = np.where(agg_dict==station)[0]
					#If we've not seen this station before, create new entries
					if station not in to_return[species]:
						to_return[species][station] = {}
						to_return[species][station]['obs'] = []
						to_return[species][station]['sim_obs'] = []
						to_return[species][station]['lat'] = []
						to_return[species][station]['lon'] = []
						to_return[species][station]['time'] = []
						if useControl:
							to_return[species][station]['control'] = []
						for field in bonus_fields:
							to_return[species][station][field] = []
					#Append to existing lists, for concatenation later
					to_return[species][station]['obs'].append(np.array(ydict['Observations'])[dictind])
					to_return[species][station]['sim_obs'].append(np.mean(np.array(ydict.iloc[:,0:nEnsemble]),axis=1)[dictind])
					to_return[species][station]['lat'].append(np.array(ydict['Latitude'])[dictind])
					to_return[species][station]['lon'].append(np.array(ydict['Longitude'])[dictind])
					to_return[species][station]['time'].append(np.array(ydict['time'])[dictind])
					if useControl:
						to_return[species][station]['control'].append(np.array(ydict['Control'])[dictind])
					for field in bonus_fields:
						to_return[species][station][field].append(np.array(ydict[field])[dictind])
		#If this is a map, we're done. If this is a point, we need to concatenate and sort by date.
		if to_return[species]['interpret_as'] == 'points':
			for station in to_return[species]:
				if station == 'interpret_as': #this is the one entry that is not a station in this dictionary
					continue
				else:
					for field in to_return[species][station]: #Go through each field and concatenate
						to_return[species][station][field] = np.concatenate(to_return[species][station][field])
					sorted_time_inds = to_return[species][station]['time'].argsort() #Get sorted indices by time. 
					for field in to_return[species][station]: #Go through each field and sort in time order.
						to_return[species][station][field] = to_return[species][station][field][sorted_time_inds]
	#Done!
	return to_return

