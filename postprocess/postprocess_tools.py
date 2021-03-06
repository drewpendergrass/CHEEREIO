import numpy as np 
import xarray as xr
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter
from datetime import datetime
import sys
import pickle
import pandas as pd
sys.path.append('../core')
from HIST_Ens import HIST_Ens 
import tropomi_tools as tt 

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

def combineScaleFactors(ensemble_dir,output_dir):
	subdirs,dirnames,subdir_numbers = globDirs(ensemble_dir,removeNature=True)
	path_to_sfs = glob(f'{subdirs[0]}*_SCALEFACTOR.nc')
	path_to_sfs.sort()
	sf_names = [pts.split('/')[-1] for pts in path_to_sfs]
	files_split_by_sf = []
	for name in sf_names:
		pathstocombine = [path+name for path in subdirs]
		ds_files = []
		for path in pathstocombine:
			ds_files.append(xr.open_dataset(path))
		ds = xr.concat(ds_files,'Ensemble')
		ds.assign_coords({'Ensemble':np.array(subdir_numbers)})
		ds.to_netcdf(output_dir+'/'+name)

def makeDatasetForDirectory(hist_dir,species_names,timeperiod=None,hourlysub = 6,subset_rule = 'SURFACE', fullpath_output_name = None):
	specconc_list = globSubDir(hist_dir,timeperiod,hourlysub)
	concstrings = [f'SpeciesConc_{name}' for name in species_names]
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
	subdirs,dirnames,subdir_numbers = globDirs(ensemble_dir,includeOutputDir=True)
	array_list = []
	for subdir in subdirs:
		print(f'Processing {subdir}')
		array_list.append(makeDatasetForDirectory(subdir,species_names,timeperiod,hourlysub,subset_rule))
	ds = xr.concat(array_list,'Ensemble')
	ds.assign_coords({'Ensemble':np.array(subdir_numbers)})
	if fullpath_output_name:
		ds.to_netcdf(fullpath_output_name)
	return ds

def makeDatasetForDirectoryLevelEdge(hist_dir,timeperiod=None,hourlysub = 6, fullpath_output_name = None, average_levels = True):
	edgeconc_list = globSubDirLevelEdge(hist_dir,timeperiod,hourlysub)
	concstring = 'Met_PEDGE'
	ds = xr.open_mfdataset(edgeconc_list,concat_dim='time',combine="nested",data_vars='minimal', coords='minimal', compat='override')
	ds = ds[concstring]
	if fullpath_output_name:
		ds.to_netcdf(fullpath_output_name)
	return ds

def makeDatasetForEnsembleLevelEdge(ensemble_dir,timeperiod=None,hourlysub = 6,fullpath_output_name = None):
	subdirs,dirnames,subdir_numbers = globDirs(ensemble_dir,includeOutputDir=True)
	array_list = []
	for subdir in subdirs:
		print(f'Processing {subdir}')
		array_list.append(makeDatasetForDirectoryLevelEdge(subdir,timeperiod,hourlysub))
	ds = xr.concat(array_list,'Ensemble')
	ds.assign_coords({'Ensemble':np.array(subdir_numbers)})
	if fullpath_output_name:
		ds.to_netcdf(fullpath_output_name)
	return ds

def makeYEachAssimPeriod(timestamp_list, useLevelEdge=False,useStateMet = False,useArea=False, use_numav = True, use_albedo=True, fullpath_output_name = None):
	masterY = {}
	for timestamp in timestamp_list:
		print(f'Processing the Y dictionary for time {timestamp}')
		hist = HIST_Ens(timestamp=timestamp,useLevelEdge=useLevelEdge,useStateMet = useStateMet,useArea=useArea,testing=False,saveAlbedo=use_albedo)
		bigy = hist.bigYDict 
		for spec in list(bigy.keys()):
			t = [np.datetime64(int(tt),'ns') for tt in bigy[spec][4]]
			t = np.array(t,dtype='datetime64[us]')
			colnum = np.shape(bigy[spec][0])[1]
			colnames = []
			for i in range(colnum):
				colnames.append(f"Ens{str(i+1).zfill(3)}")
			df = pd.DataFrame(bigy[spec][0], columns = colnames)
			df['Satellite'] = bigy[spec][1]
			df['Latitude'] = bigy[spec][2]
			df['Longitude'] = bigy[spec][3]
			if use_numav:
				df['Num_Averaged'] = bigy[spec][5]
				if use_albedo:
					df['Albedo_SWIR'] = bigy[spec][6]
					df['Albedo_NIR'] = bigy[spec][7]
					df['Blended_Albedo'] = bigy[spec][8]
			else:
				df['Num_Averaged'] = None
				if use_albedo:
					df['Albedo_SWIR'] = bigy[spec][5]
					df['Albedo_NIR'] = bigy[spec][6]
					df['Blended_Albedo'] = bigy[spec][7]
			df['time'] = t
			bigy[spec] = df
		masterY[timestamp] = bigy
	if fullpath_output_name:
		f = open(fullpath_output_name,"wb")
		pickle.dump(masterY,f)
		f.close()
	return masterY


def makeYWholePeriod(timestamp,hourlysub=6,use_numav = False, fullpath_output_name = None):
	hist = HIST_Ens(timestamp=timestamp,useLevelEdge=True,fullperiod=True,interval=hourlysub,testing=False)
	bigy = hist.bigYDict 
	for spec in list(bigy.keys()):
		t = [datetime.strptime(tt,"%Y-%m-%dT%H:%M:%S.%fZ") for tt in bigy[spec][4]]
		t = np.array(t,dtype='datetime64[us]')
		colnum = np.shape(bigy[spec][0])[1]
		colnames = []
		for i in range(colnum):
			colnames.append(f"Ens{str(i+1).zfill(3)}")
		df = pd.DataFrame(bigy[spec][0], columns = colnames)
		df['Satellite'] = bigy[spec][1]
		df['Latitude'] = bigy[spec][2]
		df['Longitude'] = bigy[spec][3]
		if use_numav:
			df['Num_Averaged'] = bigy[spec][4]
		else:
			df['Num_Averaged'] = np.ones(len(bigy[spec][3]))
		df['Time'] = t
		bigy[spec] = df
	if fullpath_output_name:
		f = open(fullpath_output_name,"wb")
		pickle.dump(bigy,f)
		f.close()
	return bigy


def plotSurfaceCellEnsMeanNorm(ds,species_name,latind,lonind,outfile=None,unit='ppm'):
	if unit=='ppm':
		multiplier = 1e6
	elif unit=='ppb':
		multiplier = 1e9
	elif unit=='ppt':
		multiplier = 1e12
	else:
		raise ValueError('Unit not recognized.')
	da = np.array(ds[f'SpeciesConc_{species_name}'])
	time = np.array(ds['time'])
	ens = da[:,:,latind,lonind]*multiplier
	ensmean = np.mean(ens,axis=0)
	enssd = np.std(ens,axis=0)
	tsPlot(time,ensmean-ensmean,enssd,species_name,unit,outfile=outfile)

def plotSurfaceCell(ds,species_name,latind,lonind,outfile=None,unit='ppt',includesNature=False):
	if unit=='ppm':
		multiplier = 1e6
	elif unit=='ppb':
		multiplier = 1e9
	elif unit=='ppt':
		multiplier = 1e12
	else:
		raise ValueError('Unit not recognized.')
	da = ds[f'SpeciesConc_{species_name}'].isel(lat=latind,lon=lonind)
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
	da = ds[f'SpeciesConc_{species_name}'].mean(axis=(2,3))
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

def tsPlotSatCompare(bigY,species,numens,unit='ppb',satellite_name='TROPOMI',outfile=None):
	ensmeans = []
	ensstds = []
	satmeans = []
	datestrs = list(bigY.keys())
	datevals = [datetime.strptime(dateval,'%Y%m%d_%H%M') for dateval in datestrs]
	for date in datestrs:
		conc2D=np.array(bigY[date][species].iloc[:,1:(numens+1)])
		assimperiodensmean = np.mean(conc2D,axis=0) #One average for each ensemble member
		ensmean = np.mean(assimperiodensmean) #Ensemble mean for total average
		enssd = np.std(assimperiodensmean)
		satcol=np.array(bigY[date][species]['Satellite'])
		satmean = np.mean(satcol)
		ensmeans.append(ensmean)
		ensstds.append(enssd)
		satmeans.append(satmean)
	ensmeans = np.array(ensmeans)
	ensstds = np.array(ensstds)
	satmeans = np.array(satmeans)
	plt.rcParams.update({'font.size': 16})
	plt.figure(figsize=(6,4))
	plt.plot(datevals,ensmeans,color='b',label='Ensemble mean')
	plt.plot(datevals,ensmeans+ensstds,':',color='b')
	plt.plot(datevals,ensmeans-ensstds,':',color='b')
	plt.plot(datevals,satmeans,color='g',label=satellite_name)
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


def tsPlotSatCompareFullRange(df,species,numens,freq='H',unit='ppb',satellite_name='TROPOMI',outfile=None):
	df = df.groupby(pd.Grouper(key='Time',freq=freq)).mean()
	df.reset_index(inplace=True)
	conc2D=np.array(df.iloc[:,1:(numens+1)])
	satcol=np.array(df['Satellite'])
	sattime=np.array(df['Time'])
	ensmean = np.mean(conc2D,axis=1)
	enssd = np.std(conc2D,axis=1)
	plt.rcParams.update({'font.size': 16})
	plt.figure(figsize=(6,4))
	plt.plot(sattime,ensmean,color='b',label='Ensemble mean')
	plt.plot(sattime,ensmean+enssd,':',color='b')
	plt.plot(sattime,ensmean-enssd,':',color='b')
	plt.plot(sattime,satcol,color='g',label=satellite_name)
	plt.legend()
	plt.xlabel('Time')
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
	plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=5))
	plt.ylabel(f'{species} ({unit})')
	plt.gcf().autofmt_xdate()
	plt.gcf().tight_layout()
	if outfile:
		plt.savefig(outfile)
	else:
		plt.show()



def tsPlot(time,ensmean,enssd,species_name,unit,nature=None,outfile=None):
	plt.rcParams.update({'font.size': 16})
	plt.figure(figsize=(6,4))
	plt.plot(time,ensmean,color='b',label='Ensemble mean')
	plt.plot(time,ensmean+enssd,':',color='b')
	plt.plot(time,ensmean-enssd,':',color='b')
	if nature is not None:
		plt.plot(time,nature,color='g',label='Nature')
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

def emisPlot(time,ensmean,enssd,name,outfile=None):
	plt.rcParams.update({'font.size': 16})
	plt.figure(figsize=(6,4))
	plt.plot(time,ensmean,color='b')
	plt.plot(time,ensmean+enssd,':',color='b')
	plt.plot(time,ensmean-enssd,':',color='b')
	plt.xlabel('Time')
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
	plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=3))
	plt.ylabel(f'{name}')
	plt.gcf().autofmt_xdate()
	plt.gcf().tight_layout()
	if outfile:
		plt.savefig(outfile)
	else:
		plt.show()

def plotEmissionsCell(ds_file,latind,lonind,outfile=None):
	ds = xr.open_dataset(ds_file)
	time = np.array(ds['time'])
	da = np.array(ds['Scalar'])
	ens = da[:,:,latind,lonind]
	ensmean = np.mean(ens,axis=0)
	enssd = np.std(ens,axis=0)
	emisPlot(time,ensmean,enssd,'_'.join(ds_file.split('/')[-1].split('_')[0:-1]),outfile)

