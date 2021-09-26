import numpy as np 
import xarray as xr
from glob import glob
import matplotlib.pyplot as plt

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

def makeDatasetForDirectory(hist_dir,species_names,fullpath_output_name = None):
	specconc_list = glob(f'{hist_dir}/GEOSChem.SpeciesConc*.nc4')
	specconc_list.sort()
	concstrings = [f'SpeciesConc_{name}' for name in species_names]
	ds = xr.open_mfdataset(specconc_list,concat_dim='time',combine="nested",data_vars='minimal', coords='minimal', compat='override')
	ds = ds[concstrings]
	if fullpath_output_name:
		ds.to_netcdf(fullpath_output_name)
	return ds

def makeDatasetForEnsemble(ensemble_dir,species_names,fullpath_output_name = None):
	subdirs,dirnames,subdir_numbers = globDirs(ensemble_dir,includeOutputDir=True)
	array_list = []
	for subdir in subdirs:
		print(f'Processing {subdir}')
		array_list.append(makeDatasetForDirectory(subdir,species_names))
	ds = xr.concat(array_list,'Ensemble')
	ds.assign_coords({'Ensemble':np.array(subdir_numbers)})
	if fullpath_output_name:
		ds.to_netcdf(fullpath_output_name)
	return ds

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
	ens = da[:,:,0,latind,lonind]*multiplier
	ensmean = np.mean(ens,axis=0)
	enssd = np.std(ens,axis=0)
	tsPlot(time,ensmean-ensmean,enssd,species_name,unit,outfile=outfile)

def plotSurfaceCell(ds,species_name,latind,lonind,outfile=None,unit='ppm',includesNature=False,nature_error=None,natureErrType='relative'):
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
	if includesNature:
		nature = da[0,:,0,latind,lonind]*multiplier
		if natureErrType=="relative":
			naterr = nature*nature_error
		elif natureErrType=="absolute":
			naterr = np.repeat(nature_error*multiplier,len(nature))
		else:
			raise ValueError('Nature error must be relative or absolute.')
		ens = da[1::,:,0,latind,lonind]*multiplier
		ensmean = np.mean(ens,axis=0)
		enssd = np.std(ens,axis=0)
		tsPlot(time,ensmean,enssd,species_name,unit,nature,naterr,outfile=outfile)
	else:
		ens = da[:,:,0,latind,lonind]*multiplier
		ensmean = np.mean(ens,axis=0)
		enssd = np.std(ens,axis=0)
		tsPlot(time,ensmean,enssd,species_name,unit,outfile=outfile)

def tsPlot(time,ensmean,enssd,species_name,unit,nature=None,naterr=None,outfile=None):
	plt.figure(figsize=(10,9))
	if nature:
		plt.plot(time,nature,color='g',label='Nature run')
		plt.plot(time,nature+naterr,':',color='g')
		plt.plot(time,nature-naterr,':',color='g')
		plt.plot(time,ensmean,color='b',label='Ensemble runs')
		plt.plot(time,ensmean+enssd,':',color='b')
		plt.plot(time,ensmean-enssd,':',color='b')
		plt.xlabel('Time')
		plt.ylabel(f'{species_name} ({unit})')
		plt.legend()
	else:
		plt.plot(time,ensmean,color='b')
		plt.plot(time,ensmean+enssd,':',color='b')
		plt.plot(time,ensmean-enssd,':',color='b')
		plt.xlabel('Time')
		plt.ylabel(f'{species_name} ({unit})')
	if outfile:
		plt.savefig(outfile)

def emisPlot(time,ensmean,enssd,name,outfile=None):
	plt.figure(figsize=(10,9))
	plt.plot(time,ensmean,color='b')
	plt.plot(time,ensmean+enssd,':',color='b')
	plt.plot(time,ensmean-enssd,':',color='b')
	plt.xlabel('Time')
	plt.ylabel(f'{name}')
	if outfile:
		plt.savefig(outfile)

def plotEmissionsCell(ds_file,latind,lonind,outfile=None):
	ds = xr.open_dataset(ds_file)
	time = np.array(ds['time'])
	da = np.array(ds['Scalar'])
	ens = da[:,:,latind,lonind]
	ensmean = np.mean(ens,axis=0)
	enssd = np.std(ens,axis=0)
	emisPlot(time,ensmean,enssd,ds_file.split('/')[-1].split('_')[0:-2],outfile)
