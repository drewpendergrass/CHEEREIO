import numpy as np
import pandas as pd
import xarray as xr
from datetime import date
import toolbox as tx
import settings_interface as si 
from glob import glob
import sys

teststr = str(sys.argv[1])
if (teststr=="TESTING") or (teststr == "TESTPROD"):
	testbool = True
else:
	testbool = False

spc_config = si.getSpeciesConfig()

parent_dir = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/ensemble_runs"
subdirs = glob(f"{parent_dir}/*/")
subdirs.remove(f"{parent_dir}/logs/")
dirnames = [d.split('/')[-2] for d in subdirs]
subdir_numstring = [n.split('_')[-1] for n in dirnames]
subdir_nums = [int(n.split('_')[-1]) for n in dirnames]

maxdir = max(subdir_nums)
meanval = (maxdir+1)/2

#Load in the relevant parameters from ens_config.json.

emis_scaling_factors = spc_config['CONTROL_VECTOR_EMIS'].keys()
mask_ocean_bool = spc_config['MaskOceanScaleFactor']
mask_coast_bool = spc_config['MaskCoastsGT25pctOcean']
mask_arctic_bool = spc_config['Mask60NScaleFactor']
mask_antarctic_bool = spc_config['Mask60SScaleFactor']
perttype = spc_config["pertType"]
perturbation = spc_config["pPERT"]
minsf = spc_config['MinimumScalingFactorAllowed']
maxsf = spc_config['MaximumScalingFactorAllowed']
correlatedInitialScalings = spc_config['correlatedInitialScalings']
speedyCorrelationApprox = spc_config['speedyCorrelationApprox'] == 'True'
corrDistances = spc_config['corrDistances']


timestamp = str(sys.argv[2]) #Time for scaling factor time dimension. Format assumed to be YYYYMMDD
timestamp = timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]

if len(spc_config["REGION"])==0:
	gridlabel = spc_config["RES"]
else:
	gridlabel = f'{spc_config["REGION"]}_{spc_config["met_name"]}'

lon,lat,mask = tx.makeLatLonGridWithMask(gridlabel)

#Check that the user-inputed perturbations are sensible given user-supplied interpretation rule
for emis in emis_scaling_factors:
	pt = perttype[emis]
	p = float(perturbation[emis])
	corrbool = correlatedInitialScalings[emis]
	if (corrbool == 'True') and (pt != "std"):
		pt = "std"
		print(f'WARNING: Correlated initial scalings require the "std" setting. Overriding your setting of {pt} for emission {emis}.')
	if pt == "exp":
		if (p <= 1): #perturbation, max positive amount. i.e. if it is 4 scaling factors will range between 0.25 and 4. Uniform distribution used, no correlation.
			raise ValueError('Exponential perturbation must be at least 1.')
	elif pt == "percent":
		if (p <= 0): #perturbation, fraction. i.e. if it is 0.5 scaling factors will range between 0.5 and 1.5. Uniform distribution used, no correlation.
			raise ValueError('Percent perturbation must be positive.')
		elif (p > 1): #perturbation, fraction. i.e. if it is 0.5 scaling factors will range between 0.5 and 1.5. Uniform distribution used, no correlation.
			raise ValueError('Percent perturbation must be 1 or less.')
	elif pt == "std":
		if (p <= 0): #perturbation, standard deviation from a normal distribution. i.e. if it is 0.5 scaling factors will be centered at 1 with standard deviation 0.5.
			raise ValueError('Standard deviation perturbation must be positive.')
	else:
		raise ValueError("Perturbation type unrecognized.")

scaling_factor_cube = np.zeros((len(subdir_numstring), len(emis_scaling_factors),len(lat),len(lon)))
subdircount = 0
speciescount = 0

if ('True' in list(correlatedInitialScalings.values())) and (not speedyCorrelationApprox):
	distmat = tx.getDistMat(gridlabel)

for stringnum,num in zip(subdir_numstring,subdir_nums): #Loop through the non-nature directories
	if num == 0:
		continue
	for emis_name in emis_scaling_factors: #Loop through the species we want scaling factors for
		maskoceanboolval=mask_ocean_bool[emis_name]
		maskarcticboolval=mask_arctic_bool[emis_name]
		maskantarcticboolval=mask_antarctic_bool[emis_name]
		pt=perttype[emis_name]
		p=float(perturbation[emis_name])
		minval=float(minsf[emis_name])
		maxval=float(maxsf[emis_name])
		corrbool=correlatedInitialScalings[emis_name]
		corrdist=float(corrDistances[emis_name])
		#Generate random uniform scaling factors. If testing, just generate uniform field of same percentage below/above mean as restarts, offset by configurable parameter
		if testbool:
			offset = 1
			scale = 0
			scaling_factors = (scale*np.random.rand(1,len(lat),len(lon)))+offset
			scaling_factors *= ((num/meanval)+float(spc_config['TESTBIAS']))
		else:
			if corrbool == "True": #Will sample a normal with correlation
				if speedyCorrelationApprox:
					scaling_factors = tx.speedySample(corrdist,lat[10]-lat[9],p, (len(lat),len(lon)))
				else:
					cov = tx.makeCovMat(distmat,corrdist)
					scaling_factors = tx.sampleCorrelatedStructure(corrdist,cov,p, (len(lat),len(lon)))
			else:
				if pt == "exp":
					scaling_factor_exp = (2*np.random.rand(len(lat),len(lon)))-1
					scaling_factors = p**scaling_factor_exp
				elif pt == "percent":
					scaling_factors = (2*p*np.random.rand(len(lat),len(lon)))-p+1
				elif pt == "std":
					scaling_factors = np.random.normal(loc=1,scale=p,size=[len(lat),len(lon)])
		if ~np.isnan(minval): #Enforce minimum sf.
			scaling_factors[scaling_factors<minval] = minval
		if ~np.isnan(maxval): #Enforce maximum sf.
			scaling_factors[scaling_factors>maxval] = maxval
		if maskoceanboolval=='True':
			scaling_factors[mask[1],mask[0]] = 1
		if maskarcticboolval=='True':
			latwhere = np.where(lat>=60)[0]
			scaling_factors[latwhere,:] = 1
		if maskantarcticboolval=='True':
			latwhere = np.where(lat<=-60)[0]
			scaling_factors[latwhere,:] = 1
		scaling_factor_cube[subdircount,speciescount,:,:] = scaling_factors[:,:]
		speciescount+=1
	subdircount+=1
	speciescount = 0

#Utility function to save NetCDF
def saveNetCDF(scalefactor2D, stringnum, emis_name):
	scaling_factors = np.expand_dims(scalefactor2D, axis=0) # add time dimension
	name = f'{emis_name}_SCALEFACTOR'
	outdir = f"{parent_dir}/{spc_config['RUN_NAME']}_{stringnum}"
	endtime = spc_config['END_DATE']
	endtime = endtime[0:4]+'-'+endtime[4:6]+'-'+endtime[6:8]
	ds = xr.Dataset(
		{"Scalar": (("time","lat","lon"), scaling_factors,{"long_name": "Scaling factor", "units":"1"})},
		coords={
			"time": (["time"], np.array([0.0]), {"long_name": "time", "calendar": "standard", "units":f"hours since {timestamp} 00:00:00"}),
			"lat": (["lat"], lat,{"long_name": "Latitude", "units":"degrees_north"}),
			"lon": (["lon"], lon,{"long_name": "Longitude", "units":"degrees_east"})
		},
		attrs={
			"Title":"CHEEREIO scaling factors",
			"Conventions":"COARDS",
			"Format":"NetCDF-4",
			"Model":"GENERIC",
			"NLayers":"1",
			"History":f"Originally auto-generated by the CHEEREIO utility on {str(date.today())}",
			"Start_Date":f"{timestamp}",
			"Start_Time":"0",
			"End_Date":f"{endtime}",
			"End_Time":"0"
		}
	)
	ds.to_netcdf(f"{outdir}/{name}.nc")
	print(f"Scaling factors \'{name}.nc\' in folder {spc_config['RUN_NAME']}_{stringnum} initialized successfully!")

#subtract mean to avoid biased initial conditions (e.g. mean one), save out initial std for entire ensemble, save out scalefactors for each ensemble member
for ind,emis_name in enumerate(emis_scaling_factors):
	scaling_factor_mean = np.mean(scaling_factor_cube[:,ind,:,:],axis=0)
	scaling_factor_cube[:,ind,:,:] -= (scaling_factor_mean-1) #transform to mean of 1
	scaling_factor_sd = np.std(scaling_factor_cube[:,ind,:,:],axis=0)
	np.save(f'{parent_dir}/{emis_name}_SCALEFACTOR_INIT_STD.npy',scaling_factor_sd) #save out initial std for entire ensemble
	for dirind, stringnum in enumerate(subdir_numstring):
		saveNetCDF(scaling_factor_cube[dirind,ind,:,:],stringnum, emis_name) #save out scalefactors for each ensemble member


