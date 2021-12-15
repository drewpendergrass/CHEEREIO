import numpy as np
import pandas as pd
import xarray as xr
from datetime import date
import toolbox as tx
from glob import glob
import sys

teststr = str(sys.argv[1])
if teststr=="TESTING":
	testbool = True
	spc_config = tx.getSpeciesConfig(testing=True)
elif teststr=="TESTPROD":
	testbool = True
	spc_config = tx.getSpeciesConfig(testing=False) #Production run, but do deterministic scalings.
else:
	testbool = False
	spc_config = tx.getSpeciesConfig(testing=False)



parent_dir = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}/ensemble_runs"
subdirs = glob(f"{parent_dir}/*/")
subdirs.remove(f"{parent_dir}/logs/")
dirnames = [d.split('/')[-2] for d in subdirs]
subdir_numstring = [n.split('_')[-1] for n in dirnames]
subdir_nums = [int(n.split('_')[-1]) for n in dirnames]

maxdir = max(subdir_nums)
meanval = (maxdir+1)/2

emis_scaling_factors = spc_config['CONTROL_VECTOR_EMIS'].keys()
mask_ocean_bool = spc_config['MaskOceanScaleFactor']
mask_arctic_bool = spc_config['MaskArcticCircleScaleFactor']
mask_antarctic_bool = spc_config['Mask60SScaleFactor']

timestamp = str(sys.argv[2]) #Time for scaling factor time dimension. Format assumed to be YYYYMMDD
timestamp = timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]

if len(spc_config["REGION"])==0:
	gridlabel = spc_config["RES"]
else:
	gridlabel = f'{spc_config["REGION"]}_{spc_config["met_name"]}'

subhalf = False
subquarter = False
#Create traditional GEOS-Chem longitude and latitude centers, as specified by the settings in ens_config.json
if gridlabel == '4.0x5.0':
	lon = np.arange(-180.0,176.0, 5.0)
	lat = np.concatenate([[-89.0],np.arange(-86.0,87.0, 4.0), [89.0]])
	mask = pd.read_csv('../templates/landmask_4x5_gcgrid.csv',header=None)
	mask = np.array(np.transpose(mask))
	mask = np.where(mask==0)
elif gridlabel == '2.0x2.5':
	lon = np.arange(-180.0,178.0, 2.5)
	lat = np.concatenate([[-89.5],np.arange(-88.0,89.0, 2.0), [89.5]])
	mask = pd.read_csv('../templates/landmask_2x2p5_gcgrid.csv',header=None)
	mask = np.array(np.transpose(mask))
	mask = np.where(mask==0)
elif gridlabel == '1x1':
	lon = np.arange(-179.5,180.0, 1.0)
	lat = np.arange(-89.5,90, 1.0)
	mask = pd.read_csv('../templates/landmask_1x1_gcgrid.csv',header=None)
	mask = np.array(np.transpose(mask))
	mask = np.where(mask==0)
elif (gridlabel == '0.5x0.625') | (gridlabel == 'MERRA2'): #MERRA2 NATIVE GRID
	lon = -180.0 + (0.625*np.arange(0.0,576.0,1.0))
	lat = -90.0 + (0.5*np.arange(0.0,361.0,1.0))
	mask = pd.read_csv('../templates/landmask_0p5x0p625_gcgrid.csv',header=None)
	mask = np.array(np.transpose(mask))
	mask = np.where(mask==0)
elif (gridlabel == 'AS_MERRA2'): #ASIA NESTED GRID FOR MERRA2
	lon = np.arange(60.0,150.01, 0.625)
	lat = np.arange(-11.0,55.01, 0.5)
	subhalf = True
elif (gridlabel == 'EU_MERRA2'): #EUROPE NESTED GRID FOR MERRA2
	lon = np.arange(-30.0,50.01, 0.625)
	lat = np.arange(30.0,70.01, 0.5)
	subhalf = True
elif (gridlabel == 'NA_MERRA2'): #NORTH AMERICA NESTED GRID FOR MERRA2
	lon = np.arange(-140.0,-39.99, 0.625)
	lat = np.arange(10.0,70.01, 0.5)
	subhalf = True
elif (gridlabel == '0.25x0.3125') | (gridlabel == 'GEOSFP'):
	lon = -180.0 + (0.3125*np.arange(0.0,1152.0,1.0))
	lat = -90.0 + (0.25*np.arange(0.0,721.0,1.0))
	mask = pd.read_csv('../templates/landmask_0p25x0p3125_gcgrid.csv',header=None)
	mask = np.array(np.transpose(mask))
	mask = np.where(mask==0)
elif (gridlabel == 'CH_GEOSFP'): #CHINA NESTED GRID FOR GEOS-FP
	lon = np.arange(70.0,140.01, 0.3125)
	lat = np.arange(15.0,55.01, 0.25)
	subquarter = True
elif (gridlabel == 'EU_GEOSFP'): #EU NESTED GRID FOR GEOS-FP
	lon = np.arange(-15.0,40.01, 0.3125)
	lat = np.arange(32.75,61.26, 0.25)
	subquarter = True
elif (gridlabel == 'NA_GEOSFP'): #NA NESTED GRID FOR GEOS-FP
	lon = np.arange(-130.0,-59.99, 0.3125)
	lat = np.arange(9.75,60.01, 0.25)
	subquarter = True
else:
	raise ValueError('Scaling factor initialization utility does not recognize grid specification.')

if subhalf or subquarter:
	if subhalf:
		fulllon = -180.0 + (0.625*np.arange(0.0,576.0,1.0))
		fulllat = -90.0 + (0.5*np.arange(0.0,361.0,1.0))
		mask = pd.read_csv('../templates/landmask_0p5x0p625_gcgrid.csv',header=None)
	elif subquarter:
		fulllon = -180.0 + (0.3125*np.arange(0.0,1152.0,1.0))
		fulllat = -90.0 + (0.25*np.arange(0.0,721.0,1.0))
		mask = pd.read_csv('../templates/landmask_0p25x0p3125_gcgrid.csv',header=None)
	lonok = np.where((fulllon>=np.min(lon)) & (fulllon<=np.max(lon)))[0]
	latok = np.where((fulllat>=np.min(lat)) & (fulllat<=np.max(lat)))[0]
	mask = np.array(np.transpose(mask))
	mask = mask[lonok,latok]
	mask = np.where(mask==0)

perturbation = float(spc_config["pPERT"]) #perturbation, max positive amount. i.e. if it is 4 scaling factors will range between 0.25 and 4.

if (perturbation <= 1):
	raise ValueError('Perturbation must be at least 1.')

for stringnum,num in zip(subdir_numstring,subdir_nums): #Loop through the non-nature directories
	if num == 0:
		continue
	for emis_name,maskoceanboolval,maskarcticboolval,maskantarcticboolval in zip(emis_scaling_factors,mask_ocean_bool,mask_arctic_bool,mask_antarctic_bool): #Loop through the species we want scaling factors for
		#Generate random uniform scaling factors. If testing, just generate uniform field of same percentage below/above mean as restarts, offset by configurable parameter
		if testbool:
			offset = 1
			scale = 0
			scaling_factors = (scale*np.random.rand(1,len(lat),len(lon)))+offset
			scaling_factors *= ((num/meanval)+float(spc_config['TESTBIAS']))
		else:
			scaling_factor_exp = (2*np.random.rand(1,len(lat),len(lon)))-1
			scaling_factors = perturbation**scaling_factor_exp
		if maskoceanboolval=='True':
			scaling_factors[0,mask[1],mask[0]] = 1
		if maskarcticboolval=='True':
			latwhere = np.where(lat>66.55)[0]
			scaling_factors[0,latwhere,:] = 1
		if maskantarcticboolval=='True':
			latwhere = np.where(lat<=-60)[0]
			scaling_factors[0,latwhere,:] = 1
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
