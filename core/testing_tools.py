from Assimilator import Assimilator
import numpy as np
import pandas as pd
import json
import settings_interface as si 
import xarray as xr

#Makes an assimilator object
def makeAssimilator():
	a = Assimilator('20190108_0000',2,1) #we aren't using ens or core num here.
	return a

#Overrides settings without modifying ens_config. If overwrite is true, it deletes previous adjustments
def overrideSettings(settings_to_override, overwrite = False):
	with open('../settings_to_override.json') as f:
		over_data = json.load(f)
	if overwrite:
		over_data = settings_to_override #Flat overwrite
	else: #If not overwriting, loop through and add keys individually
		for key in list(settings_to_override.keys()):
			if key != 'override':
				over_data[key] = settings_to_override[key]
	over_data['override'] = "True" #Activate overriding capability
	with open('../settings_to_override.json', 'w') as f:
		json.dump(over_data, f, ensure_ascii=False, indent=4)

def setupPytestSettings(simtype='methane'):
	assimdir = si.getSpeciesConfig()['ASSIM_PATH'] #Assim path, so that we can overwrite text in override json file.
	if simtype=='methane':
		with open('../tests/data_for_tests/METHANE_TEST/methane_settings_to_override.json') as f:
			settings_to_override = json.load(f)
	else:
		raise ValueError('Cheereio Pytest setting not supported')
	for key in settings_to_override:
		if '{ASSIM_PATH}' in settings_to_override[key]: #Overwrite {ASSIM_PATH} in settings to override file with real path.
			settings_to_override[key] = settings_to_override[key].replace('{ASSIM_PATH}',assimdir)
		if key == "TROPOMI_dirs": #Recurse down another level
			for key2 in settings_to_override[key]: #Overwrite {ASSIM_PATH} in settings to override file with real path.
				settings_to_override[key][key2] = settings_to_override[key][key2].replace('{ASSIM_PATH}',assimdir)
	overrideSettings(settings_to_override,overwrite=True)

def turnOffOverride():
	with open('../settings_to_override.json') as f:
		over_data = json.load(f)
	over_data['override'] = "False" #Turn off overriding capability
	with open('../settings_to_override.json', 'w') as f:
		json.dump(over_data, f, ensure_ascii=False, indent=4)


def makeMiniFakeDataSet(nlat = 3, nlon = 3, nlev = 3, ntime = 3):
	lat = np.linspace(0,10, nlat)
	lon = np.linspace(0,10, nlon)
	lev = np.linspace(0,10, nlev)
	time = pd.date_range(start='2022-08-01',end='2022-08-08',periods=ntime).values
	values = np.arange(ntime*nlev*nlat*nlon).reshape((ntime,nlev,nlat,nlon)) #0 through nlatnlon on timestamp 0, and so on.
	ds = xr.Dataset(
		{"SpeciesConc_TEST": (("time","lev","lat","lon"), values,{"long_name": "Test data", "units":"1"}),
		"Met_PEDGE": (("time","lev","lat","lon"), values,{"long_name": "Test data", "units":"1"})},
		coords={
			"time": (["time"], time, {"long_name": "time", "calendar": "standard"}),
			"lev": (["lev"], lev,{"long_name": "Levels"}),
			"lat": (["lat"], lat,{"long_name": "Latitude", "units":"degrees_north"}),
			"lon": (["lon"], lon,{"long_name": "Longitude", "units":"degrees_east"})
		},
		attrs={
			"Title":"Test dataset",
		}
	)
	return ds

def makeMiniFakeObsData(latlocs,lonlocs,ntime):
	obsdata = {'latitude':np.array([]),'longitude':np.array([]),'utctime':np.array([])}
	for latloc,lonloc in zip(latlocs,lonlocs):
		obsdata['latitude'].append(np.repeat(latloc,ntime))
		obsdata['longitude'].append(np.repeat(lonloc,ntime))
		obsdata['utctime'].append(pd.date_range(start='2022-08-01',end='2022-08-08',periods=ntime).values)
	return obsdata

#Walks through with extensive print statements an assimilation cycle
def walkThroughAssimilation(assim,latind=65,lonind=24): #default is a point in northern California, arbitrary 
	print('*************************************************')
	print('*******ASSIMILATION WALKTHROUGH, VERBOSE*********')
	print('*************************************************')
	assim.verbose=2
	print(f'Coordinates for this walkthrough: {(assim.gt[1].getLat()[latind],assim.gt[1].getLon()[lonind])}')
	assim.prepareMeansAndPerts(latind,lonind)
	print(f'ybar has value {assim.ybar_background}')
	print(f'ydiff has value {assim.ydiff}')
	assim.makeR(latind,lonind)
	sqrtdiag = np.sqrt(np.diag(assim.R))
	print(f'The square root of the diagonal of R is {sqrtdiag}')
	print(f'Relative to ybar, that diagonal is value {sqrtdiag/assim.ybar_background}')
	assim.makeC()
	assim.makePtildeAnalysis()
	assim.makeWAnalysis()
	assim.makeWbarAnalysis()
	assim.adjWAnalysis()
	assim.makeAnalysisCombinedEnsemble()
	analysisSubset,backgroundSubset,analysisPertSubset,backgroundPertSubset = assim.getAnalysisAndBackgroundColumn(latind,lonind,doBackground=True,doPerts=True)
	meananalysis = np.mean(analysisSubset,axis=1)
	meanbackground = np.mean(backgroundSubset,axis=1)
	print(f'Background column has dimension {np.shape(backgroundSubset)} and ensemble mean value {meanbackground}')
	print(f'Analysis column has dimension {np.shape(analysisSubset)} and ensemble mean value {meananalysis}')
	print(f'Ensemble mean analysis minus ens mean background has value {meananalysis-meanbackground}')
	print(f'This represents a percent difference of {100*((meananalysis-meanbackground)/meanbackground)}%')
	dofs = assim.calculateDOFS(analysisPertSubset,backgroundPertSubset)
	print(f'DOFS has value {dofs}')
	analysisSubsetAdjusted = assim.applyAnalysisCorrections(analysisSubset,backgroundSubset,latind,lonind)
	meananalysisSubsetAdjusted = np.mean(analysisSubsetAdjusted,axis=1)
	print(f'Old analysis ensemble mean for column was {meananalysis}.')
	print(f'After postprocessing, new analysis ensemble mean for column is {meananalysisSubsetAdjusted}.')
	print(f'This represents a percent difference of {100*((meananalysisSubsetAdjusted-meananalysis)/meananalysis)}%')


