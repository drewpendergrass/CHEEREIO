import json
import geopy.distance as gd
import numpy as np

#Load in settings from ens_config.json. If settings_to_override is active (entry override is True)
#the overwrite entries provided
def getSpeciesConfig():
	with open('../ens_config.json') as f:
		data = json.load(f)
	with open('../settings_to_override.json') as f:
		over_data = json.load(f)
	if over_data['override'] == "True":
		for key in list(over_data.keys()):
			if key != 'override':
				data[key] = over_data[key]
	return data

#Get the latitude and longitude list for a particular core (indexed by ensemble and core)
def getLatLonList(ensnum,corenum):
	data = getSpeciesConfig()
	with open(f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch/latlon_par.json") as f:
		gridsplit = json.load(f)
	return [gridsplit[f'{ensnum}'][f'{corenum}']['lat'],gridsplit[f'{ensnum}'][f'{corenum}']['lon']]

#Get the lat lon values for the grid from the JSON
def getLatLonVals(data=None):
	if not data:
		data = getSpeciesConfig()
	with open(f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch/latlon_vals.json") as f:
		ll_data = json.load(f)
	return [ll_data['lat'],ll_data['lon']]

#Inputs are in degrees
def calcDist_km(lat1,lon1,lat2,lon2):
	coords_1 = (lat1,lon1)
	coords_2 = (lat2,lon2)
	d = gd.great_circle(coords_1, coords_2).km #This approximation is fine.
	return d

def makeLatLonGridWithMask(gridlabel,mask_coast_bool="True"):
	mask = None
	subhalf = False
	subquarter = False
	#Create traditional GEOS-Chem longitude and latitude centers, as specified by the settings in ens_config.json
	if gridlabel == '4.0x5.0':
		lon = np.arange(-180.0,176.0, 5.0)
		lat = np.concatenate([[-89.0],np.arange(-86.0,87.0, 4.0), [89.0]])
		if mask_coast_bool=="True":
			mask = np.transpose(np.genfromtxt('../templates/landmask_4x5_gcgrid_landfraction_gt_75pct.csv',delimiter=','))
		else:
			mask = np.transpose(np.genfromtxt('../templates/landmask_4x5_gcgrid.csv',delimiter=','))
		mask = np.where(mask==0)
	elif gridlabel == '2.0x2.5':
		lon = np.arange(-180.0,178.0, 2.5)
		lat = np.concatenate([[-89.5],np.arange(-88.0,89.0, 2.0), [89.5]])
		if mask_coast_bool=="True":
			mask = np.transpose(np.genfromtxt('../templates/landmask_2x2p5_gcgrid_landfraction_gt_75pct.csv',delimiter=','))
		else:
			mask = np.transpose(np.genfromtxt('../templates/landmask_2x2p5_gcgrid.csv',delimiter=','))
		mask = np.where(mask==0)
	elif gridlabel == '1x1':
		lon = np.arange(-179.5,180.0, 1.0)
		lat = np.arange(-89.5,90, 1.0)
		mask = np.transpose(np.genfromtxt('../templates/landmask_1x1_gcgrid.csv',delimiter=','))
		mask = np.where(mask==0)
	elif (gridlabel == '0.5x0.625') | (gridlabel == 'MERRA2'): #MERRA2 NATIVE GRID
		lon = -180.0 + (0.625*np.arange(0.0,576.0,1.0))
		lat = -90.0 + (0.5*np.arange(0.0,361.0,1.0))
		mask = np.transpose(np.genfromtxt('../templates/landmask_0p5x0p625_gcgrid.csv',delimiter=','))
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
		mask = np.transpose(np.genfromtxt('../templates/landmask_0p25x0p3125_gcgrid.csv',delimiter=','))
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
	#Handle the high res nested grid case.
	if subhalf or subquarter:
		if subhalf:
			fulllon = -180.0 + (0.625*np.arange(0.0,576.0,1.0))
			fulllat = -90.0 + (0.5*np.arange(0.0,361.0,1.0))
			mask = np.transpose(np.genfromtxt('../templates/landmask_0p5x0p625_gcgrid.csv',delimiter=','))
		elif subquarter:
			fulllon = -180.0 + (0.3125*np.arange(0.0,1152.0,1.0))
			fulllat = -90.0 + (0.25*np.arange(0.0,721.0,1.0))
			mask = np.transpose(np.genfromtxt('../templates/landmask_0p25x0p3125_gcgrid.csv',delimiter=','))
		lonok = np.where((fulllon>=np.min(lon)) & (fulllon<=np.max(lon)))[0]
		latok = np.where((fulllat>=np.min(lat)) & (fulllat<=np.max(lat)))[0]
		mask = mask[lonok,latok]
		mask = np.where(mask==0)
	return [lon,lat,mask]

def makeDistMat(instruction = 'file', verbose = 1, makeDistMat = True):
	#In this case, use the saved lat/lon data in the scratch folder
	if instruction == 'file':
		lat,lon = getLatLonVals(data)
	#Otherwise, generate the lat lon data according to the supplied gridlabel
	else:
		lon,lat,_ = makeLatLonGridWithMask(instruction)
	X,Y = np.meshgrid(lon,lat)
	XY = np.column_stack((np.ndarray.flatten(X),np.ndarray.flatten(Y)))
	if makeDistMat:
		numpoints = np.shape(XY)[0]
		distmat = np.zeros((numpoints,numpoints))
		for i in range(numpoints):
			if verbose >= 1:
				print(f'Calculating row {i} of {numpoints}')
			for j in range(i, numpoints):
				distval = calcDist_km(XY[i,1],XY[i,0],XY[j,1],XY[j,0])
				distmat[i,j] = distval
				distmat[j,i] = distval
		return [distmat,XY]
	#Option to just get lon/lat coordinates from distmat file
	else:
		return XY

#Get index values within the localization range
#If negate is true, then get index values outside the localization range
def getIndsOfInterest(latind,lonind,negate=False):
	data = getSpeciesConfig()
	lat,lon = getLatLonVals(data)
	latval = lat[latind]
	lonval = lon[lonind]
	loc_rad = float(data['LOCALIZATION_RADIUS_km'])
	distgrid = np.zeros((len(lat),len(lon)))
	for i in range(len(lon)):
		dists_col = np.array([calcDist_km(latval,lonval,a,lon[i]) for a in lat])
		distgrid[:,i] = dists_col
	if negate:
		valid_inds = np.where(distgrid>loc_rad)
	else:
		valid_inds = np.where(distgrid<=loc_rad)
	return valid_inds