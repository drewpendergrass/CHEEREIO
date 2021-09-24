import json
import geopy.distance as gd
import numpy as np

def getSpeciesConfig(testing=False):
	if testing:
		with open('../testing/test_config.json') as f:
			data = json.load(f)
	else:
		with open('../ens_config.json') as f:
			data = json.load(f)
	return data

#Get the latitude and longitude list for a particular core (indexed by ensemble and core)
def getLatLonList(ensnum,corenum,testing=False):
	data = getSpeciesConfig(testing)
	with open(f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch/latlon_par.json") as f:
		gridsplit = json.load(f)
	return [gridsplit[f'{ensnum}'][f'{corenum}']['lat'],gridsplit[f'{ensnum}'][f'{corenum}']['lon']]

#Get the lat lon values for the grid from the JSON
def getLatLonVals(data=None,testing=False):
	if not data:
		data = getSpeciesConfig(testing)
	with open(f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch/latlon_vals.json") as f:
		ll_data = json.load(f)
	return [ll_data['lat'],ll_data['lon']]

#Inputs are in degrees
def calcDist_km_fromscratch(lat1,lon1,lat2,lon2):
	coords_1 = (lat1,lon1)
	coords_2 = (lat2,lon2)
	d = gd.distance(coords_1, coords_2).km
	return d

#Inputs are in degrees
def calcDist_km(lat1,lon1,lat2,lon2,testing=False):
	data = getSpeciesConfig(testing)
	lat,lon = getLatLonVals(data=data)
	dist_tensor = np.load(f"{data['MY_PATH']}/{data['RUN_NAME']}/dist_tensor.npy")
	latinds = [np.where(lat==lat1)[0],np.where(lat==lat2)[0]]
	loninds = [np.where(lon==lon1)[0],np.where(lon==lon2)[0]]
	return dist_tensor[min(latinds),min(loninds),max(latinds),max(loninds)]


def makeDistTensor(testing=False):
	data = getSpeciesConfig(testing)
	lat,lon = getLatLonVals(data=data)
	dist_tensor = np.zeros((len(lat),len(lon),len(lat),len(lon)))
	for i in range(len(lat)):
		lat1 = lat[i]
		for j in range(len(lon)):
			lon1 = lon[j]
			for k in range(i,len(lat)):
				lat2 = lat[k]
				for l in range(j,len(lon)):
					lon2 = lon[l]
					dist_tensor[i,j,k,l] = calcDist_km_fromscratch(lat1,lon1,lat2,lon2)
	np.save(f"{data['MY_PATH']}/{data['RUN_NAME']}/dist_tensor.npy",dist_tensor)

#Get index values within the localization range
#If negate is true, then get index values outside the localization range
def getIndsOfInterest(latind,lonind,negate=False,testing=False):
	data = getSpeciesConfig(testing)
	lat,lon = getLatLonVals(data)
	latval = lat[latind]
	lonval = lon[lonind]
	loc_rad = float(data['LOCALIZATION_RADIUS_km'])
	distgrid = np.zeros((len(lat),len(lon)))
	for i in range(len(lon)):
		dists_col = np.array([calcDist_km(latval,lonval,a,lon[i],testing) for a in lat])
		distgrid[:,i] = dists_col
	if negate:
		valid_inds = np.where(distgrid>loc_rad)
	else:
		valid_inds = np.where(distgrid<=loc_rad)
	return valid_inds