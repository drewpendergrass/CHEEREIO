import json

def getSpeciesConfig():
	with open('../ens_config.json') as f:
		data = json.load(f)
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

#Get index values within the localization range
def getIndsOfInterest(latind,lonind):
	data = getSpeciesConfig()
	lat,lon = getLatLonVals()
	latval = lat[latind]
	lonval = lon[lonind]

	pass