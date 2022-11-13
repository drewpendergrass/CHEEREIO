import json

#Load in settings from ens_config.json. If settings_to_override is active (entry override is True)
#the overwrite entries provided
def getSpeciesConfig():
	with open('../ens_config.json') as f:
		data = json.load(f)
	#Load extensions and add if turned on
	extensions = data['Extensions']
	for key in list(extensions.keys()):
		if extensions[key] == 'True':
			data = addExtension(data,f'../extensions/{key}_extension.json')
	#Override settings --- mostly used by the automated testing utility.
	with open('../settings_to_override.json') as f:
		over_data = json.load(f)
	if over_data['override'] == "True":
		for key in list(over_data.keys()):
			if key != 'override':
				data[key] = over_data[key]
	return data

#Add extension setting data to main data list
def addExtension(data,file):
	with open(file) as f:
		extdata = json.load(f)
	for key in list(extdata.keys()):
		data[key] = extdata[key]
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
