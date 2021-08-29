import json

def getSpeciesConfig():
	with open('../ens_config.json') as f:
		data = json.load(f)
	return data

def getLatLonList(ensnum):
	data = getSpeciesConfig()
	with open(f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch/latlon_par.json") as f:
		gridsplit = json.load(f)
	return [gridsplit[f'{ensnum}']['lat'],gridsplit[f'{ensnum}']['lon']]
