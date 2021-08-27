import json

def getSpeciesConfig():
	with open('../ens_config.json') as f:
		data = json.load(f)
	return data
