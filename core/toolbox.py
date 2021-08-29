import json
import numpy as np
import xarray as xr
from glob import glob

def getSpeciesConfig():
	with open('../ens_config.json') as f:
		data = json.load(f)
	return data

def initialSplitEnsColumns():
	data = getSpeciesConfig()
	path_to_sim = f'{data['MY_PATH']}/{data['RUN_NAME']}/'
	subdirs = glob(f"{path_to_sim}ensemble_runs/*/")
	subdirs.remove(f"{path_to_sim}ensemble_runs/logs/")
	dirnames = [d.split('/')[-2] for d in subdirs]
	subdir_numbers = [int(n.split('_')[-1]) for n in dirnames]
	rst_dataset = xr.open_dataset(data['RESTART_FILE']) 
	lat = np.array(rst_dataset['lat'])
	lat_inds = np.arange(0,len(lat))
	lon = np.array(rst_dataset['lon'])
	lon_inds = np.arange(0,len(lon))
	lat_full_list = np.repeat(lat_inds,len(lon))
	lon_full_list = np.tile(lon_inds,len(lat))
	total_cells = len(lat)*len(lon)
	nens = data['nEnsemble']
	min_cells = np.floor(total_cells/nens)
	remainder_cells = total_cells%nens
	column_count_array = np.repeat(min_cells,nens)
	column_count_array[0:remainder_cells] = min_cells+1
	array_endpoints = np.insert(np.cumsum(column_count_array),0,0)
	split_lat_inds = []
	split_lon_inds = []
	for i in range(len(array_endpoints)-1):
		split_lat_inds.append(lat_full_list[i:(i+1)])
		split_lon_inds.append(lon_full_list[i:(i+1)])
	dict_to_save = {}
	for i in range(len(dirnames)):
		num = subdir_numbers[i]
		latlist = split_lat_inds[i]
		lonlist = split_lon_inds[i]
		dict_to_save[num] = {'lat':latlist,'lon':lonlist}
	# the json file where the output must be stored
	out_file = open(f"{path_to_sim}scratch/latlon_par.json", "w")
	json.dump(dict_to_save, out_file, indent = 6)
	out_file.close()

def getLatLonList(ensnum):
	data = getSpeciesConfig()
	with open(f'{data['MY_PATH']}/{data['RUN_NAME']}/scratch/latlon_par.json') as f:
		gridsplit = json.load(f)
	return [gridsplit[ensnum]['lat'],gridsplit[ensnum]['lon']]



