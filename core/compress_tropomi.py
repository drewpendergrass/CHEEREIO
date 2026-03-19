import xarray as xr
import numpy as np
from glob import glob
from tropomi_tools import *
import argparse
import os

parser = argparse.ArgumentParser(description='Combine SpeciesConc data into an ensemble mean.')
parser.add_argument('-i', '--path_to_tropomi_input', type=str, help='Path to directory containing TROPOMI input files.')
parser.add_argument('-o', '--path_to_tropomi_output', type=str, help='Path to directory to save compressed TROPOMI files, following identical format and .')
parser.add_argument('-s', '--species', type=str,help='Species to compress.')
args = parser.parse_args()

obs_list = glob(f'{args.path_to_tropomi_input}/**/S5P_*.nc', recursive=True)
obs_list.sort()

def getShapeName(data,dimlen):
	if dimlen==data["qa_value"].shape[0]:
		shape = "obs"
	elif dimlen==data["column_AK"].shape[1]:
		shape = "level"
	elif dimlen==data["pressures"].shape[1]:
		shape = "level_edge"
	else:
		raise ValueError("Did not recognize shape")
	return shape


for obs in obs_list:
	print(f'Loading {obs}')
	data = read_tropomi(obs, args.species, filterinfo=None, includeObsError = True)
	to_save = {}
	for key in data:
		if len(data[key].shape)==1:
			shape = getShapeName(data,data[key].shape[0])
		elif len(data[key].shape)==2:
			shape = (getShapeName(data,data[key].shape[0]),getShapeName(data,data[key].shape[1]))
		else:
			raise ValueError(f"Did not recognize shape for shape {key}")
		to_save[key] = (shape,data[key])
	ds = xr.Dataset(
		to_save,
		coords={
			"obs": np.arange(len(data["qa_value"])),
			"level": np.arange(data["column_AK"].shape[1]),
			"level_edge": np.arange(data["pressures"].shape[1])
		}
	)
	rel_path = os.path.relpath(obs, args.path_to_tropomi_input)
	out_path = os.path.join(args.path_to_tropomi_output, rel_path)
	os.makedirs(os.path.dirname(out_path), exist_ok=True)
	print(f'Saving {out_path}')
	ds.to_netcdf(out_path)
