import argparse
import postprocess_tools as pt
import numpy as np
import pickle
import sys
sys.path.append('../core')
import settings_interface as si 

def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description='CHEEREIO runtime snapshot')
parser.add_argument("--gossipgirl", type=str2bool, nargs='?',const=True, default=False,help="Activate Gossip Girl mode.")
args = parser.parse_args()
gg_mode = args.gossipgirl

if gg_mode:
	print('')
	print('')
	print('Hello Upper East Siders. Gossip Girl here. And I have the biggest news ever.')
	print('One of my many sources, SLURM, sends us this: “Spotted on the computational cluster, bags in hand: CHEEREIO.”')
	print('Was it only a year ago our It Girl mysteriously disappeared for "data assimilation”? And just as suddenly, she’s back.')
	print('Don’t believe me? See for yourselves. We will now be showing a snapshot of CHEEREIO at runtime, saving a few figures and movies to postprocessing.')
	print('')
	print('')
else:
	print('')
	print('')
	print('Begin CHEEREIO runtime snapshot. A few diagnostic figures and movies will be saved to postprocessing.')
	print('')
	print('')

data = si.getSpeciesConfig()
pp_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/postprocess"
ens_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs"
ASSIM_TIME = data['ASSIM_TIME']

print('Loading simulated observation and observation dictionaries...')
bigy = pt.makeYEachAssimPeriod(path_to_bigy_subsets=f"{pp_dir}/bigy",assim_time=int(ASSIM_TIME))
print('Simulated observation and observation dictionaries loaded.')

nEnsemble = int(data['nEnsemble'])
observed_species = data['OBSERVED_SPECIES']
obs_units = data['OBSERVATION_UNITS']
useControl=data['DO_CONTROL_RUN']=="true"
controlInEns = data['DO_CONTROL_WITHIN_ENSEMBLE_RUNS']=="true"

if useControl and controlInEns:
	control_dir = f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs/{data['RUN_NAME']}_0000"
else:
	control_dir = None

print('Plotting timeseries of observations against assimilation')

for spec in observed_species:
	print('')
	print('')
	if control_dir is not None:
		outfile = f'{pp_dir}/SNAPSHOT_observations_ts_compare_{spec}_w_control.png'
		pt.tsPlotSatCompare(bigy,spec,nEnsemble,unit=obs_units[spec],observer_name=data['OBS_TYPE'][spec],useControl=True,outfile=f'{pp_dir}/SNAPSHOT_observations_ts_compare_{spec}_w_control.png')
	else:
		outfile = f'{pp_dir}/SNAPSHOT_observations_ts_compare_{spec}.png'
		pt.tsPlotSatCompare(bigy,spec,nEnsemble,unit=obs_units[spec],observer_name=data['OBS_TYPE'][spec],useControl=False,outfile=outfile)
	print(f'Saved a snapshot at {outfile}')

print('')
print('')

print('Timeseries plots complete.')

print('Aggregating scale factors from across the ensemble...')
pt.combineScaleFactors(ens_dir,pp_dir,flag_snapshot=True)
print('Scale factor aggregation complete.')



