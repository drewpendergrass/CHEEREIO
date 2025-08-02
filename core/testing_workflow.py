#This code automates the walkThroughAssimilation() fxn in testing_tools for ease-of-use (e.g. batch submission).

import pickle
import sys
import argparse
from Assimilator import Assimilator 
import testing_tools as tt

def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser(description='Lightweight wrapper of postprocess script for bigy plots (OmF, obs totals).')
parser.add_argument('-o', '--file_out', type=str, default="../../scratch/assimilator.pkl", help='Name (and optional path) of pkl file to save of Assimilator. By default, saves in scratch in file "assimilator.pkl".')
parser.add_argument('-d', '--date', type=str, help='Date for feeding into Assimilator constructor (e.g. 20160108_0000).')
parser.add_argument('-rd', '--rip_date', type=str, default=None, help='RIP date for feeding into Assimilator constructor (e.g. 20160108_0000). If not using RIP, neglect')
parser.add_argument('-ens', '--ens_num', type=int,default=1, help='Ensemble number for Assimilator.')
parser.add_argument('-core', '--core_num', type=int,default=1, help='Core number for Assimilator.')
parser.add_argument('-lat', '--lat_ind', type=int, help='Lat ind for assimilation walkthrough.')
parser.add_argument('-lon', '--lon_ind', type=int, help='Lon ind for assimilation walkthrough.')

args = parser.parse_args()

assim = Assimilator(args.date,args.ens_num,args.core_num, args.rip_date)

#Pickle assimilator
with open(args.file_out, 'wb') as f:
    pickle.dump(assim, f)

print('Saved Assimilator object! You can load this and manipulate for further testing. Use the following code:')
print('')
print('')
print(f"import pickle")
print('')
print(f"with open('{args.file_out}', 'rb') as f:")
print(f"     assim = pickle.load(f)")
print('')
print('')
tt.walkThroughAssimilation(assim,latind=args.lat_ind,lonind=args.lon_ind)

