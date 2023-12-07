#This utility removes links to scalefactors from a given HEMCO Config file in a run directory (path supplied by user).
#This is accomplished by replacing scalefactor nc files with ones (no scaling). 
#Used to turn a directory into a control/spinup run.

import sys

path_to_HEMCO = str(sys.argv[1])

with open(f'{path_to_HEMCO}/HEMCO_Config.rc') as f:
	lines = f.readlines()

for i in range(len(lines)):
	if lines[i].startswith('7') and ('ASSIM' in lines[i]):
		code, label = lines[i].split()[0:2]
		lines[i] = f'{code} {label} 1.0 - - - - 1 1' #this just replaces with a constant scalefactor of 1

with open(f'{path_to_HEMCO}/HEMCO_Config.rc', 'w') as f:
	for line in lines:
		f.write(line)
