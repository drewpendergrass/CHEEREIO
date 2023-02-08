import settings_interface as si 
import os.path
from datetime import datetime,timedelta
import sys
import numpy as np

periodstr = str(sys.argv[1])
spc_config = si.getSpeciesConfig()

parent_dir = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}"
path_to_scratch = f"{parent_dir}/scratch"
ens_dir = f"{parent_dir}/ensemble_runs"
ASSIM_TIME = spc_config['ASSIM_TIME']
START_DATE = spc_config['START_DATE']
ENS_END_DATE = spc_config['END_DATE']
ASSIM_START_DATE = spc_config['ASSIM_START_DATE']
ENS_SPINUP_START = spc_config['ENS_SPINUP_START']
ENS_SPINUP_END = spc_config['ENS_SPINUP_END']
ENS_END_DATE_datetime = datetime.strptime(ENS_END_DATE, "%Y%m%d")

with open(f"{path_to_scratch}/ACTUAL_RUN_IN_PLACE_ASSIMILATION_WINDOW") as f:
    lines = f.readlines()

actual_aw = float(lines[0])
do_rip_aw = False

if not np.isnan(actual_aw):
	actual_aw = int(actual_aw)
	do_rip_aw = True

#This tool also handles whether we scale at the end of the burn in period, producing signal files to trigger the appropriate processes
DO_BURN_IN = spc_config['DO_BURN_IN'] == "true"
if DO_BURN_IN:
	BURN_IN_END = spc_config['BURN_IN_END']
	BURN_IN_END_datetime = datetime.strptime(BURN_IN_END, "%Y%m%d")

if periodstr=="FIRST":
	start_string = f"{START_DATE} 000000"
	start_datetime = datetime.strptime(START_DATE, "%Y%m%d")
	end_string = f"{ASSIM_START_DATE} 000000"
	end_datetime = datetime.strptime(ASSIM_START_DATE, "%Y%m%d")
elif periodstr=="SPINUP":
	start_string = f"{ENS_SPINUP_START} 000000"
	start_datetime = datetime.strptime(ENS_SPINUP_START, "%Y%m%d")
	end_string = f"{ENS_SPINUP_END} 000000"
	end_datetime = datetime.strptime(ENS_SPINUP_END, "%Y%m%d")
else:
	#Load INPUT GEOS TEMP
	with open(f"{parent_dir}/scratch/INPUT_GEOS_TEMP") as f:
		lines = f.readlines()
	ig_startstring = lines[0].rstrip()
	ig_start_datetime = datetime.strptime(ig_startstring, "%Y%m%d %H%M%S")
	ig_endstring = lines[1].rstrip()
	ig_end_datetime = datetime.strptime(ig_endstring, "%Y%m%d %H%M%S")
	if do_rip_aw and (periodstr=="ASSIM"): #If we are just after the first assimilation period (POSTFIRST), don't worry about RIP.
		start_datetime = ig_start_datetime+timedelta(hours=int(actual_aw)) #Advance start by assimilation window
		start_string = start_datetime.strftime("%Y%m%d %H%M%S")
	else:
		start_datetime = ig_end_datetime
		start_string = ig_endstring
	delta = timedelta(hours=int(ASSIM_TIME))
	end_datetime = start_datetime+delta
	end_string = end_datetime.strftime("%Y%m%d %H%M%S")

#Check if the upcoming assimilation period will end  after the burn in period ends
if DO_BURN_IN:
	#If burn in period will be over at end of next GC cycle
	if end_datetime >= BURN_IN_END_datetime:
		#If we have already processed the burn in period in GC
		if os.path.exists(f"{parent_dir}/scratch/BURN_IN_PERIOD_PROCESSED"):
			#We only want to scale concentrations once after the burn in. If the BURN_IN_SCALING_COMPLETE file is missing but the BURN_IN_PERIOD_PROCESSED is present
			#then CHEEREIO has already done simple scaling when this script is called. Create a signal file to turn this thing off
			if not os.path.exists(f"{parent_dir}/scratch/BURN_IN_SCALING_COMPLETE"):
				with open(f"{parent_dir}/scratch/BURN_IN_SCALING_COMPLETE", "w") as f:
					f.write(f"GC has run for the entirety of the burn-in period, and CHEEREIO has done simple scaling at end of burn-in already.\n")
					f.close()
		else: #if we haven't already marked the burn in period processed, create a signal file. With just this file present, CHEEREIO will do the simple scaling
			with open(f"{parent_dir}/scratch/BURN_IN_PERIOD_PROCESSED", "w") as f:
				f.write(f"GC has run for the entirety of the burn-in period. If this file is present, and BURN_IN_SCALING_COMPLETE is missing, CHEEREIO will do simple scaling.\n")
				f.close()

#Check if we have finished the total ensemble run
if start_datetime >= ENS_END_DATE_datetime:
	with open(f"{parent_dir}/scratch/ENSEMBLE_COMPLETE", "w") as j:
		f.write("Ensemble completed; delete this file if you want to re-run.\n") #If so, save flag file to ensemble folder

with open(f"{parent_dir}/scratch/INPUT_GEOS_TEMP", "w") as f:
	f.write(f"{start_string}\n")
	f.write(f"{end_string}\n")
	f.close()

