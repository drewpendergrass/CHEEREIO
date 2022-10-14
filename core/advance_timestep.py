import settings_interface as si 
import os.path
from datetime import datetime,timedelta
import sys

periodstr = str(sys.argv[1])
spc_config = si.getSpeciesConfig()

parent_dir = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}"
ens_dir = f"{parent_dir}/ensemble_runs"
ASSIM_TIME = spc_config['ASSIM_TIME']
ENS_END_DATE = spc_config['END_DATE']
ASSIM_START_DATE = spc_config['ASSIM_START_DATE']
ENS_SPINUP_END = spc_config['ENS_SPINUP_END']
ENS_END_DATE_datetime = datetime.strptime(ENS_END_DATE, "%Y%m%d")

#This tool also handles whether we scale at the end of the burn in period, producing signal files to trigger the appropriate processes
SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD = spc_config['SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD'] == "true"
if SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD:
	BURN_IN_END = spc_config['BURN_IN_END']
	BURN_IN_END_datetime = datetime.strptime(BURN_IN_END, "%Y%m%d")

with open(f"{parent_dir}/scratch/CURRENT_DATE_TIME") as f:
    start_string = f.readlines()[0].rstrip()

start_datetime = datetime.strptime(start_string, "%Y%m%d %H%M%S")
if periodstr=="FIRST":
	end_string = f"{ASSIM_START_DATE} 000000"
	end_datetime = datetime.strptime(ASSIM_START_DATE, "%Y%m%d")
elif periodstr=="SPINUP":
	end_string = f"{ENS_SPINUP_END} 000000"
	end_datetime = datetime.strptime(ENS_SPINUP_END, "%Y%m%d")
else:
	delta = timedelta(hours=int(ASSIM_TIME))
	end_datetime = start_datetime+delta
	end_string = end_datetime.strftime("%Y%m%d %H%M%S")

#Check if the upcoming assimilation period will end  after the burn in period ends
if SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD:
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

