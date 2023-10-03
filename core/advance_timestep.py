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
DIFFERENT_RUN_IN_PLACE_FOR_BURN_IN = spc_config['DIFFERENT_RUN_IN_PLACE_FOR_BURN_IN']=='True'
ENS_END_DATE_datetime = datetime.strptime(ENS_END_DATE, "%Y%m%d")
DO_RERUN = spc_config["DO_VARON_RERUN"] == "True"
if DO_RERUN:
	number_of_windows_to_rerun = int(spc_config["number_of_windows_to_rerun"])
	do_approx = False
	APPROXIMATE_VARON_RERUN = spc_config["APPROXIMATE_VARON_RERUN"] == "True"
	if APPROXIMATE_VARON_RERUN:
		with open(f"{path_to_scratch}/APPOXIMATION_STAGE") as f:
			lines = f.readlines()
			if lines[0] == 'true':
				do_approx = True

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
	if DO_RERUN:
		delta = timedelta(hours=int(ASSIM_TIME))
		if APPROXIMATE_VARON_RERUN: 
			if do_approx: #If we approximated, then the restart is now extrapolated. Advance one period ahead.
				start_datetime = ig_end_datetime
				start_string = ig_endstring
				end_datetime = start_datetime+delta
				end_string = end_datetime.strftime("%Y%m%d %H%M%S")
			else: #If we are approximating but didn't this time, then we did the assimilation. Next run will rerun previous period to simulate conc changes.
				if periodstr=="POSTFIRST": #This is the only time we don't rerun the previous window. Advance forward one window to kick things off.
					start_datetime = ig_end_datetime
					start_string = ig_endstring
					end_datetime = start_datetime+delta
					end_string = end_datetime.strftime("%Y%m%d %H%M%S") #Will run LETKF next.
				else:
					start_datetime = ig_start_datetime
					start_string = ig_startstring
					end_datetime = ig_end_datetime
					end_string = ig_endstring
		else: #if we are doing the Varon rerun but not approximating, we start n assimilation periods before (default 1, set by user) and end one after current end slot.
			start_datetime = ig_end_datetime-(number_of_windows_to_rerun*delta)
			start_string = start_datetime.strftime("%Y%m%d %H%M%S")
			end_datetime = ig_end_datetime+delta
			end_string = end_datetime.strftime("%Y%m%d %H%M%S")
	else: 
		if do_rip_aw and (periodstr=="ASSIM"): #Do RIP if we are in a normal phase.
			start_datetime = ig_start_datetime+timedelta(hours=int(actual_aw)) #Advance start by assimilation window
			start_string = start_datetime.strftime("%Y%m%d %H%M%S")
		elif do_rip_aw and DIFFERENT_RUN_IN_PLACE_FOR_BURN_IN and (periodstr=="POSTBURN"): #If we are just after the burn in period, still advance using the burn in RIP.
			aw = int(spc_config['rip_burnin_update_time'])
			start_datetime = ig_start_datetime+timedelta(hours=int(aw)) #Advance start by burn in assimilation window
			start_string = start_datetime.strftime("%Y%m%d %H%M%S")
		else: #If we are just after the first assimilation period (POSTFIRST), don't worry about RIP. Same if there is no RIP to begin with.
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

#Toggle approximation stage if rerunning with approximation; don't do this for the first run period (POSTFIRST)
if DO_RERUN and APPROXIMATE_VARON_RERUN and periodstr!="POSTFIRST":
	if do_approx:
		with open(f"{path_to_scratch}/APPOXIMATION_STAGE", "w") as f:
			f.write("false") #Approximated this time; next time we won't
	else:
		with open(f"{path_to_scratch}/APPOXIMATION_STAGE", "w") as f:
			f.write("true") #Assimilated this time; next time we approximate

with open(f"{parent_dir}/scratch/INPUT_GEOS_TEMP", "w") as f:
	f.write(f"{start_string}\n")
	f.write(f"{end_string}\n")
	f.close()

