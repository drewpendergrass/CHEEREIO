from glob import glob
import settings_interface as si 
import sys 
from GT_Container import GT_Container
import time
import numpy as np
import os
from datetime import datetime,timedelta

timestamp = str(sys.argv[1]) #Time to assimilate. Expected in form YYYYMMDD_HHMM, UTC time.

data = si.getSpeciesConfig()
path_to_scratch = f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch"


DO_ADDL_INFLATION = False
if data['Activate_Relaxation_To_Prior_Spread'].lower()=='true':
	for species in data['species_not_in_statevec_to_RTPS']:
		#Don't inflate species in state_vector_conc.
		if species not in data['STATE_VECTOR_CONC']:
			DO_ADDL_INFLATION=True
			break

DO_RERUN = data["DO_VARON_RERUN"] == "True"
if DO_RERUN:
	number_of_windows_to_rerun = int(data["number_of_windows_to_rerun"])
	APPROXIMATE_VARON_RERUN = data["APPROXIMATE_VARON_RERUN"] == "True"

with open(f"{path_to_scratch}/ACTUAL_RUN_IN_PLACE_ASSIMILATION_WINDOW") as f:
    lines = f.readlines()

actual_aw = float(lines[0])
do_rip_aw = False

if not np.isnan(actual_aw):
	actual_aw = int(actual_aw)
	do_rip_aw = True

#Calculate time to use for restart load if we are doing run in place or rerun; different from timestamp
if DO_RERUN:
	if APPROXIMATE_VARON_RERUN: #never pass along a different restart time if approximating, regardless of approximation stage.
		timestamp_restart = None
		if DO_ADDL_INFLATION:
			#Addl inflation requires timestamp of background concentrations. For standard simulations, this is one window before the restart.
			timestamp_datetime = datetime.strptime(timestamp, "%Y%m%d_%H%M")
			delta = timedelta(hours=int(data['ASSIM_TIME']))
			timestamp_background_dt = timestamp_datetime-delta #Background restart is 1 timestep behind current timestamp.
			timestamp_background = timestamp_background_dt.strftime("%Y%m%d_%H%M") 
	else:
		timestamp_datetime = datetime.strptime(timestamp, "%Y%m%d_%H%M")
		delta = timedelta(hours=int(data['ASSIM_TIME']))
		timestamp_restart_dt = timestamp_datetime-(number_of_windows_to_rerun*delta) #Restart we are assimilating is n timesteps behind current timestamp (default 1, because rerunning from n assimilation periods behind).
		timestamp_restart = timestamp_restart_dt.strftime("%Y%m%d_%H%M") 
		if DO_ADDL_INFLATION:
			#Addl inflation requires timestamp of background concentrations. For rerun, this is one window before the restart.
			timestamp_background_dt = timestamp_datetime-((number_of_windows_to_rerun+1)*delta) #Background restart is n+1 timesteps behind current timestamp.
			timestamp_background = timestamp_background_dt.strftime("%Y%m%d_%H%M") 
elif do_rip_aw:
	ASSIM_TIME = int(data['ASSIM_TIME'])
	backwards = ASSIM_TIME-actual_aw #How many hours backwards should we look from current timestamp for restart to use in building state vector?
	timestamp_datetime = datetime.strptime(timestamp, "%Y%m%d_%H%M")
	delta = timedelta(hours=int(backwards))
	timestamp_restart_dt = timestamp_datetime-delta
	timestamp_restart = timestamp_restart_dt.strftime("%Y%m%d_%H%M")
	if DO_ADDL_INFLATION:
		#Addl inflation requires timestamp of background concentrations. For run in place, this is the assimilation window before the restart.
		timestamp_background_dt = timestamp_datetime-timedelta(hours=ASSIM_TIME)
		timestamp_background = timestamp_background_dt.strftime("%Y%m%d_%H%M")
else:
	timestamp_restart = None
	if DO_ADDL_INFLATION:
		#Addl inflation requires timestamp of background concentrations. For standard simulations, this is one window before the restart.
		timestamp_datetime = datetime.strptime(timestamp, "%Y%m%d_%H%M")
		delta = timedelta(hours=int(data['ASSIM_TIME']))
		timestamp_background_dt = timestamp_datetime-delta #Background restart is 1 timestep behind current timestamp.
		timestamp_background = timestamp_background_dt.strftime("%Y%m%d_%H%M") 

#Change timestamp for GC_Translator
if timestamp_restart is not None:
	timestamp = timestamp_restart 

SaveDOFS = data["SaveDOFS"] == "True"

dateval = timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]

print(f'One core is gathering columns to overwrite at time {dateval}.')
start = time.time()
wrapper = GT_Container(timestamp)
end = time.time()
print(f'Core gathered columns and ensemble in {end - start} seconds. Begin saving.')
start = time.time()
wrapper.reconstructAnalysisEnsemble()
wrapper.updateRestartsAndScalingFactors()
if DO_ADDL_INFLATION:
	wrapper.performAdditionalInflation(timestamp_background)
wrapper.saveRestartsAndScalingFactors()

if SaveDOFS:
	npy_dofs_files = glob(f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs/logs/dofs_scratch/*.npy")
	npy_dofs = np.stack([np.load(file) for file in npy_dofs_files])
	dofs_combined = np.nansum(npy_dofs,axis=0) #Perfectly nonoverlapping, so just nansum to combine
	np.save(f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs/logs/dofs_complete/combined_dofsgrid_{timestamp}.npy",dofs_combined)
	#Delete all scratch files now that they are combined
	for file in npy_dofs_files:
		os.remove(file)

end = time.time()
print(f'Saved updated restarts and emissions in {end - start} seconds. We can cleanup and resume GC now!')
