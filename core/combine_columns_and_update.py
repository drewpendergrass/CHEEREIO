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

with open(f"{path_to_scratch}/ACTUAL_RUN_IN_PLACE_ASSIMILATION_WINDOW") as f:
    lines = f.readlines()

actual_aw = float(lines[0])
do_rip_aw = False

if not np.isnan(actual_aw):
	actual_aw = int(actual_aw)
	do_rip_aw = True

#Calculate time to use for restart load if we are doing run in place; different from timestamp
if do_rip_aw:
	ASSIM_TIME = int(data['ASSIM_TIME'])
	backwards = ASSIM_TIME-actual_aw #How many hours backwards should we look from current timestamp for restart to use in building state vector?
	timestamp_datetime = datetime.strptime(timestamp, "%Y%m%d_%H%M")
	delta = timedelta(hours=int(backwards))
	timestamp_restart_dt = timestamp_datetime-delta
	timestamp_restart = timestamp_restart_dt.strftime("%Y%m%d_%H%M")
else:
	timestamp_restart = None

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
