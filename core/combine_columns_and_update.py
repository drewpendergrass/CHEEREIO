from glob import glob
import toolbox as tx 
import sys 
import letkf_utils as lu
import time
import os

timestamp = str(sys.argv[1]) #Time to assimilate. Expected in form YYYYMMDD_HHMM, UTC time.

cmdarg = str(sys.argv[2])
if cmdarg=="TESTING":
	testing=True
else:
	testing=False

data = tx.getSpeciesConfig(testing=testing)

SaveDOFS = data["SaveDOFS"] == "True"

dateval = timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]

print(f'One core is gathering columns to overwrite at time {dateval}.')
start = time.time()
wrapper = lu.GT_Container(timestamp,testing)
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
