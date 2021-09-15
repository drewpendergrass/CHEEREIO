from glob import glob
import toolbox as tx 
import sys 
import letkf_utils as lu
import time

timestamp = str(sys.argv[1]) #Time to assimilate. Expected in form YYYYMMDD_HHMM, UTC time.

cmdarg = str(sys.argv[2])
if cmdarg=="TESTING":
	data = tx.getSpeciesConfig(testing=True)
else:
	data = tx.getSpeciesConfig(testing=False)

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
end = time.time()
print(f'Saved updated restarts and emissions in {end - start} seconds. We can cleanup and resume GC now!')
