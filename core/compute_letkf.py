import observation_operators as obs
import letkf_utils as lu
import sys
import time

timestamp = str(sys.argv[1]) #Time to assimilate. Expected in form YYYYMMDD_HHMM, UTC time.

dateval = timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]

print('------------------BEGIN LETKF------------------')
print(f'Gathering ensemble at time {dateval}.')
start = time.time()
assimilator = lu.Assimilator(timestamp)
end = time.time()
print(f'Ensemble gathered in {end - start} seconds. Begin LETKF procedure.')
start = time.time()
assimilator.LETKF()
end = time.time()
print(f'Computation complete in {end - start} seconds. Save analysis ensemble netCDFs.')
start = time.time()
assimilator.updateRestartsAndScalingFactors()
assimilator.saveRestartsAndScalingFactors()
end = time.time()
print(f'Saving completed in {end - start} seconds.')
print(f'Data for {dateval} processed.')
print('-------------------END LETKF-------------------')
