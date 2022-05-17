import letkf_utils as lu
import sys
import time

timestamp = str(sys.argv[1]) #Time to assimilate. Expected in form YYYYMMDD_HHMM, UTC time.
ensnum = int(sys.argv[2])
corenum = int(sys.argv[3])

dateval = timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]

print(f'Core ({ensnum},{corenum}) is gathering ensemble at time {dateval}.')
start = time.time()
print(f'Assimilator call: lu.Assimilator({timestamp},{ensnum},{corenum})')
#assimilator = lu.Assimilator('20190108_0000',2,1)
assimilator = lu.Assimilator(timestamp,ensnum,corenum)
end = time.time()
print(f'Core ({ensnum},{corenum}) gathered ensemble in {end - start} seconds. Begin LETKF procedure.')
start = time.time()
assimilator.LETKF()
end = time.time()
print(f'Core ({ensnum},{corenum}) completed computation for {dateval} and saved columns in {end - start} seconds.')
print('-------------------END LETKF-------------------')
