from Assimilator import Assimilator
import sys
import time
import settings_interface as si 

data = si.getSpeciesConfig()
path_to_scratch = f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch"
timestamp = str(sys.argv[1]) #Time to assimilate. Expected in form YYYYMMDD_HHMM, UTC time.
ensnum = int(sys.argv[2])
corenum = int(sys.argv[3])
just_scale = (str(sys.argv[4])=='true') and (data['SIMPLE_SCALE_FOR_FIRST_ASSIM_PERIOD'] == "true") #If this is the first run AND user specifies, we'll just scale (no assimilation)
if just_scale:
	label_str = 'scaling'
else:
	label_str = 'LETKF'
dateval = timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]

#If we are just scaling, only ens member 1 and core 1 will do assimilation
if (not just_scale) or (just_scale and (ensnum == 1) and (corenum==1) ):
	print(f'Core ({ensnum},{corenum}) is gathering ensemble at time {dateval}.')
	start = time.time()
	print(f'Assimilator call: Assimilator({timestamp},{ensnum},{corenum})')
	#a = Assimilator('20190108_0000',2,1)
	a = Assimilator(timestamp,ensnum,corenum)
	end = time.time()
	print(f'Core ({ensnum},{corenum}) gathered ensemble in {end - start} seconds. Begin {label_str} procedure.')
	start = time.time()
	if just_scale:
		a.scaleRestarts()
		with open(f"{path_to_scratch}/ASSIMILATION_COMPLETE", "w") as f:
			f.write("Done.\n") #If so, save flag file to ensemble folder
	else:
		a.LETKF()
	end = time.time()
	print(f'Core ({ensnum},{corenum}) completed computation for {dateval} and saved columns in {end - start} seconds.')
else:
	print(f'Just scaling restarts. Core ({ensnum},{corenum}) is not needed and will hang.')
print('-------------------END LETKF-------------------')
