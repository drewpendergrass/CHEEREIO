import sys 
import letkf_utils as lu
import toolbox as tx 

#Called like so: python diff_col.py TIMESTAMP TESTSTRING LATIND LONIND

timestamp = str(sys.argv[1]) #Time to assimilate. Expected in form YYYYMMDD_HHMM, UTC time.

cmdarg = str(sys.argv[2])
if cmdarg=="TESTING":
	testbool = True
else:
	testbool = False

data = tx.getSpeciesConfig(testing=testbool)
latind = int(sys.argv[3])
lonind = int(sys.argv[4])
dateval = timestamp[0:4]+'-'+timestamp[4:6]+'-'+timestamp[6:8]

control_conc = data['CONTROL_VECTOR_CONC']
control_emis = list(data['CONTROL_VECTOR_EMIS'].keys())

wrapper = lu.GT_Container(timestamp,testbool,constructStateVecs=True)

print(f'AUTOTESTING CONTROL VECTOR SPECIES AT {(latind,lonind)}')
for species in control_conc:
	wrapper.compareSpeciesConc(species,latind,lonind)

print(f'AUTOTESTING CONTROL VECTOR EMISSIONS AT {(latind,lonind)}')
for species in control_emis:
	wrapper.compareSpeciesEmis(species,latind,lonind)
