import toolbox as tx
from datetime import datetime,timedelta
import sys

teststr = str(sys.argv[1])
periodstr = str(sys.argv[2])
if teststr=="TESTING":
	spc_config = tx.getSpeciesConfig(testing=True)
else:
	spc_config = tx.getSpeciesConfig(testing=False)

parent_dir = f"{spc_config['MY_PATH']}/{spc_config['RUN_NAME']}"
ens_dir = f"{parent_dir}/ensemble_runs"
ASSIM_TIME = spc_config['ASSIM_TIME']
ENS_END_DATE = spc_config['END_DATE']
ASSIM_START_DATE = spc_config['ASSIM_START_DATE']
ENS_END_DATE_datetime = datetime.strptime(ENS_END_DATE, "%Y%m%d")

with open(f"{parent_dir}/scratch/CURRENT_DATE_TIME") as f:
    start_string = f.readlines()[0].rstrip()

start_datetime = datetime.strptime(start_string, "%Y%m%d %H%M%S")
if periodstr=="FIRST":
	end_string = f"{ASSIM_START_DATE} 000000"
else:
	delta = timedelta(hours=int(ASSIM_TIME))
	end_datetime = start_datetime+delta
	end_string = end_datetime.strftime("%Y%m%d %H%M%S")

#Check if we have finished the total ensemble run
if start_datetime >= ENS_END_DATE_datetime:
	with open(f"{parent_dir}/scratch/ENSEMBLE_COMPLETE", "w") as j:
		f.write("Ensemble completed; delete this file if you want to re-run.\n") #If so, save flag file to ensemble folder

with open(f"{parent_dir}/scratch/INPUT_GEOS_TEMP", "w") as f:
	f.write(f"{start_string}\n")
	f.write(f"{end_string}\n")
	f.close()

