import settings_interface as si 
from glob import glob

data = si.getSpeciesConfig()
path_to_ensemble = f"{data['MY_PATH']}/{data['RUN_NAME']}/ensemble_runs"
path_to_scratch = f"{data['MY_PATH']}/{data['RUN_NAME']}/scratch"

subdirs = glob(f"{path_to_ensemble}/*/")
subdirs.remove(f"{path_to_ensemble}/logs/")

expected_count = len(subdirs)
finished_count = 0

for subdir in subdirs:
	with open(f"{subdir}GC.log") as f:
		for line in f:
			pass
		last_line = line
	if last_line.startswith('**************   E N D'):
		finished_count+=1 #We got one!

if expected_count == finished_count:
	with open(f"{path_to_scratch}/ALL_RUNS_COMPLETE", "w") as f:
		f.write("All ensemble runs are complete according to their log file.\n")
		f.write("Proceeding to data assimilation phase.\n")
		f.write("done\n")
		f.close()
