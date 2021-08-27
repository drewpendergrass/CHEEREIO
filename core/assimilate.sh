#!/bin/bash

#Assimilation script. Currently not parallelized to test flow.

#Figure out timestamp for assimilation
MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"

end_timestamp="$(tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP)"
end_timestamp="${end_timestamp%??}" #Clear last two characters
end_timestamp="${end_timestamp// /_}" #Replace space with underscore
    
source activate $(jq -r ".CondaEnv" ../ens_config.json) #Activate conda environment.
python compute_letkf.py ${end_timestamp} >> ${MY_PATH}/${RUN_NAME}/ensemble_runs/logs/letkf.out
py_exit_status=$?
source deactivate #Exit Conda environment

printf "done" > ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE #This file's presence will take us to cleanup phase

#If python does not exit with exit code one, make file that will break loop
if [ $py_exit_status != 0 ]; then
	printf "Python assimilation script exited without code 0\n" > ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS #This file's presence will break loop
fi