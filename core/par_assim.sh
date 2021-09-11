#!/bin/bash

#Assimilation script. Currently not parallelized to test flow.

#Figure out timestamp for assimilation
if [ ${1} ]; then
	MY_PATH="$(jq -r ".MY_PATH" ../testing/test_config.json)"
	RUN_NAME="$(jq -r ".RUN_NAME" ../testing/test_config.json)"
	TESTSTR="True"
else
	MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
	RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
	TESTSTR="False"
fi

end_timestamp="$(tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP)"
end_timestamp="${end_timestamp%??}" #Clear last two characters
end_timestamp="${end_timestamp// /_}" #Replace space with underscore
    
source activate $(jq -r ".CondaEnv" ../ens_config.json) #Activate conda environment.
python par_letkf.py ${end_timestamp} ${2} ${3} ${TESTSTR} >> ${MY_PATH}/${RUN_NAME}/ensemble_runs/logs/letkf_${1}_${2}.out
py_exit_status=$?
source deactivate #Exit Conda environment

#A cleanup script will check if everyone is done
printf "done" > ${MY_PATH}/${RUN_NAME}/scratch/complete_subprocess_${1}_${2}.CHECK

#If python does not exit with exit code one, make file that will break loop
if [ $py_exit_status != 0 ]; then
	printf "Python assimilation script exited without code 0 in ensemble ${1} and core ${2} \n" > ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS #This file's presence will break loop
fi