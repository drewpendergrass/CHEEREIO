#!/bin/bash

#Assimilation script. Currently not parallelized to test flow.

#Figure out timestamp for assimilation
if [ "${1}" = true ]; then
	MY_PATH="$(jq -r ".MY_PATH" ../testing/test_config.json)"
	RUN_NAME="$(jq -r ".RUN_NAME" ../testing/test_config.json)"
	CONDA_ENV=$(jq -r ".CondaEnv" ../testing/test_config.json)
	TESTSTR="True"
else
	MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
	RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
	CONDA_ENV=$(jq -r ".CondaEnv" ../ens_config.json)
	TESTSTR="False"
fi

end_timestamp="$(tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP)"
end_timestamp="${end_timestamp%??}" #Clear last two characters
end_timestamp="${end_timestamp// /_}" #Replace space with underscore
    
source activate ${CONDA_ENV} #Activate conda environment.
#if testing, output to console as well as file
if [ ${1} ]; then
	python -u par_letkf.py ${end_timestamp} ${2} ${3} ${TESTSTR} | tee ${MY_PATH}/${RUN_NAME}/ensemble_runs/logs/letkf_${2}_${3}.out
else
	python -u par_letkf.py ${end_timestamp} ${2} ${3} ${TESTSTR} >> ${MY_PATH}/${RUN_NAME}/ensemble_runs/logs/letkf_${2}_${3}.out
fi
py_exit_status=$?
source deactivate #Exit Conda environment

#If python does not exit with exit code one, make file that will break loop
if [ $py_exit_status != 0 ]; then
	printf "Python assimilation script exited without code 0 in ensemble ${2} and core ${3} \n" > ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS #This file's presence will break loop
fi