#!/bin/bash

#Check if all columns are saved. If so, combine and overwrite

MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
CONDA_ENV=$(jq -r ".CondaEnv" ../ens_config.json)
SIMPLE_SCALE_FOR_FIRST_ASSIM_PERIOD="$(jq -r ".SIMPLE_SCALE_FOR_FIRST_ASSIM_PERIOD" ../ens_config.json)"
TESTSTR='PRODUCTION'

end_timestamp="$(tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP)"
end_timestamp="${end_timestamp%??}" #Clear last two characters
end_timestamp="${end_timestamp// /_}" #Replace space with underscore

#If we are not doing simple scaling, we need to do the actual work to combine columns.
if [ "${1}" = false ]; then

	conda run -n $(jq -r ".CondaEnv" ../ens_config.json) python check_for_all_columns.py

	if [ -f ${MY_PATH}/${RUN_NAME}/scratch/ALL_COLUMNS_FOUND ]; then
		conda run -n $(jq -r ".CondaEnv" ../ens_config.json) python -u combine_columns_and_update.py ${end_timestamp} >> ${MY_PATH}/${RUN_NAME}/ensemble_runs/logs/letkf_master.out
		py_exit_status=$?
		if [ $py_exit_status != 0 ]; then
			printf "Python combine columns script exited without code 0 \n" > ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS #This file's presence will break loop
		else
			echo 'Done' > ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE
		fi
	fi

fi

