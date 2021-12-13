#!/bin/bash

#Clean up scratch directory and advance timestep for entire ensemble.

if [ "${1}" = true ]; then
	MY_PATH="$(jq -r ".MY_PATH" ../testing/test_config.json)"
	RUN_NAME="$(jq -r ".RUN_NAME" ../testing/test_config.json)"
else
	MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
	RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
fi

bash update_current_time.sh ${1} #Advance timestep forward
if [ "${1}" = false ]; then
	bash update_input_geos.sh "ASSIM" #Overwrite the input.geos files.
fi

#Remove signal files
rm ${MY_PATH}/${RUN_NAME}/scratch/ALL_RUNS_COMPLETE
rm ${MY_PATH}/${RUN_NAME}/scratch/ALL_COLUMNS_FOUND
rm ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE

#Remove columns
find ${MY_PATH}/${RUN_NAME}/scratch/ -name "*.npy" -type f -delete