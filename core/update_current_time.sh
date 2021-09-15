#!/bin/bash

#This script updates CURRENT_DATE_TIME in the scratch folder
#based on the end of the current input.geos file.

if [ ${1} ]; then
	MY_PATH="$(jq -r ".MY_PATH" ../testing/test_config.json)"
	RUN_NAME="$(jq -r ".RUN_NAME" ../testing/test_config.json)"
else
	MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
	RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
fi

tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP > ${MY_PATH}/${RUN_NAME}/scratch/CURRENT_DATE_TIME