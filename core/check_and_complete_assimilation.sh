#!/bin/bash

#Check if all columns are saved. If so, combine and overwrite

if [ ${1} ]; then
	MY_PATH="$(jq -r ".MY_PATH" ../testing/test_config.json)"
	RUN_NAME="$(jq -r ".RUN_NAME" ../testing/test_config.json)"
	TESTSTR='TESTING'
else
	MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
	RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
	TESTSTR='PRODUCTION'
fi

python check_for_all_columns.py ${TESTSTR}

if [ -f ${MY_PATH}/${RUN_NAME}/scratch/ALL_COLUMNS_FOUND ]; then
	python combine_columns_and_update.py ${TESTSTR}
	echo 'Done' > ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE
fi
