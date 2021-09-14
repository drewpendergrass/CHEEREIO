#!/bin/bash

#If command arg is SETUP, setup the testdir; otherwise, just clean up

if [ "$1" = "SETUP" ]; then
	printf "Reinstalling test directory...\n"
	bash setup_test_dir.sh
else
	printf "Cleaning test directory...\n"
	RUN_NAME="$(jq -r ".RUN_NAME" test_config.json)"
	MY_PATH="$(jq -r ".MY_PATH" test_config.json)"
	rm -rf ${MY_PATH}/${RUN_NAME}/ensemble_runs/logs/*
	rm -rf ${MY_PATH}/${RUN_NAME}/scratch/*.CHECK
	rm -rf ${MY_PATH}/${RUN_NAME}/scratch/*.npy
	rm ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS
fi

printf "Test directory ready. Begin autotest for one core.\n"
bash autotest_one_core.sh
printf "Autotest for one core complete.\n"
