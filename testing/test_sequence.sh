#!/bin/bash

#If command arg is SETUP, setup the testdir; otherwise, just clean up

if [ $1 -eq "SETUP" ]; then
	bash setup_test_dir.sh
else
	RUN_NAME="$(jq -r ".RUN_NAME" test_config.json)"
	MY_PATH="$(jq -r ".MY_PATH" test_config.json)"
	rm -rf ${MY_PATH}/${RUN_NAME}/ensemble_runs/logs/*
	rm -rf ${MY_PATH}/${RUN_NAME}/scratch/*.CHECK
	rm -rf ${MY_PATH}/${RUN_NAME}/scratch/*.npy
	rm ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS
fi

bash autotest_one_core.sh