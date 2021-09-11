#!/bin/bash

#If command arg is SETUP, setup the testdir; otherwise, just clean up

if [ $1 -eq "SETUP" ]; then
	bash setup_test_dir.sh
else
	RUN_NAME="$(jq -r ".RUN_NAME" test_config.json)"
	MY_PATH="$(jq -r ".MY_PATH" test_config.json)"
	rm -rf ${RUN_NAME}/${MY_PATH}/ensemble_runs/logs/*
	rm -rf ${RUN_NAME}/${MY_PATH}/scratch/*.CHECK
	rm -rf ${RUN_NAME}/${MY_PATH}/scratch/*.npy
	rm ${RUN_NAME}/${MY_PATH}/KILL_ENS
fi

bash autotest_one_core.sh