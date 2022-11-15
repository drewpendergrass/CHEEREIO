#!/bin/bash

#This script updates HISTORY.rc in each ensemble member 
#so that the frequency and duration under the restart collection
#matches ASSIM_TIME in ens_config.json

source activate $(jq -r ".CondaEnv" ../ens_config.json)
python update_history.py "UPDATEDURFREQ"
conda deactivate

MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
ASSIM_TIME=$(jq -r ".ASSIM_TIME" ../ens_config.json)
nEnsemble=$(jq -r ".nEnsemble" ../ens_config.json)

dcr="$(jq -r ".DO_CONTROL_RUN" ../ens_config.json)"
dcer="$(jq -r ".DO_CONTROL_WITHIN_ENSEMBLE_RUNS" ../ens_config.json)" #if true, we make a run directory without assimilation within the ensemble runs structure.

if [[ ("${dcr}" = "true" && "${dcer}" = "true") ]]; then
  DO_CONTROL_WITHIN_ENSEMBLE_RUNS=true
else
  DO_CONTROL_WITHIN_ENSEMBLE_RUNS=false
fi

cd ${MY_PATH}/${RUN_NAME}/ensemble_runs

# Initialize (x=0 is nature run (if used); x=1 is ensemble member 1; etc.)
if [ "${DO_CONTROL_WITHIN_ENSEMBLE_RUNS}" = true ]; then
  x=0
else
  x=1
fi

# Copy HISTORY.rc into each run directory
while [ $x -le $nEnsemble ];do


  ### Add zeros to string name
  if [ $x -lt 10 ]; then
      xstr="000${x}"
  elif [ $x -lt 100 ]; then
      xstr="00${x}"
  elif [ $x -lt 1000 ]; then
      xstr="0${x}"
  else
      xstr="${x}"
  fi

  ### Define the run directory name
  name="${RUN_NAME}_${xstr}"


  ### Copy and point to the necessary data
  cp -r ../template_run/HISTORY.rc  ${name}

  ### Increment
  x=$[$x+1]

done
