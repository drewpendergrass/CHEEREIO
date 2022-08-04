#!/bin/bash

#This script updates HEMCO_Config.rc in each ensemble member 
#so that the frequency of diagnostic write-out 
#matches ASSIM_TIME in ens_config.json


MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
ASSIM_TIME=$(jq -r ".ASSIM_TIME" ../ens_config.json)
nEnsemble=$(jq -r ".nEnsemble" ../ens_config.json)
SIMULATE_NATURE=$(jq -r ".SIMULATE_NATURE" ../ens_config.json)

(( hours=ASSIM_TIME%24 ))
(( days=ASSIM_TIME/24 ))

### Add zeros to string name
if [ $hours -lt 10 ]; then
  hoursstr="0${hours}"
else
  hoursstr="${hours}"
fi

if [ $days -lt 10 ]; then
  daysstr="0${days}"
else
  daysstr="${days}"
fi

timestr="000000${daysstr} ${hoursstr}0000"
freqstr="DiagnFreq\:                   ${timestr}"

cd ${MY_PATH}/${RUN_NAME}/ensemble_runs

# Initialize (x=0 is nature run (if used); x=1 is ensemble member 1; etc.)
if [ "${SIMULATE_NATURE}" = true ]; then
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


  ### Modify the desired line
  sed -i "/DiagnFreq/c\\${freqstr}" ${name}/HEMCO_Config.rc

  ### Increment
  x=$[$x+1]

done
