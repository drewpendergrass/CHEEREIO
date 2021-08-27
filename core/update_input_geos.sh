#!/bin/bash

#This script updates input_geos in each ensemble member 
#based on the time stored in the scratch folder.

python advance_timestep.py

MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
ASSIM_TIME=$(jq -r ".ASSIM_TIME" ../ens_config.json)
nEnsemble=$(jq -r ".nEnsemble" ../ens_config.json)
SIMULATE_NATURE=$(jq -r ".SIMULATE_NATURE" ../ens_config.json)


counter=1
while IFS=' ' read -r date time
do 
    # Initialize (x=0 is nature run (if used), i.e. no perturbation; x=1 is ensemble member 1; etc.)
    if [ $SIMULATE_NATURE ]; then
           x=0
    else
           x=1
    fi
    # Fix input geos in each ensemble folder
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
         name="${RUN_NAME}_${xstr}"
         if [ $counter -eq 1 ]; then
            cp ${MY_PATH}/${RUN_NAME}/template_run/input.geos ${MY_PATH}/${RUN_NAME}/ensemble_runs/${name}/input.geos
         fi
         sed -i -e "s:{DATE${counter}}:${date}:g" \
                -e "s:{TIME${counter}}:${time}:g" ${MY_PATH}/${RUN_NAME}/ensemble_runs/${name}/input.geos
         x=$[$x+1]
       done
       #Increment so we do end time
       counter=$[$counter+1]
done <${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP

printf "\ninput.geos updated for ensemble.\n"
