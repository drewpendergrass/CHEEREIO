#!/bin/bash

#This script updates input_geos in each ensemble member 
#based on the time stored in the scratch folder.

source activate $(jq -r ".CondaEnv" ../ens_config.json)
python advance_timestep.py "${1}"

MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
ASSIM_TIME=$(jq -r ".ASSIM_TIME" ../ens_config.json)
nEnsemble=$(jq -r ".nEnsemble" ../ens_config.json)
GC_VERSION="$(jq -r ".GC_VERSION" ../ens_config.json)"
ACTIVATE_OBSPACK="$(jq -r ".ACTIVATE_OBSPACK" ../ens_config.json)"
gc_major_version="${GC_VERSION:0:2}"

  if [ $gc_major_version = "13" ]; then
    filename='input.geos'
  elif [ $gc_major_version = "14" ]; then
    filename='geoschem_config.yml'
    fi

dcr="$(jq -r ".DO_CONTROL_RUN" ../ens_config.json)"
dcer="$(jq -r ".DO_CONTROL_WITHIN_ENSEMBLE_RUNS" ../ens_config.json)" #if true, we make a run directory without assimilation within the ensemble runs structure.

if [[ ("${dcr}" = "true" && "${dcer}" = "true") ]]; then
  DO_CONTROL_WITHIN_ENSEMBLE_RUNS=true
else
  DO_CONTROL_WITHIN_ENSEMBLE_RUNS=false
fi

counter=1
while IFS=' ' read -r date time
do 
    # Initialize (x=0 is nature run (if used), i.e. no perturbation; x=1 is ensemble member 1; etc.)
    if [ "$DO_CONTROL_WITHIN_ENSEMBLE_RUNS" = true ]; then
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
            cp ${MY_PATH}/${RUN_NAME}/template_run/${filename} ${MY_PATH}/${RUN_NAME}/ensemble_runs/${name}/${filename}
         fi
         sed -i -e "s:{DATE${counter}}:${date}:g" \
                -e "s:{TIME${counter}}:${time}:g" ${MY_PATH}/${RUN_NAME}/ensemble_runs/${name}/${filename}
         x=$[$x+1]
         if [ "$ACTIVATE_OBSPACK" = true ]; then
            if [ "${1}" = "FIRST" ]; then
                python obspack_switch.py "${MY_PATH}/${RUN_NAME}/ensemble_runs/${name}/${filename}"
            fi
         fi
       done
       #Increment so we do end time
       counter=$[$counter+1]
done <${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP

conda deactivate

printf "\n${filename} updated for ensemble.\n"
