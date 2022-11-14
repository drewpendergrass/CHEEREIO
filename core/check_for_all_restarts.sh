#!/bin/bash

#This script checks to see if all restarts are written, meaning ensemble is complete.

MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
nEnsemble=$(jq -r ".nEnsemble" ../ens_config.json)

dcr=$(jq -r ".DO_CONTROL_RUN" ens_config.json)
dcer=$(jq -r ".DO_CONTROL_WITHIN_ENSEMBLE_RUNS" ens_config.json) #if true, we make a run directory without assimilation within the ensemble runs structure.

if [[ ("${dcr}" = "true" && "${dcer}" = "true") ]]; then
  DO_CONTROL_WITHIN_ENSEMBLE_RUNS=true
else
  DO_CONTROL_WITHIN_ENSEMBLE_RUNS=false
fi

end_timestamp="$(tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP)"
end_timestamp="${end_timestamp%??}" #Clear last two characters
end_timestamp="${end_timestamp// /_}" #Replace space with underscore
rst_filename="GEOSChem.Restart.${end_timestamp}z.nc4"
 
# Initialize (x=0 is nature run (if used), i.e. no perturbation; x=1 is ensemble member 1; etc.)
if [ "${DO_CONTROL_WITHIN_ENSEMBLE_RUNS}" = true ]; then
    x=0
    bonus=1
else
    x=1
    bonus=0
fi

counter=0 #Count the restarts
# Check if file_to_check exists.
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
         filetocheck="${MY_PATH}/${RUN_NAME}/ensemble_runs/${name}/${rst_filename}"
         if [ -f $filetocheck ]; then
            counter=$[$counter+1] #We found one!
         fi
         x=$[$x+1]
    done

numtocheck=$((nEnsemble + bonus))

if [ $counter -eq $numtocheck ]; then
	printf "\nAll ensemble folders contain file ${rst_filename}.\n"
	printf "Proceeding to data assimilation phase.\n"
	printf "done" > ${MY_PATH}/${RUN_NAME}/scratch/ALL_RUNS_COMPLETE #This file's presence will take us to assimilation phase
fi

