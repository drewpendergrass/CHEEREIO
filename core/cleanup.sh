#!/bin/bash

#Clean up scratch directory and advance timestep for entire ensemble.

MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
APPROXIMATE_VARON_RERUN=$(jq -r ".APPROXIMATE_VARON_RERUN" ../ens_config.json)
DO_VARON_RERUN=$(jq -r ".DO_VARON_RERUN" ../ens_config.json)

#If we are doing an "approximate rerun" simulation, and we just did an LETKF (approximation = false), then next run will rerun previous assimilation window.
#To ensure proper job control, delete all restarts from end of this window. 
if [[ ${APPROXIMATE_VARON_RERUN} = "True" ]] && [[ ${DO_VARON_RERUN} = "True" ]]; then 
	do_approx="$(tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/APPOXIMATION_STAGE)"
	if [[ ${do_approx} = "false" ]]; then
		#Previous simulation was an LETKF (approx = false) so clear restarts at end of assimilation window (will not change)
		end_timestamp="$(tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP)"
		end_timestamp="${end_timestamp%??}" #Clear last two characters
		end_timestamp="${end_timestamp// /_}" #Replace space with underscore
		rst_filename="GEOSChem.Restart.${end_timestamp}z.nc4"
	fi
fi

bash update_input_geos.sh ${1} #Overwrite the input.geos files.


rm ${MY_PATH}/${RUN_NAME}/ensemble_runs/${RUN_NAME}_*/${rst_filename}

#Remove signal files
rm ${MY_PATH}/${RUN_NAME}/scratch/ALL_RUNS_COMPLETE
rm ${MY_PATH}/${RUN_NAME}/scratch/ALL_COLUMNS_FOUND
rm ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE

#Remove columns
find ${MY_PATH}/${RUN_NAME}/scratch/ -name "*.npy" -type f -delete

if [[ ${do_approx} = "false" ]]; then
	echo 'woohoo'
fi
