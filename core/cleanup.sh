#!/bin/bash

#Clean up scratch directory and advance timestep for entire ensemble.

MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"

bash update_current_time.sh #Advance timestep forward
bash update_input_geos.sh #Overwrite the input.geos files.

#Remove signal files
rm ${MY_PATH}/${RUN_NAME}/scratch/ALL_RUNS_COMPLETE
rm ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE