#!/bin/bash

#So you've gotten a dreaded KILL_ENS file in scratch, eh? First, cancel the job array. 
#Then run this script to clear out signal files, restarts, columns, etc.
#Then you can resubmit the job array, after making whatever fix you need to make (if any).
#FYI: Sometimes CHEEREIO just fails because of cluster issues.

MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"

#Remove signal files
rm ${MY_PATH}/${RUN_NAME}/scratch/ALL_RUNS_COMPLETE
rm ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS
rm ${MY_PATH}/${RUN_NAME}/scratch/ALL_COLUMNS_FOUND
rm ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE
rm ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS

#Remove columns
find ${MY_PATH}/${RUN_NAME}/scratch/ -name "*.npy" -type f -delete

#Remove restart from end of run; this is because it can mess up job control
#Better to just rerun the latest GC iteration.
end_timestamp="$(tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP)"
end_timestamp="${end_timestamp%??}" #Clear last two characters
end_timestamp="${end_timestamp// /_}" #Replace space with underscore
rst_filename="GEOSChem.Restart.${end_timestamp}z.nc4"

rm ${MY_PATH}/${RUN_NAME}/ensemble_runs/${RUN_NAME}_*/${rst_filename}