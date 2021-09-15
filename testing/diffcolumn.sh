#!/bin/bash

MY_PATH="$(jq -r ".MY_PATH" test_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" test_config.json)"
CONDA_ENV=$(jq -r ".CondaEnv" test_config.json)
TESTSTR="TESTING"

end_timestamp="$(tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP)"
end_timestamp="${end_timestamp%??}" #Clear last two characters
end_timestamp="${end_timestamp// /_}" #Replace space with underscore

#Optional arguments, or default to (1,17)
latind=${1:-1}
lonind=${2:-17}

cd ../core
source activate ${CONDA_ENV} #Activate conda environment.
python -u diff_col.py ${end_timestamp} ${TESTSTR} ${latind} ${lonind}  | tee ${MY_PATH}/${RUN_NAME}/ensemble_runs/logs/diffcol_${latind}_${lonind}.out
source deactivate
