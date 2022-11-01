#!/bin/bash
eval "$(conda shell.bash hook)"

#Figure out timestamp for assimilation
MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
CONDA_ENV=$(jq -r ".CondaEnv" ../ens_config.json)

end_timestamp="$(tail -n 1 ${MY_PATH}/${RUN_NAME}/scratch/INPUT_GEOS_TEMP)"
end_timestamp="${end_timestamp%??}" #Clear last two characters
end_timestamp="${end_timestamp// /_}" #Replace space with underscore
    
source activate ${CONDA_ENV} #Activate conda environment.
python -u par_letkf.py ${end_timestamp} ${1} ${2} ${3} ${4} >> ${MY_PATH}/${RUN_NAME}/ensemble_runs/logs/letkf_${1}_${2}.out
py_exit_status=$?
conda deactivate #Exit Conda environment

#If python does not exit with exit code one, make file that will break loop
if [ $py_exit_status != 0 ]; then
	printf "Python assimilation script exited without code 0 in ensemble ${1} and core ${2} \n" > ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS #This file's presence will break loop
fi