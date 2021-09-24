#!/bin/bash

printf '\n*****Attempting to create template HEMCO_Config file*****\n'
printf 'This process can fail so check HEMCO_Config in the template run directory\n'
printf 'and add relevant tags BEFORE copying into ensemble run directories.\n'
printf 'See documentation for more information.\n'

RUN_NAME="$(jq -r ".RUN_NAME" ../ens_config.json)"
MY_PATH="$(jq -r ".MY_PATH" ../ens_config.json)"
CONDA_ENV="$(jq -r ".CondaEnv" ../ens_config.json)"

#Save default HEMCO_Config for nature/spinup run (no scaling factors read)
cp ${MY_PATH}/${RUN_NAME}/template_run/HEMCO_Config.rc ${MY_PATH}/${RUN_NAME}/template_run/HEMCO_Config_SPINUP_NATURE_TEMPLATE.rc
source activate ${CONDA_ENV}
python hemco_config_updater.py
conda deactivate