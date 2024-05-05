#!/bin/bash

# This creates the template run directory for versions 13 and higher.

mkdir -p ${MY_PATH}/${RUN_NAME}
cd ${MY_PATH}/${RUN_NAME}
mkdir -p ensemble_runs
mkdir -p ensemble_runs/logs
mkdir -p scratch
mkdir -p postprocess 
mkdir -p postprocess/bigy

#Copy CHEEREIO code
mkdir -p CHEEREIO
cp -r ${ASSIM_PATH}/core CHEEREIO/
cp -r ${ASSIM_PATH}/extensions CHEEREIO/
cp -r ${ASSIM_PATH}/postprocess CHEEREIO/
cp -r ${ASSIM_PATH}/tests CHEEREIO/
cp -r ${ASSIM_PATH}/templates CHEEREIO/
cp -r ${ASSIM_PATH}/environments CHEEREIO/
cp ${ASSIM_PATH}/*.json CHEEREIO/
cp ${ASSIM_PATH}/*.sh CHEEREIO/

if [ "${SaveDOFS}" = "True" ]; then
  mkdir -p ensemble_runs/logs/dofs_scratch
  mkdir -p ensemble_runs/logs/dofs_complete
fi

echo "GC-CHEERIO uses this directory to save out intermediate data and track its internal state. Modifying contents of this folder can lead to model failure." > scratch/README

cp ${MY_PATH}/${RUN_NAME}/CHEEREIO/templates/run_ensemble_simulations.sh ensemble_runs/
cp ${MY_PATH}/${RUN_NAME}/CHEEREIO/templates/run_ens.sh ensemble_runs/

sed -i -e "s:{RunName}:${RUN_NAME}:g" \
       -e "s:{NumCores}:${NumCores}:g" \
       -e "s:{Partition}:${Partition}:g" \
       -e "s:{Memory}:${Memory}:g" \
       -e "s:{WallTime}:${WallTime}:g" \
       -e "s:{MaxPar}:${MaxPar}:g" \
       -e "s:{ASSIM}:${MY_PATH}/${RUN_NAME}/CHEEREIO:g" ensemble_runs/run_ensemble_simulations.sh

if [ "${DO_CONTROL_WITHIN_ENSEMBLE_RUNS}" = true ]; then
  sed -i -e "s:{START}:0:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ens.sh
else
  sed -i -e "s:{START}:1:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ens.sh
fi

if [ "${DO_ENS_SPINUP}" = true ]; then
  cp ${MY_PATH}/${RUN_NAME}/CHEEREIO/templates/run_ensemble_spinup_simulations.sh ensemble_runs/
  cp ${MY_PATH}/${RUN_NAME}/CHEEREIO/templates/run_ensspin.sh ensemble_runs/

  sed -i -e "s:{RunName}:${RUN_NAME}:g" \
         -e "s:{NumCores}:${NumCores}:g" \
         -e "s:{Partition}:${Partition}:g" \
         -e "s:{Memory}:${EnsCtrlSpinupMemory}:g" \
         -e "s:{WallTime}:${EnsSpinupWallTime}:g" \
         -e "s:{MaxPar}:${MaxPar}:g" \
         -e "s:{ASSIM}:${MY_PATH}/${RUN_NAME}/CHEEREIO:g" ensemble_runs/run_ensemble_spinup_simulations.sh

  if [ "${DO_CONTROL_WITHIN_ENSEMBLE_RUNS}" = true ]; then
    sed -i -e "s:{START}:0:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ensspin.sh
  else
    sed -i -e "s:{START}:1:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ensspin.sh
  fi
fi

mkdir -p ${RUN_TEMPLATE}
