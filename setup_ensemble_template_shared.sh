#!/bin/bash

# This creates the template run directory for version 13.

mkdir -p ${MY_PATH}/${RUN_NAME}
cd ${MY_PATH}/${RUN_NAME}
mkdir -p ensemble_runs
mkdir -p ensemble_runs/logs
mkdir -p scratch
mkdir -p postprocess

if [ "${SaveDOFS}" = "True" ]; then
  mkdir -p ensemble_runs/logs/dofs_scratch
  mkdir -p ensemble_runs/logs/dofs_complete
fi

echo "GC-CHEERIO uses this directory to save out intermediate data and track its internal state. Modifying contents of this folder can lead to model failure." > scratch/README

cp ${ASSIM_PATH}/templates/run_ensemble_simulations.sh ensemble_runs/
cp ${ASSIM_PATH}/templates/run_ens.sh ensemble_runs/

sed -i -e "s:{RunName}:${RUN_NAME}:g" \
       -e "s:{NumCores}:${NumCores}:g" \
       -e "s:{Partition}:${Partition}:g" \
       -e "s:{Memory}:${Memory}:g" \
       -e "s:{WallTime}:${WallTime}:g" \
       -e "s:{MaxPar}:${MaxPar}:g" \
       -e "s:{ASSIM}:${ASSIM_PATH}:g" ensemble_runs/run_ensemble_simulations.sh

if [ "${DO_CONTROL_WITHIN_ENSEMBLE_RUNS}" = true ]; then
  sed -i -e "s:{START}:0:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ens.sh
else
  sed -i -e "s:{START}:1:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ens.sh
fi

if [ "${DO_ENS_SPINUP}" = true ]; then
  cp ${ASSIM_PATH}/templates/run_ensemble_spinup_simulations.sh ensemble_runs/
  cp ${ASSIM_PATH}/templates/run_ensspin.sh ensemble_runs/

  sed -i -e "s:{RunName}:${RUN_NAME}:g" \
         -e "s:{NumCores}:${NumCores}:g" \
         -e "s:{Partition}:${Partition}:g" \
         -e "s:{Memory}:${EnsCtrlSpinupMemory}:g" \
         -e "s:{WallTime}:${EnsSpinupWallTime}:g" \
         -e "s:{MaxPar}:${MaxPar}:g" \
         -e "s:{ASSIM}:${ASSIM_PATH}:g" ensemble_runs/run_ensemble_spinup_simulations.sh

  if [ "${DO_CONTROL_WITHIN_ENSEMBLE_RUNS}" = true ]; then
    sed -i -e "s:{START}:0:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ensspin.sh
  else
    sed -i -e "s:{START}:1:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ensspin.sh
  fi
fi

mkdir -p ${RUN_TEMPLATE}
