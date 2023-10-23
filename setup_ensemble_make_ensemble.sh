#!/bin/bash

# This script sets up the ensemble run directories. Users never run this directly.

source setup_ensemble_prelims.sh #get the relevant variables.

##=======================================================================
##  Set up Ensemble run directories
##=======================================================================
printf "${thickline}CHEERIO ENSEMBLE RUN DIRECTORY CREATION${thickline}"

# Initialize (x=0 is nature run (if used), i.e. no perturbation; x=1 is ensemble member 1; etc.)
if [ "${DO_CONTROL_WITHIN_ENSEMBLE_RUNS}" = true ]; then
  x=0
else
  x=1
fi

# Create run directory for each cluster so we can apply perturbation to each
while [ $x -le $nEnsemble ];do

cd ${MY_PATH}/${RUN_NAME}/ensemble_runs

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

### Define the run directory name
name="${RUN_NAME}_${xstr}"

### Make the directory
mkdir ${name}

### Copy and point to the necessary data
cp -r ../${RUN_TEMPLATE}/*  ${name}
cd $name

### Link to GEOS-Chem executable instead of having a copy in each rundir
rm -rf gcclassic
ln -s ../../${RUN_TEMPLATE}/gcclassic .

if [ $x -eq 0 ]; then
#Switch HEMCO_Config to base/nature one.
python ${ASSIM_PATH}/core/hemco_delink_scalefactors.py $(pwd) 
else 
#Use HEMCO_Config with updated scaling factors
sed_ie "s|template_run|ensemble_runs/${name}|"  HEMCO_Config.rc #Replace template_run with this folder in HEMCO_Config
fi

if [ "${ENS_SPINUP_FROM_BC_RESTART}" = true ]; then
    sed -i -e "s|SpeciesRst|SpeciesBC|g" HEMCO_Config.rc
fi


# Link to restart file
if [ "$DO_SPINUP" = true ] ; then
  ln -s ../../spinup_run/GEOSChem.Restart.${SPINUP_END}_0000z.nc4 GEOSChem.Restart.${START_DATE}_0000z.nc4
else
  if [ "${DO_ENS_SPINUP}" = true ]; then
    ln -s $RESTART_FILE GEOSChem.Restart.${ENS_SPINUP_START}_0000z.nc4
  else 
    ln -s $RESTART_FILE GEOSChem.Restart.${START_DATE}_0000z.nc4
  fi
fi

printf "\nRun files copied and linked for ${name}.\n"


cd ..

if [ -z "$REGION" ]; then
  pygrid="$RES"
else
  pygrid="${REGION}_${MET}"
fi

### Increment
x=$[$x+1]

### Print diagnostics
printf "${thinline}CREATED: ${name}${thinline}"

done

cd ${ASSIM_PATH}
source activate $(jq -r ".CondaEnv" ens_config.json) #Activate conda environment.

#Create initial scaling factors
cd core
if [ "${DO_ENS_SPINUP}" = true ]; then
  python initialize_scaling_factors.py "${ENS_SPINUP_START}" 
else
  python initialize_scaling_factors.py "${START_DATE}" 
fi
python prep_par.py "PRODUCTION" #Figure out who should assimilate which cores. 
python setup_obs_dates.py #Create date dictionaries for rapid reference at assimilation time.

if [ "${ACTIVATE_OBSPACK}" = true ]; then #Preprocess obspack if that's what we do.
  python preprocess_obspack.py
fi

conda deactivate #Exit Conda environment

#Store current time.
if [ "${DO_ENS_SPINUP}" = true ]; then
  bash update_input_geos.sh "SPINUP" #Update input.geos to first assimilation period.
else 
  bash update_input_geos.sh "FIRST" #Update input.geos to first assimilation period.
  bash change_hemcodiag_freq.sh "ensemble" #update hemco diagnostics frequency for ensemble
fi

#Create signal file for approximate rerun case.
if [[ ${APPROXIMATE_VARON_RERUN} = "True" ]] && [[ ${DO_VARON_RERUN} = "True" ]]; then
  echo 'false' > ${MY_PATH}/${RUN_NAME}/scratch/APPOXIMATION_STAGE
fi

### Navigate back to top-level directory
cd ${MY_PATH}/${RUN_NAME}

echo "This file's existence indicates that this is the first assimilation period." > scratch/IS_FIRST

printf "${thickline}DONE CREATING ENSEMBLE MEMBER RUN DIRECTORIES${thickline}"
