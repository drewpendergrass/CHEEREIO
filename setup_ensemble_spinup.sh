#!/bin/bash

# This script sets up the global spinup run directory. Users never run this directly.

source setup_ensemble_prelims.sh #get the relevant variables.

##=======================================================================
##  Set up spinup run directory
##=======================================================================

printf "${thickline}CHEERIO SPINUP RUN DIRECTORY CREATION${thickline}"

cd ${MY_PATH}/${RUN_NAME}

### Define the run name
spinup_name="${RUN_NAME}_Spinup"

### Make the directory
runDir="spinup_run"
mkdir -p ${runDir}

### Copy and point to the necessary data
cp -r ${RUN_TEMPLATE}/*  ${runDir}
cp -RLv ${MY_PATH}/${RUN_NAME}/CHEEREIO/templates/ensemble_run.template ${runDir}
cd $runDir

### Link to GEOS-Chem executable instead of having a copy in each rundir
rm -rf gcclassic
ln -s ../${RUN_TEMPLATE}/gcclassic .
# Link to restart file
ln -s $RESTART_FILE GEOSChem.Restart.${SPINUP_START}_0000z.nc4

#Switch HEMCO_Config to base/nature one.
python ${MY_PATH}/${RUN_NAME}/CHEEREIO/core/hemco_delink_scalefactors.py $(pwd) 

gc_major_version="${GC_VERSION:0:2}"
if [ $gc_major_version = "13" ]; then
filename='input.geos'
elif [ $gc_major_version = "14" ]; then
filename='geoschem_config.yml'
else
printf "\n ERROR: CANNOT WORK WITH GEOS-CHEM VERSION SUPPLIED; must be 13 or later."
exit 1
fi

### Update settings in input.geos
sed -i -e "s:{DATE1}:${SPINUP_START}:g" \
       -e "s:{DATE2}:${SPINUP_END}:g" \
       -e "s:{TIME1}:000000:g" \
       -e "s:{TIME2}:000000:g" ${filename}

### Create run script from template
sed -e "s:namename:${spinup_name}:g" \
    -e "s:{NumCores}:${NumCores}:g" \
    -e "s:{Partition}:${Partition}:g" \
    -e "s:{Memory}:${Memory}:g" \
    -e "s:{SpinupWallTime}:${SpinupWallTime}:g" ensemble_run.template > ${spinup_name}.run #change
chmod 755 ${spinup_name}.run
rm -f ensemble_run.template

### Print diagnostics
printf "${thinline}CREATED: ${runDir}${thinline}"

### Navigate back to top-level directory
cd ${MY_PATH}/${RUN_NAME}
