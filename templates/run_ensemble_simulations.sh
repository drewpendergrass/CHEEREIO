#!/bin/bash

#SBATCH -J {RunName}
#SBATCH -c {NumCores}
#SBATCH -N 1
#SBATCH -p {Partition}
#SBATCH --mem {Memory}
#SBATCH -t {WallTime}
#SBATCH -o logs/ensemble_slurm_%j.out    # File to which STDOUT will be written, %j inserts jobid       
#SBATCH -e logs/ensemble_slurm_%j.err    # File to which STDERR will be written, %j inserts jobid


#Source clean environment with compatible netcdf and compiler environments and packages like GNU parallel:
source {ASSIM}/environments/cheereio.env #This is specific to the Harvard cluster; rewrite for yours
eval "$(conda shell.bash hook)"

### Run directory
ENSDIR=$(pwd -P)

MY_PATH="$(jq -r ".MY_PATH" {ASSIM}/ens_config.json)"
RUN_NAME="$(jq -r ".RUN_NAME" {ASSIM}/ens_config.json)"

### Get current task ID
x=${SLURM_ARRAY_TASK_ID}

### Add zeros to the current task ID
if [ $x -lt 10 ]; then
  xstr="000${x}"
elif [ $x -lt 100 ]; then
  xstr="00${x}"
elif [ $x -lt 1000 ]; then
  xstr="0${x}"
else
  xstr="${x}"
fi

firstrun=true
scalefirst="$(jq -r ".SIMPLE_SCALE_FOR_FIRST_ASSIM_PERIOD" {ASSIM}/ens_config.json)"
trigger_burnin_scale=false
doburnin="$(jq -r ".DO_BURN_IN" {ASSIM}/ens_config.json)"
scaleburnin="$(jq -r ".SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD" {ASSIM}/ens_config.json)"
amplifyspread="$(jq -r ".AMPLIFY_ENSEMBLE_SPREAD_FOR_FIRST_ASSIM_PERIOD" {ASSIM}/ens_config.json)"

### Run GEOS-Chem in the directory corresponding to the cluster Id
cd  {RunName}_${xstr}

# Set the proper # of threads for OpenMP
# SLURM_CPUS_PER_TASK ensures this matches the number you set with NumCores in the ens_config file
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#Run GC; hang until assimilation complete (later also will do assimilation).
#This will loop until a file appears in scratch signalling assimilation is complete.
while [ ! -f ${MY_PATH}/${RUN_NAME}/scratch/ENSEMBLE_COMPLETE ]; do
# Run GEOS_Chem.  The "time" command will return CPU and wall times.
# Stdout and stderr will be directed to the "GC.log" log file
  srun -c $OMP_NUM_THREADS time -p ./gcclassic >> GC.log
  wait
  taillog="$(tail -n 1 GC.log)"
  #Check if GC finished.
  if [[ "${taillog:0:1}" != "*" ]]; then
    #This file's presence breaks loop loop
    printf "GEOS-Chem in ensemble member ${xstr} did not complete successfully\n" > ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS 
  fi
    #If there is a problem, the KILL_ENS file will be produced. Break then
  if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
    break
  fi
  #Ensemble member 1 handles checking. CD to core.
  if [ $x -eq 1 ]; then
    cd {ASSIM}/core
  fi
  #Hang until ALL_RUNS_COMPLETE found in scratch folder
  until [ -f ${MY_PATH}/${RUN_NAME}/scratch/ALL_RUNS_COMPLETE ]
  do
    #If this is ensemble member 1, look for all restarts and flag if found. Otherwise do nothing.
    if [ $x -eq 1 ]; then
      bash check_for_all_restarts.sh #Check if restarts exist; if they do, create ALL_RUNS_COMPLETE.
    fi
    #If there is a problem, the KILL_ENS file will be produced. Break then
    if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
      break 2
    fi
    sleep 1
  done
  #CD to core
  cd {ASSIM}/core
  #Check if we are in the first assimilation cycle after burn in completes
  if [ -f ${MY_PATH}/${RUN_NAME}/scratch/BURN_IN_PERIOD_PROCESSED ] && [ ! -f ${MY_PATH}/${RUN_NAME}/scratch/BURN_IN_SCALING_COMPLETE ]; then
    trigger_burnin_scale=true
  else
    trigger_burnin_scale=false
  fi
  #Check if we are doing a simple scale or a full assimilation
  if [[ ("${firstrun}" = "true" && "${scalefirst}" = "true") || ("${trigger_burnin_scale}" = "true" && "${doburnin}" = "true" && "${scaleburnin}" = "true") ]]; then
    simplescale=true
  else
    simplescale=false
  fi
  #Check if we are amplifying concentration spreads
  if [[ ("${firstrun}" = "true" && "${amplifyspread}" = "true") ]]; then
    doamplification=true
  else
    doamplification=false
  fi
  #Use GNU parallel to submit parallel sruns, except nature
  if [ $x -ne 0 ]; then
    if [ {MaxPar} -eq 1 ]; then
      bash par_assim.sh ${x} 1 ${firstrun} ${simplescale} ${doamplification}
    else
      parallel -j {MaxPar} "bash par_assim.sh ${x} {1} ${firstrun} ${simplescale} ${doamplification}" ::: {1..{MaxPar}}
    fi 
  fi
  #Hang until assimilation completes or cleanup completes (in case things go too quickly)
  until [ -f ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE ] || [ ! -f ${MY_PATH}/${RUN_NAME}/scratch/ALL_RUNS_COMPLETE ]; do
    #If this is ensemble member 1, check if assimilation is complete; if it is, do the final overwrites.
    if [ $x -eq 1 ]; then
      bash check_and_complete_assimilation.sh ${simplescale}
    fi
    #If there is a problem, the KILL_ENS file will be produced. Break then
    if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
      break 2
    fi
    sleep 1
  done
  #If there is a problem, the KILL_ENS file will be produced. Break then
  if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
    break
  fi
  #If this is ensemble member 1, and this is the first run, switch to main assimilation mode with regular intervals.
  if [ $x -eq 1 ] && [ "${firstrun}" = true ]; then
    if [ "${doburnin}" = true ]; then #If we are burning in, do burn in durfreq update to restarts
      bash change_histrst_durfreq.sh "BURNIN_DURFREQ"
    else 
      bash change_histrst_durfreq.sh "UPDATEDURFREQ"
    fi
  fi
  #If this is ensemble member one, and we have finished burn in, and we are doing burn in, switch to main assimilation mode with regular intervals.
  if [ $x -eq 1 ] && [ "${doburnin}" = true ] && [ "${trigger_burnin_scale}" = true ]; then
    bash change_histrst_durfreq.sh "UPDATEDURFREQ"
  fi
  #If this is ensemble member 1, execute cleanup. This is because we only want it to run once.
  if [ $x -eq 1 ]; then
    if [ "${firstrun}" = true ]; then
      bash cleanup.sh "POSTFIRST" #This also will break us out of this loop when assimilation complete.
    else
      bash cleanup.sh "ASSIM" #This also will break us out of this loop when assimilation complete.
    fi
  fi
  #For all runs, switch off the first run marker.
  if [ "${firstrun}" = true ]; then
    firstrun=false
  fi
  #Hang until cleanup complete, as determined by temp file deletion.
  until [ ! -f ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE ]; do
    #If there is a problem, the KILL_ENS file will be produced. Break then
    if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
      break 2
    fi
    sleep 1
  done

  #CD back to run directory
  cd ${ENSDIR}/{RunName}_${xstr}

  #If there is a problem, the KILL_ENS file will be produced. Break then
  if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
    break
  fi

  #Everything cleaned up; we can head back to the beginning.
done

#If there is a problem, the KILL_ENS file will be produced. Exit with code 1 in that case
if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
  exit 1
else
  exit 0
fi

