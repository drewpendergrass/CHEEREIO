# Path to assimilation setup
ASSIM_PATH="$(jq -r ".ASSIM_PATH" test_config.json)"

thickline="\n===========================================================\n"
thinline="\n-----------------------------------------------------------\n"

# Name for this run
RUN_NAME="$(jq -r ".RUN_NAME" test_config.json)"

# Path where you want to set up assimilation code and run directories
MY_PATH="$(jq -r ".MY_PATH" test_config.json)"
# Path to initial restart file
RESTART_FILE="$(jq -r ".RESTART_FILE" test_config.json)"

START_DATE=$(jq -r ".START_DATE" test_config.json)
ASSIM_DATE=$(jq -r ".ASSIM_DATE" test_config.json)
END_DATE=$(jq -r ".END_DATE" test_config.json)

#Run script settings
NumCores="$(jq -r ".NumCores" test_config.json)"
Partition="$(jq -r ".Partition" test_config.json)"
Memory="$(jq -r ".Memory" test_config.json)"
WallTime="$(jq -r ".WallTime" test_config.json)"
SpinupWallTime="$(jq -r ".SpinupWallTime" test_config.json)"

# Ensemble settings
nEnsemble=$(jq -r ".nEnsemble" test_config.json)
pPERT="$(jq -r ".pPERT" test_config.json)"
SIMULATE_NATURE=$(jq -r ".SIMULATE_NATURE" test_config.json) #if true, we make a run directory to simulate nature. This is for theoretical experimentation and testing

  printf "${thickline}REMOVING EXISTING TESTING DIRECTORIES AND CREATING NEW ONES${thickline}"

rm -rf ${MY_PATH}/${RUN_NAME}
mkdir -p ${MY_PATH}/${RUN_NAME}
cd ${MY_PATH}/${RUN_NAME}

mkdir ensemble_runs
mkdir ensemble_runs/logs
mkdir scratch
echo "GC-CHEERIO uses this directory to save out intermediate data and track its internal state. Modifying contents of this folder can lead to model failure." > scratch/README

cp ${ASSIM_PATH}/templates/run_ensemble_simulations.sh ensemble_runs/
cp ${ASSIM_PATH}/templates/run_ens.sh ensemble_runs/

sed -i -e "s:{RunName}:${RUN_NAME}:g" \
       -e "s:{NumCores}:${NumCores}:g" \
       -e "s:{Partition}:${Partition}:g" \
       -e "s:{Memory}:${Memory}:g" \
       -e "s:{WallTime}:${WallTime}:g" \
       -e "s:{TESTBOOL}:true:g" \
       -e "s:{ASSIM}:${ASSIM_PATH}:g" ensemble_runs/run_ensemble_simulations.sh

if [ SIMULATE_NATURE ]; then
  sed -i -e "s:{START}:0:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ens.sh
else
  sed -i -e "s:{START}:1:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ens.sh
fi

#Horizontal resolution
RES="$(jq -r ".RES" test_config.json)"
grid_res_long=$RES
if [[ $RES = "4.0x5.0" ]]; then
  grid_res='4x5'
  grid_dir=$grid_res
elif [[ $RES = "2.0x2.5" ]]; then
  grid_res='2x25'
  grid_dir='2x2.5'
elif [[ $RES = "0.5x0.625" ]]; then
  grid_res='05x0625'
  grid_dir=$grid_res_long
elif [[ $RES = "0.25x0.3125" ]]; then
  grid_res='025x03125'
  grid_dir=$grid_res_long
else
  grid_res='4x5'
  grid_dir=$grid_res
  printf "\nWarning: Invalid horizontal resolution; defaulting to 4x5.\n"
fi

NEST="$(jq -r ".NEST" test_config.json)"
REGION="$(jq -r ".REGION" test_config.json)"  # NA,AS,CH,EU,Custom

if [[ ${NEST} = "T" ]]; then
  half_polar="F"
  nested_sim="T"
  buffer_zone="3  3  3  3"
  if [[ ${REGION} = "AS" ]] || [[ ${REGION} = "CH" ]]; then
      domain_name="AS"
      if [[ ${grid_res} = "05x0625" ]]; then
        lon_range=" 60.0 150.0"
        lat_range="-11.0  55.0"
      elif [[ ${grid_res} = "025x03125" ]]; then
        lon_range=" 70.0 140.0"
        lat_range=" 15.0  55.0"
      fi
  elif [[ ${REGION} = "EU" ]]; then
    domain_name="EU"
    if [[ ${grid_res} = "05x0625" ]]; then
      lon_range="-30.0 50.0"
      lat_range=" 30.0 70.0"
    elif [[ ${grid_res} = "025x03125" ]]; then
      lon_range="-15.0  40.0"
      lat_range=" 32.75 61.25"
    fi
  elif [[ ${REGION} = "NA" ]]; then
    domain_name="NA"
    if [[ ${grid_res} = "05x0625" ]]; then
      lon_range="-140.0 -40.0"
      lat_range="  10.0  70.0"
    elif [[ ${grid_res} = "025x03125" ]]; then
      lon_range="-130.0  -60.0"
      lat_range="   9.75  60.0"
    fi
  elif [[ ${REGION} = "Custom" ]]; then
    domain_name="custom"
    lon_range="$(jq -r ".LONS" test_config.json)"
    lat_range="$(jq -r ".LATS" test_config.json)"
  else
    printf "\nInvalid horizontal grid domain. Defaulting to global.\n"
    domain_name="global"
    lon_range="-180.0 180.0"
    lat_range=" -90.0  90.0"
    nested_sim="F"
    buffer_zone="0  0  0  0"
  fi
else
  domain_name="global"
  lon_range="-180.0 180.0"
  lat_range=" -90.0  90.0"
  nested_sim="F"
  buffer_zone="0  0  0  0"
  if [[ ${met_name} = "ModelE2.1" ]] || [[ ${met_name} = "ModelE2.2" ]]; then
      if [[ "$grid_res" == "4x5" ]]; then
        half_polar="T"
      else
        half_polar="F"
      fi
  else
    half_polar="T"
  fi
fi
   
    # Initialize (x=0 is nature run (if used), i.e. no perturbation; x=1 is ensemble member 1; etc.)
    if [ $SIMULATE_NATURE ]; then
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
  cp ${RESTART_FILE} ${name}   

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

    cd ${ASSIM_PATH}/testing
    source activate $(jq -r ".CondaEnv" test_config.json) #Activate conda environment.

    #Create initial scaling factors and randomize restarts
    cd ${ASSIM_PATH}/core
    python initialize_scaling_factors.py "TESTING" "${START_DATE}"
    python randomize_restarts.py "${MY_PATH}/${RUN_NAME}/ensemble_runs" "${ASSIM_DATE}_0000"
    python prep_par.py TESTING
    source deactivate #Exit Conda environment

    #Store current time.
    printf "${START_DATE} 000000" > ${MY_PATH}/${RUN_NAME}/scratch/CURRENT_DATE_TIME
    python advance_timestep.py "TESTING" #Update input.geos to first assimilation period.

    ### Navigate back to top-level directory
    cd ${MY_PATH}/${RUN_NAME}

    printf "${thickline}DONE CREATING TESTING DIRECTORIES${thickline}"
