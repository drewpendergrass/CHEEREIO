#!/bin/bash

# This script sets up the variables that will be referenced by all setup_ensemble scripts. Users never run this directly.

eval "$(conda shell.bash hook)"

python core/validate_ensconfig.py #Do some checks to see if ens_config looks good.
ret=$?
if [ $ret -ne 0 ]; then
     #Handle failure
     exit 1
fi

GC_VERSION="$(jq -r ".GC_VERSION" ens_config.json)"

# Path to assimilation setup
ASSIM_PATH="$(jq -r ".ASSIM_PATH" ens_config.json)"

# Name for this run
RUN_NAME="$(jq -r ".RUN_NAME" ens_config.json)"

# Path where you want to set up assimilation code and run directories
MY_PATH="$(jq -r ".MY_PATH" ens_config.json)"

# Start and end date for the spinup simulation
DO_SPINUP=$(jq -r ".DO_SPINUP" ens_config.json)
SPINUP_START=$(jq -r ".SPINUP_START" ens_config.json)
SPINUP_END=$(jq -r ".SPINUP_END" ens_config.json)
DO_ENS_SPINUP=$(jq -r ".DO_ENS_SPINUP" ens_config.json)
ENS_SPINUP_FROM_BC_RESTART=$(jq -r ".ENS_SPINUP_FROM_BC_RESTART" ens_config.json)
ENS_SPINUP_START=$(jq -r ".ENS_SPINUP_START" ens_config.json)
ENS_SPINUP_END=$(jq -r ".ENS_SPINUP_END" ens_config.json)

# Start and end date for the control (no assimilation) simulations
CONTROL_START=$(jq -r ".CONTROL_START" ens_config.json)
CONTROL_END=$(jq -r ".CONTROL_END" ens_config.json)


# Start and end date for the production simulations
START_DATE=$(jq -r ".START_DATE" ens_config.json)
END_DATE=$(jq -r ".END_DATE" ens_config.json)

# Time between assimilation (in hours)
ASSIM_TIME=$(jq -r ".ASSIM_TIME" ens_config.json)

#Run script settings
NumCores="$(jq -r ".NumCores" ens_config.json)"
NumCtrlCores="$(jq -r ".NumCtrlCores" ens_config.json)"
Partition="$(jq -r ".Partition" ens_config.json)"
Memory="$(jq -r ".Memory" ens_config.json)"
WallTime="$(jq -r ".WallTime" ens_config.json)"
SpinupWallTime="$(jq -r ".SpinupWallTime" ens_config.json)"
ControlWallTime="$(jq -r ".ControlWallTime" ens_config.json)"
EnsCtrlSpinupMemory="$(jq -r ".EnsCtrlSpinupMemory" ens_config.json)"
EnsSpinupWallTime="$(jq -r ".EnsSpinupWallTime" ens_config.json)"
MaxPar="$(jq -r ".MaxPar" ens_config.json)"

SaveDOFS=$(jq -r ".SaveDOFS" ens_config.json) 


# Path to find non-emissions input data; will use if no default found
if [[ -f ${HOME}/.geoschem/config ]]; then
    source ${HOME}/.geoschem/config
    if [[ ! -d ${GC_DATA_ROOT} ]]; then
        printf "\nWarning: Default root data directory does not exist!"
        printf "\nReverting to DATA_PATH from ens_config.json.\n"
        DATA_PATH="$(jq -r ".DATA_PATH" ens_config.json)"
        export GC_DATA_ROOT=${DATA_PATH}
    else
        DATA_PATH=${GC_DATA_ROOT}
    fi
else
    printf "\nWarning: Default root data directory does not exist!"
    printf "\nReverting to DATA_PATH from ens_config.json.\n"
    DATA_PATH="$(jq -r ".DATA_PATH" ens_config.json)"
    export GC_DATA_ROOT=${DATA_PATH}
fi

CH4_HEMCO_ROOT="$(jq -r ".CH4_HEMCO_ROOT" ens_config.json)"
USE_CUSTOM_CH4="$(jq -r ".USE_CHEEREIO_TEMPLATE_CH4_HEMCO_Config" ens_config.json)"

# Path to initial restart file
RESTART_FILE="$(jq -r ".RESTART_FILE" ens_config.json)"

# Path to boundary condition files (for nested grid simulations)
# Must put backslash before $ in $YYYY$MM$DD to properly work in sed command
BC_FILES="$(jq -r ".BC_FILES" ens_config.json)"

RUN_TEMPLATE="template_run"

valid_sims=("fullchem" "aerosol" "CH4" "CO2" "Hg" "POPs" "tagCH4" "tagCO" "tagO3" "TransportTracers")
sim_name="$(jq -r ".sim_name" ens_config.json)"

if [[ ! " ${valid_sims[@]} " =~ " ${sim_name} " ]]; then
    sim_name="fullchem"
    printf "\nWarning: Invalid sim_name from ens_config.json; defaulting to fullchem.\n"
fi

#Simulation options
sim_extra_option="$(jq -r ".sim_extra_option" ens_config.json)"

if [[ ${sim_name} = "fullchem" ]]; then
  valid_chemgrids=("trop+strat" "trop_only")
  valid_extraoptions=("none" "benchmark" "complexSOA" "complexSOA_SVPOA" "marinePOA" "aciduptake" "TOMAS15" "TOMAS40" "APM" "RRTMG")
  chemgrid="$(jq -r ".chemgrid" ens_config.json)"
  if [[ ! " ${valid_chemgrids[@]} " =~ " ${chemgrid} " ]]; then
    chemgrid="trop+strat"
    printf "\nWarning: Invalid chemgrid from ens_config.json; defaulting to trop+strat.\n"
  fi
  if [[ ! " ${valid_extraoptions[@]} " =~ " ${sim_extra_option} " ]]; then
    sim_extra_option="none"
    printf "\nWarning: Invalid sim_extra_option from ens_config.json; defaulting to none.\n"
  fi

elif [[ ${sim_name} = "TransportTracers" ]]; then
  printf "\nTransportTracers sim_extra_option set to none.\n"
  sim_extra_option="none"

#Pops species specification
elif [[ ${sim_name} = "POPs" ]]; then
  if [[ ${sim_extra_option} = "BaP" ]]; then
      POP_SPC="BaP"
      POP_XMW="252.31d-3"
      POP_KOA="3.02d11"
      POP_KBC="7.94d13"
      POP_K_POPG_OH="50d-12"
      POP_K_POPG_O3A="0d0"
      POP_K_POPG_O3B="2.8d15"
      POP_HSTAR="3.10d-5"
      POP_DEL_H="-110d3"
      POP_DEL_Hw="43d0"
  elif [[ ${sim_extra_option} = "PHE" ]]; then
      POP_SPC="PHE"
      POP_XMW="178.23d-3"
      POP_KOA="4.37d7"
      POP_KBC="1.0d10"
      POP_K_POPG_OH="2.7d-11"
      POP_K_POPG_O3A="0d0"
      POP_K_POPG_O3B="2.15d15"
      POP_HSTAR="1.74d-3"
      POP_DEL_H="-74d3"
      POP_DEL_Hw="47d0"
  elif [[ ${sim_extra_option} = "PYR" ]]; then
      POP_SPC="PYR"
      POP_XMW="202.25d-3"
      POP_KOA="7.24d8"
      POP_KBC="1.0d11"
      POP_K_POPG_OH="50d-12"
      POP_K_POPG_O3A="0d0"
      POP_K_POPG_O3B="3.0d15"
      POP_HSTAR="5.37d-4"
      POP_DEL_H="-87d3"
      POP_DEL_Hw="43d0"
  else
      printf "\nWarning: invalid sim_extra_option from ens_config.json; defaulting to BaP.\n"
      POP_SPC="BaP"
      POP_XMW="252.31d-3"
      POP_KOA="3.02d11"
      POP_KBC="7.94d13"
      POP_K_POPG_OH="50d-12"
      POP_K_POPG_O3A="0d0"
      POP_K_POPG_O3B="2.8d15"
      POP_HSTAR="3.10d-5"
      POP_DEL_H="-110d3"
      POP_DEL_Hw="43d0"
  fi
fi

#Switch to true none type if option is none
if [[ sim_extra_option = "none" ]]; then
  sim_extra_option=none
fi

#Read in meteorology
met_name="$(jq -r ".met_name" ens_config.json)"
if [[ ${met_name} = "MERRA2" ]]; then
  met_name_lc="merra2"
  met_dir='MERRA2'
  met_resolution='05x0625'
  met_native='0.5x0.625'
  met_latres='05'
  met_lonres='0625'
  met_extension='nc4'
  met_cn_year='2015'
  pressure_unit='Pa '
  pressure_scale='0.01'
  offline_dust_sf='3.86e-4'
elif [[ ${met_name} = "GEOSFP" ]]; then
  met_name_lc="geosfp"
  met_dir='GEOS_FP'
  met_resolution='025x03125'
  met_native='0.25x0.3125'
  met_latres='025'
  met_lonres='03125'
  met_extension='nc'
  met_cn_year='2011'
  pressure_unit='hPa'
  pressure_scale='1.0 '
  offline_dust_sf='6.42e-5'
elif [[ ${met_name} = "ModelE2.1" ]]; then
  met_name_lc='modele2.1'
  met_dir='E21'
  met_resolution='2x2.5'
  met_native='2x2.5'
  met_latres='20'
  met_lonres='25'
  met_extension='nc4'
  met_cn_year='1950'
  pressure_unit='Pa '
  pressure_scale='0.01'
  offline_dust_sf='1.0'
else
  met_name="MERRA2"
  met_name_lc="merra2"
  met_dir='MERRA2'
  met_resolution='05x0625'
  met_native='0.5x0.625'
  met_latres='05'
  met_lonres='0625'
  met_extension='nc4'
  met_cn_year='2015'
  pressure_unit='Pa '
  pressure_scale='0.01'
  offline_dust_sf='3.86e-4'
  printf "\nWarning: Invalid meteorology setting; defaulting to MERRA2.\n"
fi

if [[ ${met_name} = "ModelE2.1" ]]; then
  #Note: no checks applied!
  scenario="$(jq -r ".scenario" ens_config.json)"
  runid="$(jq -r ".runid" ens_config.json)"
  met_avail="$(jq -r ".met_avail" ens_config.json)"
  volc_year="$(jq -r ".volc_year" ens_config.json)"
  printf "/nWarning: Did not validate ModelE2.1 settings.\n"
else
    scenario=""
    runid=""
    volc_year='$YYYY'
    met_avail="# 1980-2021"
fi

#Horizontal resolution
RES="$(jq -r ".RES" ens_config.json)"
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

NEST="$(jq -r ".NEST" ens_config.json)"
REGION="$(jq -r ".REGION" ens_config.json)"  # NA,AS,CH,EU,Custom

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
    lon_range="$(jq -r ".LONS" ens_config.json)"
    lat_range="$(jq -r ".LATS" ens_config.json)"
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

#-----------------------------------------------------------------
# Is Int'l Date Line an edge or midpoint?
#-----------------------------------------------------------------
center_180="T" # All GEOS products
if [[ ${met_name} = "ModelE2.1" ]] || [[ ${met_name} = "ModelE2.2" ]] ; then
    case "$grid_res" in
  "2x25" ) center_180="F" ;; # Native GISS fine resolution
        * ) center_180="T" ;; # Flex-grid re-gridded resolutions
    esac  
fi

#-----------------------------------------------------------------
# Horizontal resolution-dependent settings
#-----------------------------------------------------------------
dead_tf="-999.0e0" # Default DEAD-dust scaling factor for online emissions
if [[ ${met_name} = "ModelE2.1" ]]; then
    if [[ "$runid" == "E213f10aF40oQ40nudge" ]]; then
        case "$grid_res" in 
            "4x5" ) dead_tf="0.00474046"; giss_res="F40"  ;;
            "2x25" ) dead_tf="0.00243979"; giss_res="F40"  ;;
            "05x0625" ) dead_tf="0.00276896"; giss_res="F40"  ;;
            "025x03125" ) dead_tf="0.00254319"; giss_res="F40"  ;;
    esac
    else
        case "$grid_res" in 
            "4x5" ) dead_tf="0.03564873"; giss_res="F40"  ;;
            "2x25" ) dead_tf="0.01050036"; giss_res="F40"  ;;
            "05x0625" ) dead_tf="0.01340854"; giss_res="F40"  ;;
            "025x03125" ) dead_tf="0.01066495"; giss_res="F40"  ;;
  esac
    fi
fi

#Number of levels
grid_lev="$(jq -r ".LEVS" ens_config.json)"

# Ensemble settings
nEnsemble=$(jq -r ".nEnsemble" ens_config.json)

dcr=$(jq -r ".DO_CONTROL_RUN" ens_config.json)
dcer=$(jq -r ".DO_CONTROL_WITHIN_ENSEMBLE_RUNS" ens_config.json) #if true, we make a run directory without assimilation within the ensemble runs structure.

if [[ ("${dcr}" = "true" && "${dcer}" = "true") ]]; then
  DO_CONTROL_WITHIN_ENSEMBLE_RUNS=true
else
  DO_CONTROL_WITHIN_ENSEMBLE_RUNS=false
fi


# Names of emissions species that we are assimilating, loaded from configuration file
EMIS_TO_ASSIM=( $(jq -r ".CONTROL_VECTOR_EMIS[]" ens_config.json) )

GCC_RUN_FILES="${ASSIM_PATH}/GCClassic/src/GEOS-Chem/run/GCClassic" #Called srcrundir in createRunDir.sh
srcrundir=$GCC_RUN_FILES
gcdir="${ASSIM_PATH}/GCClassic/src/GEOS-Chem"
wrapperdir="${ASSIM_PATH}/GCClassic"

# Load file with utility functions to setup configuration files
. ${gcdir}/run/shared/setupConfigFiles.sh

# Define separator lines
thickline="\n===========================================================\n"
thinline="\n-----------------------------------------------------------\n"
