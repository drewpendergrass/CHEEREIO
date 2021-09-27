#!/bin/bash

# This script will set up asynchronous localized ensemble transform kalman filter (4D-LETKF assimilations/inversions) with GEOS-Chem.
# This involves creating an ensemble of run directories configured for the specific needs of the user and the 4D-LETKF procedure.
# Organization and much of the code are heavily inspired or directly lifted from both the createRunDir utility and the excellent 
# CH4 analytical inversion workflow repository; both were initially written by M. Sulprizio.
#
# Usage: bash setup_ensemble.sh. Be sure to customize the switches below and prepare ens_config.json before running.
#
# See README.md for details (dcp, 29/06/2021)

##=======================================================================
## Read user settings. Modify ens_config.json (not this script) if at all possible!
##=======================================================================

# Turn on/off different steps. This will allow you to come back to this
# script and set up different stages later.
SetupTemplateRundir=false
CompileTemplateRundir=false
SetupSpinupRun=false
SetupEnsembleRuns=true

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

# Start and end date for the production simulations
START_DATE=$(jq -r ".START_DATE" ens_config.json)
END_DATE=$(jq -r ".END_DATE" ens_config.json)

# Time between assimilation (in hours)
ASSIM_TIME=$(jq -r ".ASSIM_TIME" ens_config.json)

#Run script settings
NumCores="$(jq -r ".NumCores" ens_config.json)"
Partition="$(jq -r ".Partition" ens_config.json)"
Memory="$(jq -r ".Memory" ens_config.json)"
WallTime="$(jq -r ".WallTime" ens_config.json)"
SpinupWallTime="$(jq -r ".SpinupWallTime" ens_config.json)"
MaxPar="$(jq -r ".MaxPar" ens_config.json)"

printf " \n"
printf "   ___    _  _     ___     ___     ___     ___     ___     ___   \n"
printf "  / __|  | || |   | __|   | __|   | _ \   | __|   |_ _|   / _ \  \n"
printf " | (__   | __ |   | _|    | _|    |   /   | _|     | |   | (_) | \n"
printf "  \___|  |_||_|   |___|   |___|   |_|_\   |___|   |___|   \___/  \n"
printf '_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""| \n'
printf '"`-0-0-`"`-0-0-`"`-0-0-`"`-0-0-`"`-0-0-`"`-0-0-`"`-0-0-`"`-0-0-` \n'
printf " \n"

# Path to find non-emissions input data; will use if no default found
if [[ -f ${HOME}/.geoschem/config ]]; then
    source ${HOME}/.geoschem/config
    if [[ ! -d ${GC_DATA_ROOT} ]]; then
        printf "\nWarning: Default root data directory does not exist!"
        printf "\nReverting to DATA_PATH from ens_config.json.\n"
        DATA_PATH="$(jq -r ".DATA_PATH" ens_config.json)"
        export GC_DATA_ROOT=${DATA_PATH}
    else
        printf "\nOverriding ens_config.json with default root data directory."
        DATA_PATH=${GC_DATA_ROOT}
    fi
else
    printf "\nWarning: Default root data directory does not exist!"
    printf "\nReverting to DATA_PATH from ens_config.json.\n"
    DATA_PATH="$(jq -r ".DATA_PATH" ens_config.json)"
    export GC_DATA_ROOT=${DATA_PATH}
fi

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
pPERT="$(jq -r ".pPERT" ens_config.json)"
SIMULATE_NATURE=$(jq -r ".SIMULATE_NATURE" ens_config.json) #if true, we make a run directory to simulate nature. This is for theoretical experimentation and testing

# Names of emissions species that we are assimilating, loaded from configuration file
EMIS_TO_ASSIM=( $(jq -r ".CONTROL_VECTOR_EMIS[]" ens_config.json) )

GCC_RUN_FILES="${ASSIM_PATH}/GCClassic/src/GEOS-Chem/run/GCClassic" #Called srcrundir in createRunDir.sh
gcdir="${ASSIM_PATH}/GCClassic/src/GEOS-Chem"
wrapperdir="${ASSIM_PATH}/GCClassic"

# Load file with utility functions to setup configuration files
. ${gcdir}/run/shared/setupConfigFiles.sh

# Define separator lines
thickline="\n===========================================================\n"
thinline="\n-----------------------------------------------------------\n"

##=======================================================================
## Set up template run directory
##=======================================================================
if "$SetupTemplateRundir"; then

# Copy run directory files directly from GEOS-Chem repository

printf "${thickline}CHEERIO TEMPLATE RUN DIRECTORY CREATION${thickline}"

mkdir -p ${MY_PATH}/${RUN_NAME}
cd ${MY_PATH}/${RUN_NAME}
mkdir -p ensemble_runs
mkdir -p ensemble_runs/logs
mkdir -p scratch
mkdir -p postprocess
echo "GC-CHEERIO uses this directory to save out intermediate data and track its internal state. Modifying contents of this folder can lead to model failure." > scratch/README

cp ${ASSIM_PATH}/templates/run_ensemble_simulations.sh ensemble_runs/
cp ${ASSIM_PATH}/templates/run_ens.sh ensemble_runs/

sed -i -e "s:{RunName}:${RUN_NAME}:g" \
       -e "s:{NumCores}:${NumCores}:g" \
       -e "s:{Partition}:${Partition}:g" \
       -e "s:{Memory}:${Memory}:g" \
       -e "s:{WallTime}:${WallTime}:g" \
       -e "s:{TESTBOOL}:false:g" \
       -e "s:{MaxPar}:${MaxPar}:g" \
       -e "s:{ASSIM}:${ASSIM_PATH}:g" ensemble_runs/run_ensemble_simulations.sh

if [ "${SIMULATE_NATURE}" = true ]; then
  sed -i -e "s:{START}:0:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ens.sh
else
  sed -i -e "s:{START}:1:g" -e "s:{END}:${nEnsemble}:g" ensemble_runs/run_ens.sh
fi

mkdir -p ${RUN_TEMPLATE}
cp ${gcdir}/run/shared/cleanRunDir.sh ${RUN_TEMPLATE}
cp ${gcdir}/run/shared/download_data.py ${RUN_TEMPLATE}
cp ${GCC_RUN_FILES}/getRunInfo ${RUN_TEMPLATE}
cp ${GCC_RUN_FILES}/archiveRun.sh ${RUN_TEMPLATE}
cp ${GCC_RUN_FILES}/README ${RUN_TEMPLATE}
cp ${GCC_RUN_FILES}/gitignore ${RUN_TEMPLATE}/.gitignore
cp ${GCC_RUN_FILES}/input.geos.templates/input.geos.${sim_name} ${RUN_TEMPLATE}/input.geos
cp ${GCC_RUN_FILES}/HISTORY.rc.templates/HISTORY.rc.${sim_name} ${RUN_TEMPLATE}/HISTORY.rc
cp ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/HEMCO_Config.rc.${sim_name} ${RUN_TEMPLATE}/HEMCO_Config.rc
# Some simulations (like tagO3) do not have a HEMCO_Diagn.rc file,
# so skip copying it unless the file exists (bmy, 12/11/20)
if [[ -f ${GCC_RUN_FILES}/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.${sim_name} ]]; then
    cp ${GCC_RUN_FILES}/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.${sim_name}  ${RUN_TEMPLATE}/HEMCO_Diagn.rc
fi

# Some simulations should default to online natural emissions for dust, seasalt, soil NO and BVOCs
if [[ ${met_name} = "ModelE2.1" ]] || [[ ${met_name} = "ModelE2.2" ]]; then
    if [[ -f ${GCC_RUN_FILES}/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.${sim_name}.onlineE ]]; then
      cp ./HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.${sim_name}.onlineE  ${RUN_TEMPLATE}/HEMCO_Diagn.rc
    fi
fi

if [[ "x${sim_name}" == "xfullchem" || "x${sim_name}" == "xCH4" ]]; then
    cp -r ${gcdir}/run/shared/metrics.py  ${RUN_TEMPLATE}
    chmod 744 ${RUN_TEMPLATE}/metrics.py
fi

cp -r ${GCC_RUN_FILES}/runScriptSamples ${RUN_TEMPLATE}
mkdir ${RUN_TEMPLATE}/OutputDir
# Set permissions
chmod 744 ${RUN_TEMPLATE}/cleanRunDir.sh
chmod 744 ${RUN_TEMPLATE}/archiveRun.sh
chmod 744 ${RUN_TEMPLATE}/runScriptSamples/*

# Copy species database; append APM or TOMAS species if needed
# Also copy APM input files to the run directory
cp -r ${gcdir}/run/shared/species_database.yml   ${RUN_TEMPLATE}
if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    cat ${gcdir}/run/shared/species_database_tomas.yml >> ${RUN_TEMPLATE}/species_database.yml
elif [[ ${sim_extra_option} =~ "APM" ]]; then
    cat ${gcdir}/run/shared/species_database_apm.yml >> ${RUN_TEMPLATE}/species_database.yml
    cp ${gcdir}/run/shared/apm_tmp.dat ${RUN_TEMPLATE}/apm_tmp.dat
    cp ${gcdir}/run/shared/input.apm   ${RUN_TEMPLATE}/input.apm
fi 

# Create symbolic link to code directory
ln -s ${wrapperdir} ${RUN_TEMPLATE}/CodeDir

# Create build directory
mkdir ${RUN_TEMPLATE}/build
printf "To build GEOS-Chem type:\n   cmake ../CodeDir\n   cmake . -DRUNDIR=..\n   make -j\n   make install\n" >> ${RUN_TEMPLATE}/build/README

cd $RUN_TEMPLATE
mkdir -p OutputDir
mkdir -p Restarts

# Specify meteorology
if [[ ${met_name} = "ModelE2.1" ]]; then
    sed_ie "/### BEGIN SECTION SETTINGS/r ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/header.gcap2"                    HEMCO_Config.rc
    sed_ie "/# --- Meteorology fields for FlexGrid ---/r ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/met_fields.gcap2" HEMCO_Config.rc
else
    sed_ie "/### BEGIN SECTION SETTINGS/r ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/header.gmao"                     HEMCO_Config.rc
    sed_ie "/# --- Meteorology fields for FlexGrid ---/r ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/met_fields.gmao"  HEMCO_Config.rc
fi

# Replace token strings in certain files
sed_ie "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   input.geos
sed_ie "s|{MET}|${met_name}|"             input.geos
sed_ie "s|{SIM}|${sim_name}|"             input.geos
sed_ie "s|{RES}|${grid_res_long}|"        input.geos
sed_ie "s|{NLEV}|${grid_lev}|"            input.geos
sed_ie "s|{LON_RANGE}|${lon_range}|"      input.geos
sed_ie "s|{LAT_RANGE}|${lat_range}|"      input.geos
sed_ie "s|{HALF_POLAR}|${half_polar}|"    input.geos
sed_ie "s|{CENTER_180}|${center_180}|"    input.geos
sed_ie "s|{NESTED_SIM}|${nested_sim}|"    input.geos
sed_ie "s|{BUFFER_ZONE}|${buffer_zone}|"  input.geos
sed_ie "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   HEMCO_Config.rc
sed_ie "s|{GRID_DIR}|${grid_dir}|"        HEMCO_Config.rc
sed_ie "s|{MET_DIR}|${met_dir}|"          HEMCO_Config.rc
sed_ie "s|{NATIVE_RES}|${met_native}|"    HEMCO_Config.rc
sed_ie "s|{LATRES}|${met_latres}|"        HEMCO_Config.rc
sed_ie "s|{LONRES}|${met_lonres}|"        HEMCO_Config.rc
sed_ie "s|{DUST_SF}|${offline_dust_sf}|"  HEMCO_Config.rc
sed_ie "s|{DEAD_TF}|${dead_tf}|"          HEMCO_Config.rc
sed_ie "s|{MET_AVAIL}|${met_avail}|"      HEMCO_Config.rc

# Assign appropriate vertical resolution files for a given simulation
if [[ ${met_name} = "ModelE2.1" ]]; then
    sed_ie          's|{Br_GC}|* Br_GC          $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_Br         2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie         's|{BrO_GC}|* BrO_GC         $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_BrO        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie    's|{GLOBAL_ACTA}|* GLOBAL_ACTA    $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_ACTA       2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_Cl}|* GLOBAL_Cl      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_Cl         2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_ClO}|* GLOBAL_ClO     $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_ClO        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_HCl}|* GLOBAL_HCl     $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_HCl        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie   's|{GLOBAL_HCOOH}|* GLOBAL_HCOOH   $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_HCOOH      2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie    's|{GLOBAL_HNO3}|* GLOBAL_HNO3    $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_HNO3       2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_HO2}|* GLOBAL_HO2     $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_HO2        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie    's|{GLOBAL_HOCl}|* GLOBAL_HOCl    $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_HOCl       2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_NIT}|* GLOBAL_NIT     $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_NIT        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_NO}|* GLOBAL_NO      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_NO         2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_NO2}|* GLOBAL_NO2     $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_NO2        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_NO3}|* GLOBAL_NO3     $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_NO3        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_O3}|* GLOBAL_O3      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_O3         2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_OH}|* GLOBAL_OH      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_OH         2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie       's|{AERO_SO4}|* AERO_SO4       $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_SO4        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie       's|{AERO_NH4}|* AERO_NH4       $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_NH4        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie       's|{AERO_NIT}|* AERO_NIT       $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_NIT        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{AERO_BCPI}|* AERO_BCPI      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_BCPI       2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{AERO_BCPO}|* AERO_BCPO      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_BCPO       2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{AERO_OCPI}|* AERO_OCPI      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_OCPI       2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{AERO_OCPO}|* AERO_OCPO      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_OCPO       2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{AERO_DST1}|* AERO_DST1      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.SpeciesConc.2010-2019.$MM.{VERTRES}L.nc4 SpeciesConc_DST1       2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc 
    sed_ie      's|{GLOBAL_OA}|* GLOBAL_OA      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.AerosolMass.2010-2019.$MM.{VERTRES}L.nc4 TotalOA                2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie        's|{PCO_CH4}|* PCO_CH4        $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.ProdLoss.2010-2019.$MM.{VERTRES}L.nc4    ProdCOfromCH4          2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{PCO_NMVOC}|* PCO_NMVOC      $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.ProdLoss.2010-2019.$MM.{VERTRES}L.nc4    ProdCOfromNMVOC        2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie          's|{PH2O2}|* PH2O2          $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.ProdLoss.2010-2019.$MM.{VERTRES}L.nc4    Prod_H2O2              2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie        's|{O3_PROD}|* O3_PROD        $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.ProdLoss.2010-2019.$MM.{VERTRES}L.nc4    Prod_Ox                2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie        's|{O3_LOSS}|* O3_LOSS        $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.ProdLoss.2010-2019.$MM.{VERTRES}L.nc4    Loss_Ox                2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie           's|{JBrO}|* JBrO           $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.JValues.2010-2019.$MM.{VERTRES}L.nc4     Jval_BrO               2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie          's|{JH2O2}|* JH2O2          $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.JValues.2010-2019.$MM.{VERTRES}L.nc4     Jval_H2O2              2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie           's|{JNO2}|* JNO2           $ROOT/GCAP2/OFFLINE_FIELDS/13.0.0/GEOSChem.JValues.2010-2019.$MM.{VERTRES}L.nc4     Jval_NO2               2015/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie       's|{CH4_LOSS}|* CH4_LOSS       $ROOT/GCAP2/GMI/v2015-02/{VERTRES}L/gmi.clim.CH4.geos5.2x25.nc                      loss                   2005/1-12/1/0 C xyz s-1      * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{CO2_COPROD}|* CO2_COPROD     $ROOT/GCAP2/CO2/v2019-02/CHEM/CO2_prod_rates.2x25.{VERTRES}L.nc                     LCO               2004-2009/1-12/1/0 C xyz kgC/m3/s * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{Br_TOMCAT}|* Br_TOMCAT      $ROOT/GCAP2/MERCURY/v2014-09/BrOx/BrOx.GMI.geos5.2x25.{VERTRES}L.nc4                LBRO2N                 1985/1-12/1/0 C xyz pptv     * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{BrO_TOMCAT}|* BrO_TOMCAT     $ROOT/GCAP2/MERCURY/v2014-09/BrOx/BrOx.GMI.geos5.2x25.{VERTRES}L.nc4                LBRO2H                 1985/1-12/1/0 C xyz pptv     * - 1 1|' HEMCO_Config.rc 
    sed_ie      's|{GLOBAL_OC}|1002 GLOBAL_OC   $ROOT/GCAP2/POPs/v2015-08/OCPO.4x5.{VERTRES}L.nc4                                   OCPO              2005-2009/1-12/1/0 C xyz kg       * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_BC}|1002 GLOBAL_BC   $ROOT/GCAP2/POPs/v2015-08/BCPO.4x5.{VERTRES}L.nc4                                   BCPO              2005-2009/1-12/1/0 C xyz kg       * - 1 1|' HEMCO_Config.rc
    sed_ie  's|{TES_CLIM_CCL4}|* TES_CLIM_CCL4  $ROOT/GCAP2/RRTMG/v2018-11/species_clim_profiles.2x25.{VERTRES}L.nc4                CCl4                   2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie 's|{TES_CLIM_CFC11}|* TES_CLIM_CFC11 $ROOT/GCAP2/RRTMG/v2018-11/species_clim_profiles.2x25.{VERTRES}L.nc4                CFC11                  2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie 's|{TES_CLIM_CFC12}|* TES_CLIM_CFC12 $ROOT/GCAP2/RRTMG/v2018-11/species_clim_profiles.2x25.{VERTRES}L.nc4                CFC12                  2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie 's|{TES_CLIM_CFC22}|* TES_CLIM_CFC22 $ROOT/GCAP2/RRTMG/v2018-11/species_clim_profiles.2x25.{VERTRES}L.nc4                CFC22                  2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie   's|{TES_CLIM_CH4}|* TES_CLIM_CH4   $ROOT/GCAP2/RRTMG/v2018-11/species_clim_profiles.2x25.{VERTRES}L.nc4                CH4                    2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie   's|{TES_CLIM_N2O}|* TES_CLIM_N2O   $ROOT/GCAP2/RRTMG/v2018-11/species_clim_profiles.2x25.{VERTRES}L.nc4                N2O                    2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie    's|{GMI_LOSS_CO}|* GMI_LOSS_CO    $ROOT/GCAP2/GMI/v2015-02/{VERTRES}L/gmi.clim.CO.geos5.2x25.nc                       loss                   2005/1-12/1/0 C xyz s-1     CO - 1 1|' HEMCO_Config.rc
    sed_ie    's|{GMI_PROD_CO}|* GMI_PROD_CO    $ROOT/GCAP2/GMI/v2015-02/{VERTRES}L/gmi.clim.CO.geos5.2x25.nc                       prod                   2005/1-12/1/0 C xyz v/v/s   CO - 1 1|' HEMCO_Config.rc    
    sed_ie     's|{OCEAN_MASK}|1000 OCEAN_MASK  $METDIR/TOPO                                                                        focean                    */1/1/0 C xy 1 1  -180/-90/180/90|' HEMCO_Config.rc
    sed_ie        's|{Bry_DIR}|GCAP2/Bry/v2015-01/{VERTRES}L|'                                                                                                                                    HEMCO_Config.rc
    sed_ie        's|{GMI_DIR}|GCAP2/GMI/v2015-02/{VERTRES}L|'                                                                                                                                    HEMCO_Config.rc
    sed_ie        's|{UCX_DIR}|GCAP2/UCX/v2018-02|'                                                                                                                                               HEMCO_Config.rc
    sed_ie         "s|{GFEDON}|off|"                                                                                                                                                              HEMCO_Config.rc
    sed_ie         "s|{ONLINE}|on |"                                                                                                                                                              HEMCO_Config.rc
    sed_ie       "s|{UCX_LEVS}|${grid_lev}|"                                                                                                                                                      HEMCO_Config.rc
    sed_ie       "s|{SCENARIO}|${scenario}|"                                                                                                                                                      HEMCO_Config.rc
    sed_ie          "s|{RUNID}|${runid}|"                                                                                                                                                         HEMCO_Config.rc
    sed_ie      "s|{VOLC_YEAR}|${volc_year}|g"                                                                                                                                                    HEMCO_Config.rc
    sed_ie        "s|{DEAD_TF}|${dead_tf}|"                                                                                                                                                       HEMCO_Config.rc
    sed_ie        "s|{GISSRES}|${giss_res}|g"                                                                                                                                                     HEMCO_Config.rc
    sed_ie        "s|{VERTRES}|${grid_lev}|g"                                                                                                                                                     HEMCO_Config.rc
else
    sed_ie    's|{GLOBAL_ACTA}|* GLOBAL_ACTA    $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_ACTA  2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie          's|{Br_GC}|* Br_GC          $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_Br    2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie         's|{BrO_GC}|* BrO_GC         $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_BrO   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_Cl}|* GLOBAL_Cl      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_Cl    2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_ClO}|* GLOBAL_ClO     $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_ClO   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_HCl}|* GLOBAL_HCl     $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_HCl   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie   's|{GLOBAL_HCOOH}|* GLOBAL_HCOOH   $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_HCOOH 2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie    's|{GLOBAL_HNO3}|* GLOBAL_HNO3    $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_HNO3  2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_HO2}|* GLOBAL_HO2     $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_HO2   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie    's|{GLOBAL_HOCl}|* GLOBAL_HOCl    $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_HOCl  2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_NIT}|* GLOBAL_NIT     $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_NIT   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_NO}|* GLOBAL_NO      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_NO    2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_NO2}|* GLOBAL_NO2     $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_NO2   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{GLOBAL_NO3}|* GLOBAL_NO3     $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_NO3   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_O3}|* GLOBAL_O3      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_O3    2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_OH}|* GLOBAL_OH      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_OH    2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie       's|{AERO_SO4}|* AERO_SO4       $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_SO4   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie       's|{AERO_NH4}|* AERO_NH4       $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_NH4   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie       's|{AERO_NIT}|* AERO_NIT       $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_NIT   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{AERO_BCPI}|* AERO_BCPI      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_BCPI  2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{AERO_BCPO}|* AERO_BCPO      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_BCPO  2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{AERO_OCPI}|* AERO_OCPI      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_OCPI  2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{AERO_OCPO}|* AERO_OCPO      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_OCPO  2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{AERO_DST1}|* AERO_DST1      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.SpeciesConc.$YYYY$MM$DD_0000z.nc4      SpeciesConc_DST1  2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc 
    sed_ie      's|{GLOBAL_OA}|* GLOBAL_OA      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.AerosolMass.$YYYY$MM$DD_0000z.nc4      TotalOA           2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie        's|{PCO_CH4}|* PCO_CH4        $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.ProdLoss.$YYYY$MM$DD_0000z.nc4         ProdCOfromCH4     2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{PCO_NMVOC}|* PCO_NMVOC      $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.ProdLoss.$YYYY$MM$DD_0000z.nc4         ProdCOfromNMVOC   2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie          's|{PH2O2}|* PH2O2          $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.ProdLoss.$YYYY$MM$DD_0000z.nc4         Prod_H2O2         2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie        's|{O3_PROD}|* O3_PROD        $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.ProdLoss.$YYYY$MM$DD_0000z.nc4         Prod_Ox           2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie        's|{O3_LOSS}|* O3_LOSS        $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.ProdLoss.$YYYY$MM$DD_0000z.nc4         Loss_Ox           2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie           's|{JBrO}|* JBrO           $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.JValues.$YYYY$MM$DD_0000z.nc4          Jval_BrO          2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie          's|{JH2O2}|* JH2O2          $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.JValues.$YYYY$MM$DD_0000z.nc4          Jval_H2O2         2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie           's|{JNO2}|* JNO2           $ROOT/GCClassic_Output/13.0.0/$YYYY/GEOSChem.JValues.$YYYY$MM$DD_0000z.nc4          Jval_NO2          2010-2019/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie       's|{CH4_LOSS}|* CH4_LOSS       $ROOT/CH4/v2014-09/4x5/gmi.ch4loss.geos5_47L.4x5.nc                                 CH4loss           1985/1-12/1/0      C xyz s-1      * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{CO2_COPROD}|* CO2_COPROD     $ROOT/CO2/v2019-02/CHEM/CO2_prod_rates.GEOS5.2x25.47L.nc                            LCO               2004-2009/1-12/1/0 C xyz kgC/m3/s * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{Br_TOMCAT}|* Br_TOMCAT      $ROOT/MERCURY/v2014-09/BrOx/BrOx.GMI.geos5.2x25.nc                                  LBRO2N                 1985/1-12/1/0 C xyz pptv     * - 1 1|' HEMCO_Config.rc
    sed_ie     's|{BrO_TOMCAT}|* BrO_TOMCAT     $ROOT/MERCURY/v2014-09/BrOx/BrOx.GMI.geos5.2x25.nc                                  LBRO2H                 1985/1-12/1/0 C xyz pptv     * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_OC}|1002 GLOBAL_OC   $ROOT/POPs/v2015-08/OCPO.$met.$RES.nc                                               OCPO              2005-2009/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie      's|{GLOBAL_BC}|1002 GLOBAL_BC   $ROOT/POPs/v2015-08/BCPO.$met.$RES.nc                                               BCPO              2005-2009/1-12/1/0 C xyz 1        * - 1 1|' HEMCO_Config.rc
    sed_ie  's|{TES_CLIM_CCL4}|* TES_CLIM_CCL4  $ROOT/RRTMG/v2018-11/species_clim_profiles.2x25.nc                                  CCl4                   2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie 's|{TES_CLIM_CFC11}|* TES_CLIM_CFC11 $ROOT/RRTMG/v2018-11/species_clim_profiles.2x25.nc                                  CFC11                  2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie 's|{TES_CLIM_CFC12}|* TES_CLIM_CFC12 $ROOT/RRTMG/v2018-11/species_clim_profiles.2x25.nc                                  CFC12                  2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie 's|{TES_CLIM_CFC22}|* TES_CLIM_CFC22 $ROOT/RRTMG/v2018-11/species_clim_profiles.2x25.nc                                  CFC22                  2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie   's|{TES_CLIM_CH4}|* TES_CLIM_CH4   $ROOT/RRTMG/v2018-11/species_clim_profiles.2x25.nc                                  CH4                    2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie   's|{TES_CLIM_N2O}|* TES_CLIM_N2O   $ROOT/RRTMG/v2018-11/species_clim_profiles.2x25.nc                                  N2O                    2000/1/1/0    C xyz ppbv     * - 1 1|' HEMCO_Config.rc
    sed_ie    's|{GMI_LOSS_CO}|* GMI_LOSS_CO    $ROOT/GMI/v2015-02/gmi.clim.CO.geos5.2x25.nc                                        loss                   2005/1-12/1/0 C xyz s-1     CO - 1 1|' HEMCO_Config.rc
    sed_ie    's|{GMI_PROD_CO}|* GMI_PROD_CO    $ROOT/GMI/v2015-02/gmi.clim.CO.geos5.2x25.nc                                        prod                   2005/1-12/1/0 C xyz v/v/s   CO - 1 1|' HEMCO_Config.rc 
    sed_ie     's|{OCEAN_MASK}|1000 OCEAN_MASK  $METDIR/$CNYR/01/$MET.$CNYR0101.CN.$RES.$NC                                         FROCEAN           2000/1/1/0 C xy 1 1       -180/-90/180/90|' HEMCO_Config.rc
    sed_ie        's|{Bry_DIR}|STRAT/v2015-01/Bry|'                                                                                                                                               HEMCO_Config.rc
    sed_ie        's|{GMI_DIR}|GMI/v2015-02|'                                                                                                                                                     HEMCO_Config.rc
    sed_ie        's|{UCX_DIR}|UCX/v2018-02|'                                                                                                                                                     HEMCO_Config.rc
    sed_ie         "s|{GFEDON}|on |"                                                                                                                                                              HEMCO_Config.rc
    sed_ie         "s|{ONLINE}|off|"                                                                                                                                                              HEMCO_Config.rc
    sed_ie       's|{UCX_LEVS}|72|'                                                                                                                                                               HEMCO_Config.rc
    sed_ie      "s|{VOLC_YEAR}|${volc_year}|g"                                                                                                                                                    HEMCO_Config.rc
    sed_ie        's|{VERTRES}|47|'                                                                                                                                                               HEMCO_Config.rc
fi

printf "\n  -- This template run directory has no start or end date in input.geos."
printf "\n     This will be modified in ensemble run directories.\n"


###Update HISTORY.rc to save a NetCDF every assimilation cycle
###Default frequency is 1 hour, change if you'd like!
dur_days=$((ASSIM_TIME % 24))
dur_hours=$((ASSIM_TIME / 24))
  ### Add zeros to string name
  if [ $dur_days -lt 10 ]; then
      dur_days_str="0${dur_days}"
  else
      dur_days_str="${dur_days}"
  fi
  ### Add zeros to string name
  if [ $dur_hours -lt 10 ]; then
      dur_hours_str="0${dur_hours}"
  else
      dur_hours_str="${dur_hours}"
  fi
sed_ie "s|{FREQUENCY}|00000000 010000|"  HISTORY.rc
sed_ie "s|{DURATION}|000000${dur_days_str} ${dur_hours_str}0000|"   HISTORY.rc

printf "\n  -- The default frequency of diagnostics is set to hourly"
printf "\n     with a duration of once per assimlation cycle."
printf "\n     You may modify these settings in HISTORY.rc and HEMCO_Config.rc.\n"

# Call function to setup configuration files with settings common between
# GEOS-Chem Classic and GCHP.
if [[ "x${sim_name}" == "xfullchem" ]]; then
    set_common_settings ${sim_extra_option}
fi

# Modify input files for benchmark that are specific to GEOS-Chem Classic
if [[ "x${sim_extra_option}" == "xbenchmark" ]]; then
    replace_colon_sep_val "Use GC classic timers?"   T    input.geos
fi

# Modify input files for TOMAS that are specific to GEOS-Chem Classic
if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    replace_colon_sep_val "Tran/conv timestep [sec]" 1800 input.geos
    replace_colon_sep_val "Chem/emis timestep [sec]" 3600 input.geos
fi

# Set online dust emission mass tuning factor according to met field
# and resolution. These values were obtained from hcox_dustdead_mod.F90.
if [[ "x${sim_name}" == "xfullchem" || "x${sim_name}" == "xaerosol" ]]; then
    if [[ "x${met_name}" == "xGEOSFP" && "x${grid_res}" == "x4x5" ]]; then
  replace_colon_sep_val "--> Mass tuning factor" 8.3286e-4 HEMCO_Config.rc
    fi
    if [[ "x${met_name}" == "xGEOSFP" && "x${grid_res}" == "x2x25" ]]; then
  replace_colon_sep_val "--> Mass tuning factor" 5.0416e-4 HEMCO_Config.rc
    fi
    if [[ "x${met_name}" == "xMERRA2" && "x${grid_res}" == "x4x5" ]]; then
  replace_colon_sep_val "--> Mass tuning factor" 7.8533e-4 HEMCO_Config.rc
    fi
    if [[ "x${met_name}" == "xMERRA2" && "x${grid_res}" == "x2x25" ]]; then
  replace_colon_sep_val "--> Mass tuning factor" 4.7586e-4 HEMCO_Config.rc
    fi
fi

# Modify input files for troposphere-only chemistry grids
if [[ "x${chemgrid}" == "xtrop_only" ]]; then
    replace_colon_sep_val "=> Set init. strat. H2O"  F input.geos
    replace_colon_sep_val "Settle strat. aerosols"   F input.geos
    replace_colon_sep_val "Online PSC AEROSOLS"      F input.geos
    replace_colon_sep_val "Perform PSC het. chem.?"  F input.geos
    replace_colon_sep_val "Calc. strat. aero. OD?"   F input.geos
    replace_colon_sep_val "Use UCX strat. chem?"     F input.geos
    replace_colon_sep_val "Active strat. H2O?"       F input.geos
    replace_colon_sep_val "--> STATE_PSC"        false HEMCO_Config.rc
    replace_colon_sep_val "--> GMI_PROD_LOSS"    false HEMCO_Config.rc
    replace_colon_sep_val "--> UCX_PROD_LOSS"     true HEMCO_Config.rc
    sed_ie "s|'Chem_StatePSC|#'Chem_StatePSC|"      HISTORY.rc
fi

# Modify input files for nested-grid simulations
if [[ "x${nested_sim}" == "xT" ]]; then
    replace_colon_sep_val "--> GC_BCs" true HEMCO_Config.rc
    #Point to user-specified boundary conditions
    insert_text "(((GC_BCs" "{BC_fill}" HEMCO_Config.rc
    remove_text "* BC_"                 HEMCO_Config.rc
    sed_ie        "s|{BC_fill}|* BC_  ${BC_FILES} SpeciesBC_?ADV?  1980-2021/1-12/1-31/* EFY xyz 1 * - 1 1|" HEMCO_Config.rc
    if [[ "x${domain_name}" == "xNA" ]]; then
      replace_colon_sep_val "--> NEI2011_MONMEAN" false HEMCO_Config.rc
      replace_colon_sep_val "--> NEI2011_HOURLY"  true  HEMCO_Config.rc
    fi
fi

# Modify input files for POPs simulations
if [[ ${sim_name} =~ "POPs" ]]; then
    sed_ie "s|{POPs_SPC}|${POP_SPC}|"               input.geos
    sed_ie "s|{POPs_XMW}|${POP_XMW}|"               input.geos
    sed_ie "s|{POPs_KOA}|${POP_KOA}|"               input.geos
    sed_ie "s|{POPs_KBC}|${POP_KBC}|"               input.geos
    sed_ie "s|{POPs_K_POPG_OH}|${POP_K_POPG_OH}|"   input.geos
    sed_ie "s|{POPs_K_POPG_O3A}|${POP_K_POPG_O3A}|" input.geos
    sed_ie "s|{POPs_K_POPG_O3B}|${POP_K_POPG_O3B}|" input.geos
    sed_ie "s|{POPs_HSTAR}|${POP_HSTAR}|"           input.geos
    sed_ie "s|{POPs_DEL_H}|${POP_DEL_H}|"           input.geos
    sed_ie "s|{POPs_DEL_Hw}|${POP_DEL_Hw}|"         input.geos
    sed_ie "s|{POPs_SPC}|${POP_SPC}|"               HEMCO_Config.rc
    sed_ie "s|{POPs_SPC}|${POP_SPC}|"               HEMCO_Config.rc
    sed_ie "s|{POPs_SPC}|${POP_SPC}|"               HEMCO_Config.rc
    sed_ie "s|{POPs_SPC}|${POP_SPC}|"               HEMCO_Diagn.rc
fi

#--------------------------------------------------------------------
# Change timesteps for nested-grid simulations
# Transport should be 300s (5min); chemistry should be 600s (10min)
#--------------------------------------------------------------------
if [[ "x${domain_name}" == "xAS"     ]] || \
       [[ "x${domain_name}" == "xEU"     ]] || \
       [[ "x${domain_name}" == "xNA"     ]] || \
       [[ "x${domain_name}" == "xcustom" ]]; then
    cmd='s|\[sec\]: 600|\[sec\]: 300|'
    sed_ie "$cmd" input.geos
    cmd='s|\[sec\]: 1200|\[sec\]: 600|'
    sed_ie "$cmd" input.geos
fi

# Modify default settings for GCAP simulations
if [[ ${met_name} = "ModelE2.1" ]]; then
    replace_colon_sep_val "--> CEDS"                  false HEMCO_Config.rc
    replace_colon_sep_val "--> POET_EOH"              false HEMCO_Config.rc
    replace_colon_sep_val "--> TZOMPASOSA_C2H6"       false HEMCO_Config.rc
    replace_colon_sep_val "--> XIAO_C3H8"             false HEMCO_Config.rc
    replace_colon_sep_val "--> AEIC"                  false HEMCO_Config.rc
    replace_colon_sep_val "--> CEDS_SHIP"             false HEMCO_Config.rc
    replace_colon_sep_val "--> OFFLINE_DUST"          false HEMCO_Config.rc
    replace_colon_sep_val "--> OFFLINE_BIOGENICVOC"   false HEMCO_Config.rc
    replace_colon_sep_val "--> OFFLINE_SEASALT"       false HEMCO_Config.rc
    replace_colon_sep_val "--> OFFLINE_SOILNOX"       false HEMCO_Config.rc
    replace_colon_sep_val "--> GMD_SFC_CH4"           false HEMCO_Config.rc
    replace_colon_sep_val "--> QFED2"                 false HEMCO_Config.rc
    replace_colon_sep_val "--> SfcVMR"                false HEMCO_Config.rc
    replace_colon_sep_val "--> CMIP6_SFC_BC"          true  HEMCO_Config.rc
    replace_colon_sep_val "--> CMIP6_SFC_LAND_ANTHRO" true  HEMCO_Config.rc
    replace_colon_sep_val "--> CMIP6_AIRCRAFT"        true  HEMCO_Config.rc
    replace_colon_sep_val "--> CMIP6_SHIP"            true  HEMCO_Config.rc
    replace_colon_sep_val "--> BB4MIPS"               true  HEMCO_Config.rc
fi

#Modify HEMCO_Config so that GEOS-Chem will read in assimilated scaling factors.
cd ${ASSIM_PATH}/core
bash prepare_template_hemco_config.sh
printf "${thinline}"


printf "\n  -- Template run directory created at ${MY_PATH}/${RUN_NAME}/${RUN_TEMPLATE}."
printf "\n     Modify all configuration files BEFORE creating spinup or ensemble run directories."
printf "\n     All settings from this folder will be repeated throughout the ensemble.\n"

printf "${thinline}CHEERIO TEMPLATE RUN DIRECTORY CREATED${thinline}"

### Navigate back to top-level directory
cd ${MY_PATH}/${RUN_NAME}

fi # SetupTemplateRunDir

##=======================================================================
##  Compile template run directory
##=======================================================================

### Compile GEOS-Chem and store executable in template run directory
if "$CompileTemplateRundir"; then

  printf "${thickline}COMPILING TEMPLATE RUN DIRECTORY${thickline}"

  cd ${MY_PATH}/${RUN_NAME}/${RUN_TEMPLATE}/build
  cmake ../CodeDir
  cmake . -DRUNDIR=..
  make -j
  make install

  printf "${thinline}TEMPLATE RUN DIRECTORY COMPILED${thinline}"

  ### Navigate back to top-level directory
  cd ${MY_PATH}/${RUN_NAME}

fi

##=======================================================================
##  Set up spinup run directory
##=======================================================================
if  "$SetupSpinupRun"; then

    printf "${thickline}CHEERIO SPINUP RUN DIRECTORY CREATION${thickline}"

    cd ${MY_PATH}/${RUN_NAME}
    
    ### Define the run directory name
    spinup_name="${RUN_NAME}_Spinup"

    ### Make the directory
    runDir="spinup_run"
    mkdir -p ${runDir}

    ### Copy and point to the necessary data
    cp -r ${RUN_TEMPLATE}/*  ${runDir}
    cp -RLv ${ASSIM_PATH}/templates/ensemble_run.template ${runDir}
    cd $runDir

    ### Link to GEOS-Chem executable instead of having a copy in each rundir
    rm -rf gcclassic
    ln -s ../${RUN_TEMPLATE}/gcclassic .
    # Link to restart file
    ln -s $RESTART_FILE GEOSChem.Restart.${SPINUP_START}_0000z.nc4

    #Switch HEMCO_Config to base/nature one.
    rm HEMCO_Config.rc #This one has updated scaling factors.
    mv HEMCO_Config_SPINUP_NATURE_TEMPLATE.rc HEMCO_Config.rc #This one only updates BCs.
    
    ### Update settings in input.geos
    sed -i -e "s:{DATE1}:${SPINUP_START}:g" \
           -e "s:{DATE2}:${SPINUP_END}:g" \
           -e "s:{TIME1}:000000:g" \
           -e "s:{TIME2}:000000:g" input.geos

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
    
fi # SetupSpinupRun

##=======================================================================
##  Set up Ensemble run directories
##=======================================================================
if "$SetupEnsembleRuns"; then

    printf "${thickline}CHEERIO ENSEMBLE RUN DIRECTORY CREATION${thickline}"
    
    # Initialize (x=0 is nature run (if used), i.e. no perturbation; x=1 is ensemble member 1; etc.)
    if [ "${SIMULATE_NATURE}" = true ]; then
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
    rm HEMCO_Config.rc #This one has updated scaling factors.
    mv HEMCO_Config_SPINUP_NATURE_TEMPLATE.rc HEMCO_Config.rc #This one only updates BCs.
  else 
    #Use HEMCO_Config with updated scaling factors
    rm HEMCO_Config_SPINUP_NATURE_TEMPLATE.rc
    sed_ie "s|template_run|ensemble_runs/${name}|"  HEMCO_Config.rc #Replace template_run with this folder in HEMCO_Config
  fi

  # Link to restart file
  if "$DO_SPINUP"; then
      ln -s ../../spinup_run/GEOSChem.Restart.${SPINUP_END}_0000z.nc4 GEOSChem.Restart.${START_DATE}_0000z.nc4
  else
      ln -s $RESTART_FILE GEOSChem.Restart.${START_DATE}_0000z.nc4
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
    python initialize_scaling_factors.py "PRODUCTION" "${START_DATE}"
    python prep_par.py PRODUCTION
    conda deactivate #Exit Conda environment

    #Store current time.
    printf "${START_DATE} 000000" > ${MY_PATH}/${RUN_NAME}/scratch/CURRENT_DATE_TIME
    bash update_input_geos.sh #Update input.geos to first assimilation period.

    ### Navigate back to top-level directory
    cd ${MY_PATH}/${RUN_NAME}

    printf "${thickline}DONE CREATING ENSEMBLE MEMBER RUN DIRECTORIES${thickline}"


fi  # SetupEnsembleRuns

exit 0
