#!/bin/bash

# This creates the template run directory for version 13.

##=======================================================================
## Read user settings. Modify ens_config.json (not this script) if at all possible!
##=======================================================================

source setup_ensemble_prelims.sh #get the relevant variables.

# Copy run directory files directly from GEOS-Chem repository

printf "${thickline}CHEERIO TEMPLATE RUN DIRECTORY CREATION${thickline}"

source setup_ensemble_template_shared.sh

# Initialize run directory variables
RUNDIR_VARS=""
RUNDIR_VARS+="RUNDIR_GC_MODE='GCClassic'\n"
RUNDIR_VARS+="RUNDIR_DATA_ROOT=$GC_DATA_ROOT\n"
RUNDIR_VARS+="RUNDIR_SIM_NAME=$sim_name\n"
RUNDIR_VARS+="RUNDIR_SIM_EXTRA_OPTION=$sim_extra_option\n"

# Determine settings based on simulation type
SettingsDir="${gcdir}/run/shared/settings"
if [[ ${sim_extra_option} == "BaP" ]]; then
    RUNDIR_VARS+="$(cat ${SettingsDir}/POPs_BaP.txt)\n"
elif [[ ${sim_extra_option} == "PHE" ]]; then
    RUNDIR_VARS+="$(cat ${SettingsDir}/POPs_PHE.txt)\n"
elif [[ ${sim_extra_option} == "PYR" ]]; then
    RUNDIR_VARS+="$(cat ${SettingsDir}/POPs_PYR.txt)\n"
fi

if [[ ${sim_extra_option} == "benchmark"  ]] || \
   [[ ${sim_extra_option} =~ "complexSOA" ]] || \
   [[ ${sim_extra_option} == "APM"        ]]; then
    RUNDIR_VARS+="RUNDIR_COMPLEX_SOA='true '\n"
    if [[ ${sim_extra_option} == "complexSOA_SVPOA" ]]; then
  RUNDIR_VARS+="RUNDIR_SVPOA='true '\n"
    else
  RUNDIR_VARS+="RUNDIR_SVPOA='false'\n"
    fi
else
    RUNDIR_VARS+="RUNDIR_COMPLEX_SOA='false'\n"
    RUNDIR_VARS+="RUNDIR_SVPOA='false'\n"
fi

if [[ ${sim_extra_option} == "aciduptake" ]]; then
    RUNDIR_VARS+="RUNDIR_DUSTALK_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_ACID_UPTAKE='true '\n"
else
    RUNDIR_VARS+="RUNDIR_DUSTALK_EXT='off'\n"
    RUNDIR_VARS+="RUNDIR_ACID_UPTAKE='false'\n"
fi

if [[ ${sim_extra_option} == "marinePOA" ]]; then
    RUNDIR_VARS+="RUNDIR_MARINE_POA='true '\n"
else
    RUNDIR_VARS+="RUNDIR_MARINE_POA='false'\n"
fi

if [[ ${sim_extra_option} == "RRTMG" ]]; then
    RUNDIR_VARS+="RUNDIR_RRTMG_OPTS='true '\n"
    RUNDIR_VARS+="RUNDIR_USE_RRTMG='true '\n"
else
    RUNDIR_VARS+="RUNDIR_RRTMG_OPTS='false'\n"
    RUNDIR_VARS+="RUNDIR_USE_RRTMG='false'\n"
fi

if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    RUNDIR_VARS+="RUNDIR_USE_NLPBL='false'\n"
    RUNDIR_VARS+="RUNDIR_USE_ONLINE_O3='false'\n"
else
    RUNDIR_VARS+="RUNDIR_USE_NLPBL='true '\n"
    RUNDIR_VARS+="RUNDIR_USE_ONLINE_O3='true '\n"
fi

if [[ ${met_name_lc} = "merra2" ]]; then
  met="merra2"
        shared_met_settings=${gcdir}/run/shared/settings/merra2.txt
  RUNDIR_VARS+="RUNDIR_MET_FIELD_CONFIG='HEMCO_Config.rc.gmao_metfields'\n"
elif [[ ${met_name_lc} = "geosfp" ]]; then
    met="geosfp"
        shared_met_settings=${gcdir}/run/shared/settings/geosfp.txt
  RUNDIR_VARS+="RUNDIR_MET_FIELD_CONFIG='HEMCO_Config.rc.gmao_metfields'\n"
elif [[ ${met_name} = "ModelE2.1" ]]; then
  met="ModelE2.1"
        shared_met_settings=${gcdir}/run/shared/settings/modele2.1.txt
  RUNDIR_VARS+="RUNDIR_MET_FIELD_CONFIG='HEMCO_Config.rc.gcap2_metfields'\n"
fi

#If you're using ModelE, you'll need to modify this.

RUNDIR_VARS+="RUNDIR_GCAP2_SCENARIO='not_used'\n"
RUNDIR_VARS+="RUNDIR_GISS_RES='not_used'\n"
RUNDIR_VARS+="RUNDIR_GCAP2_VERTRES='not_used'\n"
RUNDIR_VARS+="RUNDIR_GCAP2_RUNID='not_used'\n"
RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 1980-2021'\n"

RUNDIR_VARS+="RUNDIR_USE_AEIC='true'\n"

    # Define the volcano paths for the HEMCO_Config.rc file
    # NOTE: Benchmark simulations always use the climatological emissions!
    if [[ "x${sim_name}" == "xfullchem" ]]  ||  \
       [[ "x${sim_name}" == "xaerosol"  ]]; then
          RUNDIR_VARS+="RUNDIR_VOLC_CLIMATOLOGY='\$ROOT/VOLCANO/v2021-09/so2_volcanic_emissions_CARN_v202005.degassing_only.rc'\n"

          if [[ "x${sim_extra_option}" == "xbenchmark" ]]; then
              RUNDIR_VARS+="RUNDIR_VOLC_TABLE='\$ROOT/VOLCANO/v2021-09/so2_volcanic_emissions_CARN_v202005.degassing_only.rc'\n"
          else
              RUNDIR_VARS+="RUNDIR_VOLC_TABLE='\$ROOT/VOLCANO/v2021-09/\$YYYY/\$MM/so2_volcanic_emissions_Carns.\$YYYY\$MM\$DD.rc'\n"
          fi
    fi

if [[ ${grid_res} = "4x5" ]]; then
  RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/4x5.txt)\n"
elif [[ ${grid_res} = "2x25" ]]; then
  RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/2x25.txt)\n"
elif [[ ${grid_res} = "05x0625" ]]; then
  RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/05x0625.txt)\n"
elif [[ ${res_num} = "025x03125" ]]; then
  RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/025x03125.txt)\n"
fi

if [[ ${grid_res} = "05x0625" ]] || [[ ${grid_res} = "025x03125" ]]; then
  if [[ ${NEST} = "T" ]]; then
    RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/nested_grid.txt)\n"
    if [[ ${domain_name} = "AS" ]]; then
      RUNDIR_VARS+="RUNDIR_GRID_DOMAIN_NAME='AS'\n"
      grid_nest="AS"
        if [[ ${grid_res} = "05x0625" ]]; then
            RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[ 60.0, 150.0]'\n"
            RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[-11.0,  55.0]'\n"
        elif [[ ${grid_res} = "025x03125" ]]; then
            RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[ 70.0, 140.0]'\n"
            RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[ 15.0,  55.0]'\n"
        fi
    elif [[ ${domain_name} = "EU" ]]; then
            RUNDIR_VARS+="RUNDIR_GRID_DOMAIN_NAME='EU'\n"
            grid_nest="EU"
            if [[ ${grid_res} = "05x0625" ]]; then
                RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[-30.0, 50.0]'\n"
                RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[ 30.0, 70.0]'\n"
            elif [[ ${grid_res} = "025x03125" ]]; then
                      RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[-15.0,  40.0 ]'\n"
                RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[ 32.75, 61.25]'\n"
            fi
    elif [[ ${domain_name} = "NA" ]]; then
        RUNDIR_VARS+="RUNDIR_GRID_DOMAIN_NAME='NA'\n"
        grid_nest+="NA"
              if [[ ${grid_res} = "05x0625" ]]; then
                  RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[-140.0, -40.0]'\n"
                  RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[  10.0,  70.0]'\n"
              elif [[ ${grid_res} = "025x03125" ]]; then
                  RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[-130.0,  -60.0]'\n"
                  RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[   9.75,  60.0]'\n"
              fi
    elif [[ ${domain_name} = "custom" ]]; then
            grid_nest="CU"
            RUNDIR_VARS+="RUNDIR_GRID_DOMAIN_NAME='custom'\n"
            RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='MinLon MaxLon'\n"
            RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='MinLat MaxLat'\n"
            printf "\n  -- You will need to manually set longitude and latitude"
            printf "\n     bounds in the Grid Menu of geoschem_config.yml!\n"
    else
        printf "Invalid horizontal grid domain option.\n"
    fi
  else
    RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/global_grid.txt)\n"
    if [[ ${met} = "ModelE2.1" ]] || [[ ${met} = "ModelE2.2" ]]; then
        if [[ "$grid_res" == "4x5" ]]; then
            RUNDIR_VARS+="RUNDIR_GRID_HALF_POLAR='true '\n"
        else
            RUNDIR_VARS+="RUNDIR_GRID_HALF_POLAR='false'\n"
        fi
    else
      RUNDIR_VARS+="RUNDIR_GRID_HALF_POLAR='true '\n"
    fi
  fi
fi


RUNDIR_VARS+="$(cat ${shared_met_settings})\n"   # shared_met_settings needs to be included after RUNDIR_GRID_DIR is defined

# Set timesteps according to grid resolution
if [[ ${grid_res} = "05x0625" ]] || [[ ${grid_res} = "025x03125" ]]; then
    RUNDIR_VARS+="RUNDIR_TRANSPORT_TS='300'\n"
    RUNDIR_VARS+="RUNDIR_CHEMISTRY_TS='600'\n"
else
    if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
  RUNDIR_VARS+="RUNDIR_TRANSPORT_TS='1800'\n"
  RUNDIR_VARS+="RUNDIR_CHEMISTRY_TS='3600'\n"
    else
  RUNDIR_VARS+="RUNDIR_TRANSPORT_TS='600'\n"
  RUNDIR_VARS+="RUNDIR_CHEMISTRY_TS='1200'\n"
    fi
fi


#-----------------------------------------------------------------
# Is International Date Line an edge or midpoint?
#-----------------------------------------------------------------

if [[ ${met} = "ModelE2.1" ]] || [[ ${met} = "ModelE2.2" ]] ; then
    if [[ "$grid_res" == "2x25" ]]; then
  # Native GISS fine resolution
  RUNDIR_VARS+="RUNDIR_CENTER_LON_180='false'\n"
    else
        # FlexGrid re-gridded resolutions
  RUNDIR_VARS+="RUNDIR_CENTER_LON_180='true '\n"
    fi
else
    # All GMAO products
    RUNDIR_VARS+="RUNDIR_CENTER_LON_180='true '\n"
fi


#----------------------------------------------------------------
# Horizontal resolution-dependent settings
#-----------------------------------------------------------------

if [[ ${met} = "ModelE2.1" ]]; then
    if [[ "$runid" == "E213f10aF40oQ40nudge" ]]; then
        if [[ "$grid_res" ==  "4x5" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.00474046'\n"
        elif [[ "$grid_res" == "2x25" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.00243979'\n"
        elif [[ "$grid_res" == "05x0625" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.00276896'\n"
        elif [[ "$grid_res" == "025x03125" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.00254319'\n"
        else
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='-999.0e0'\n"
        fi
    else
        if [[ "$grid_res" ==  "4x5" ]]; then
            RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.03564873'\n"
        elif [[ "$grid_res" == "2x25" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.01050036'\n"
        elif [[ "$grid_res" == "05x0625" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.01340854'\n"
        elif [[ "$grid_res" == "025x03125" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.01066495'\n"
        else
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='-999.0e0'\n"
        fi
    fi
else
    RUNDIR_VARS+="RUNDIR_GISS_RES='not_used'\n"
    if [[ "x${sim_name}" == "xfullchem" || "x${sim_name}" == "xaerosol" ]]; then
      if [[ "x${met}" == "xgeosfp" && "x${grid_res}" == "x4x5" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='8.3286e-4'\n"
      elif [[ "x${met}" == "xgeosfp" && "x${grid_res}" == "x2x25" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='5.0416e-4'\n"
      elif [[ "x${met}" == "xmerra2" && "x${grid_res}" == "x4x5" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='7.8533e-4'\n"
      elif [[ "x${met}" == "xmerra2" && "x${grid_res}" == "x2x25" ]]; then
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='4.7586e-4'\n"
      else
          RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='-999.0e0'\n"
      fi
    else
      RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='-999.0e0'\n"
    fi
fi

RUNDIR_VARS+="RUNDIR_GRID_NLEV='$grid_lev'\n"

rundir=${MY_PATH}/${RUN_NAME}/${RUN_TEMPLATE}
mkdir -p ${rundir}/Restarts

# Define a subdirectory for rundir configuration files
rundir_config=${rundir}/CreateRunDirLogs
mkdir -p ${rundir_config}

# Copy run directory files and subdirectories
cp ${gcdir}/run/shared/cleanRunDir.sh       ${rundir}
cp ${gcdir}/run/shared/download_data.py     ${rundir}
cp ${gcdir}/run/shared/download_data.yml    ${rundir}
cp ${GCC_RUN_FILES}/getRunInfo              ${rundir}
cp ${GCC_RUN_FILES}/archiveRun.sh           ${rundir}
cp ${GCC_RUN_FILES}/README.md               ${rundir}
cp ${GCC_RUN_FILES}/gitignore               ${rundir}/.gitignore

# Use data downloader that points to GCAP2 restart files
if [[ ${met} = "ModelE2.1" ]] || [[ ${met} = "ModelE2.2" ]]; then
  cp ${gcdir}/run/shared/download_data.gcap2.40L.yml ${rundir}/download_data.yml
fi

if [[ "x${sim_name}" == "xfullchem" || "x${sim_name}" == "xCH4" ]]; then
    cp -r ${gcdir}/run/shared/metrics.py  ${rundir}
    chmod 744 ${rundir}/metrics.py
fi

# Set permissions
chmod 744 ${rundir}/cleanRunDir.sh
chmod 744 ${rundir}/archiveRun.sh


# Copy species database.  NOTE: For the new Hg simulation via KPP, we need
# to copy species_database_hg.yml to the rundir and rename it to
# species_database.yml.  This is because the Hg simulation has several
# inactive species that are active in the other simulations, and this
# causes a conflict.  Work out a better solution later.
#  -- Bob Yantosca, 10 Dec 2021
if [[ "x${sim_num}" == "x5" ]]; then
    cp -r ${gcdir}/run/shared/species_database_hg.yml ${rundir}/species_database.yml
else
    cp -r ${gcdir}/run/shared/species_database.yml ${rundir}
fi


# Append APM or TOMAS species to species database
# Also copy APM input files to the run directory
if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    cat ${gcdir}/run/shared/species_database_tomas.yml >> ${rundir}/species_database.yml
elif [[ ${sim_extra_option} =~ "APM" ]]; then
    cat ${gcdir}/run/shared/species_database_apm.yml >> ${rundir}/species_database.yml
    cp ${gcdir}/run/shared/apm_tmp.dat ${rundir}/apm_tmp.dat
    cp ${gcdir}/run/shared/input.apm   ${rundir}/input.apm
fi

# If benchmark simulation, put run script in directory
if [[ "x${sim_extra_option}" == "xbenchmark" ]]; then
    cp ${GCC_RUN_FILES}/runScriptSamples/geoschem.benchmark.run ${rundir}
    chmod 744 ${rundir}/geoschem.benchmark.run
fi

# Create symbolic link to code directory
ln -s ${wrapperdir} ${rundir}/CodeDir
ln -s ${wrapperdir}/run/GCHP/runScriptSamples ${rundir}/runScriptSamples

# Create build directory
mkdir ${rundir}/build
printf "To build GEOS-Chem type:\n   cmake ../CodeDir\n   cmake . -DRUNDIR=..\n   make -j\n   make install\n" >> ${rundir}/build/README

startdate='{DATE1}'
enddate='{DATE2}'
starttime='{TIME1}'
endtime='{TIME2}'

RUNDIR_VARS+="RUNDIR_SIM_START_DATE=$startdate\n"
RUNDIR_VARS+="RUNDIR_SIM_END_DATE=$enddate\n"
RUNDIR_VARS+="RUNDIR_SIM_START_TIME='$starttime'\n"
RUNDIR_VARS+="RUNDIR_SIM_END_TIME='$endtime'\n"


# Use monthly diagnostics by default
RUNDIR_VARS+="RUNDIR_HIST_TIME_AVG_DUR='00000100 000000'\n"
RUNDIR_VARS+="RUNDIR_HIST_TIME_AVG_FREQ='00000100 000000'\n"
RUNDIR_VARS+="RUNDIR_HIST_INST_DUR='00000100 000000'\n"
RUNDIR_VARS+="RUNDIR_HIST_INST_FREQ='00000100 000000'\n"
RUNDIR_VARS+="RUNDIR_HIST_MONTHLY_DIAG='1'\n"

# Turn on GEOS-Chem timers for benchmark simulations
if [[ "${sim_extra_option}" == "benchmark" ]]; then
    RUNDIR_VARS+="RUNDIR_USE_GCCLASSIC_TIMERS='true '\n"
else
    RUNDIR_VARS+="RUNDIR_USE_GCCLASSIC_TIMERS='false'\n"
fi


# Assign appropriate file paths and settings in HEMCO_Config.rc
if [[ ${met} = "ModelE2.1" ]]; then
    RUNDIR_VARS+="RUNDIR_DUSTDEAD_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_SEASALT_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_SOILNOX_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_DUST='false'\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_BIOVOC='false'\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='false'\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_SOILNOX='false'\n"
    RUNDIR_VARS+="RUNDIR_TOMAS_SEASALT='off'\n"
    RUNDIR_VARS+="RUNDIR_TOMAS_DUSTDEAD='off'\n"
    RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/gcap2_hemco.txt)\n"
else
    if [[ "${sim_extra_option}" == "benchmark" ]]; then
      RUNDIR_VARS+="RUNDIR_DUSTDEAD_EXT='on '\n"
      RUNDIR_VARS+="RUNDIR_SEASALT_EXT='on '\n"
      RUNDIR_VARS+="RUNDIR_SOILNOX_EXT='on '\n"
      RUNDIR_VARS+="RUNDIR_OFFLINE_DUST='false'\n"
      RUNDIR_VARS+="RUNDIR_OFFLINE_BIOVOC='false'\n"
      RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='false'\n"
      RUNDIR_VARS+="RUNDIR_OFFLINE_SOILNOX='false'\n"
      RUNDIR_VARS+="RUNDIR_TOMAS_SEASALT='off'\n"
      RUNDIR_VARS+="RUNDIR_TOMAS_DUSTDEAD='off'\n"
    else
      if [[ "${sim_extra_option}" == "marinePOA" ]]; then
          RUNDIR_VARS+="RUNDIR_SEASALT_EXT='on '\n"
          RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='false'\n"
          RUNDIR_VARS+="RUNDIR_TOMAS_SEASALT='off'\n"
      else
          RUNDIR_VARS+="RUNDIR_SEASALT_EXT='off'\n"
          if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
            RUNDIR_VARS+="RUNDIR_TOMAS_SEASALT='on '\n"
            RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='false'\n"
          else
            RUNDIR_VARS+="RUNDIR_TOMAS_SEASALT='off'\n"
            RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='true '\n"
          fi
      fi
      if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
          RUNDIR_VARS+="RUNDIR_TOMAS_DUSTDEAD='on '\n"
          RUNDIR_VARS+="RUNDIR_OFFLINE_DUST='false'\n"
      else
          RUNDIR_VARS+="RUNDIR_TOMAS_DUSTDEAD='off'\n"
          RUNDIR_VARS+="RUNDIR_OFFLINE_DUST='true '\n"
      fi
      RUNDIR_VARS+="RUNDIR_DUSTDEAD_EXT='off'\n"
      RUNDIR_VARS+="RUNDIR_SOILNOX_EXT='off'\n"
      RUNDIR_VARS+="RUNDIR_OFFLINE_BIOVOC='true '\n"
      RUNDIR_VARS+="RUNDIR_OFFLINE_SOILNOX='true '\n"
    fi
    RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/gmao_hemco.txt)\n"
fi

#--------------------------------------------------------------------
# Replace settings in config files with RUNDIR variables
#--------------------------------------------------------------------

# Save RUNDIR variables to a file in the rundirConfig folder
rundir_config_log=${rundir_config}/rundir_vars.txt
echo -e "$RUNDIR_VARS" > ${rundir_config_log}
#sort -o ${rundir_config_log} ${rundir_config_log}

cd ${rundir}

# Initialize run directory
bash ${srcrundir}/init_rd.sh ${rundir_config_log}

