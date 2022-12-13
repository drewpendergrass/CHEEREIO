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

cp ${gcdir}/run/shared/cleanRunDir.sh ${RUN_TEMPLATE}
cp ${gcdir}/run/shared/download_data.py ${RUN_TEMPLATE}
cp ${GCC_RUN_FILES}/getRunInfo ${RUN_TEMPLATE}
cp ${GCC_RUN_FILES}/archiveRun.sh ${RUN_TEMPLATE}
cp ${GCC_RUN_FILES}/README ${RUN_TEMPLATE}
cp ${GCC_RUN_FILES}/gitignore ${RUN_TEMPLATE}/.gitignore
cp ${GCC_RUN_FILES}/input.geos.templates/input.geos.${sim_name} ${RUN_TEMPLATE}/input.geos
cp ${GCC_RUN_FILES}/HISTORY.rc.templates/HISTORY.rc.${sim_name} ${RUN_TEMPLATE}/HISTORY.rc
if [[ ${USE_CUSTOM_CH4} = "True" ]]; then
  cp ${ASSIM_PATH}/templates/HEMCO_Config.rc.CH4_template13 ${RUN_TEMPLATE}/HEMCO_Config.rc
else
  cp ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/HEMCO_Config.rc.${sim_name} ${RUN_TEMPLATE}/HEMCO_Config.rc
fi
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

if [[ ${USE_CUSTOM_CH4} = "True" ]]; then
      # Specify meteorology
  if [[ ${met_name} = "ModelE2.1" ]]; then
      sed_ie "/# --- Meteorology fields for FlexGrid ---/r ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/met_fields.gcap2" HEMCO_Config.rc
  else
      sed_ie "/# --- Meteorology fields for FlexGrid ---/r ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/met_fields.gmao"  HEMCO_Config.rc
  fi
else
    # Specify meteorology
  if [[ ${met_name} = "ModelE2.1" ]]; then
      sed_ie "/### BEGIN SECTION SETTINGS/r ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/header.gcap2"                    HEMCO_Config.rc
      sed_ie "/# --- Meteorology fields for FlexGrid ---/r ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/met_fields.gcap2" HEMCO_Config.rc
  else
      sed_ie "/### BEGIN SECTION SETTINGS/r ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/header.gmao"                     HEMCO_Config.rc
      sed_ie "/# --- Meteorology fields for FlexGrid ---/r ${GCC_RUN_FILES}/HEMCO_Config.rc.templates/met_fields.gmao"  HEMCO_Config.rc
  fi
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
if [[ ${USE_CUSTOM_CH4} = "True" ]]; then
  sed_ie "s|{CH4_HEMCO_ROOT}|${CH4_HEMCO_ROOT}|"   HEMCO_Config.rc
  sed_ie "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   HEMCO_Config.rc
else
  sed_ie "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   HEMCO_Config.rc
fi
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


exit 0
