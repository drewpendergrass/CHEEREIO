.. _Deploying the Ensemble:

Deploying the ensemble
==========

Because CHEEREIO requires many of the ensemble settings to be globally available to many different types of scripts and programs, written in different languages and stored in in different directories, it expects a very particular kind of installation process and for user settings to be made available in a strict format (so that they can be read by scripts in multiple languages). The overall install process is not too different from the way GEOS-Chem run directories are created as of version 13.0.0, but it does require some care on the part of the user. This page explains how to install the ensemble and customize it to meet your needs. Important science goes into this step!

.. _Configuration:

The Ensemble configuration file
-------------

All ensemble configuration is set by the ``ens_config.json`` file in the main CHEEREIO code directory. This file contains all the information that the user would normally input in the process of making a GEOS-Chem run directory, in addition to other settings like cluster configuration, assimilation global variables (like the localization radius), the composition of the state and control vectors, links to observation operators, and details about which emissions can be assimilated. This file is in JSON format, which is a file type that is human readable but expects a fairly strict syntax. 

**Important:** CHEEREIO constantly references the ``ens_config.json`` file throughout runtime. Changes to this file during runtime will change the settings of the ongoing assimilation calculation. However, there are some settings that are applied at the time of template run directory creation (e.g. cluster settings like memory usage) and cannot be adjusted at runtime via the ``ens_config.json`` file. **To be safe, you should redeploy the ensemble when you make changes to the configuration file and avoid making any changes during runtime.**

Let's walk through the configuration file section-by-section, after first covering some basic formatting requirements.

JSON formatting
~~~~~~~~~~~~~

JSON stores data in a form very similar to a Python dictionary. To have a basic idea of what CHEEREIO expects, let's take a look at the first ten lines from an example ``ens_config.json`` file:

::

	{
		"comment000" : "***************************************************************************",
		"comment001" : "****************BEGIN BASIC GEOS-CHEM AND ENSEMBLE SETTINGS****************",
		"comment002" : "***************************************************************************",
		"RES" : "4.0x5.0",
		"met_name" : "GEOSFP",
		"LEVS" : "47",
		"NEST" : "F",
		"REGION" : "",
		"BUFFER" : "",
		"ASSIM_PATH" : "/n/home12/drewpendergrass/CHEEREIO",
		"RUN_NAME" : "SHORT_FULL_TEST",
		"MY_PATH" : "/n/holyscratch01/jacob_lab/dpendergrass/GC-LETKF",
		"DATA_PATH" : "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData",
		"RESTART_FILE" : "/n/holyscratch01/jacob_lab/dpendergrass/GC-LETKF/input_data/GEOSChem.Restart.20190101_0000z.nc4",

The first line of the file is a left angle bracket, while the final line is a right angle bracket; inside are a set of keys associated (left side of the colon) associated with values (right side of the colon) by the colon operator.  CHEEREIO only obtains values by key reference. For example, it will always obtain the resolution of the run (in this case 4x5) by searching for the value associated with the "RES" key. Line number does not matter, but each line must end in a comma.

JSON does not support comments. The first three key-value pairs were chosen to approximate a comment, but make no difference at run time.

Values need not be strings; some of the CHEEREIO settings are supplied in the form of arrays or sub-dictionaries. Two examples from ``ens_config.json`` are shown below (and will be discussed in greater detail later):

::

	"CONTROL_VECTOR_CONC" : [
		"NO",
		"NO2",
		"HNO3",
		"HNO4",
		"PAN",
		"MPAN",
		"N2O5"
	],
	"CONTROL_VECTOR_EMIS" : {
		"NO_AGR":"NO",
		"NO_OTHER":"NO",
		"NO2":"NO2"
	},

**Important**: There is one subtlety in this particular configuration file: colons in the ``WallTime`` and ``SpinupWallTime`` entries must be escaped with two backslashes (``\\``). The first backslash escapes the second backslash in JSON; the second backslash escapes the colon in the SED unix utility which is used in the CHEEREIO installation process. For example, to allow a one day eight hour simulation, you would write ``"WallTime" : "1-08\\:00",``.

More details on the JSON format are available on the JSON `website <https://www.json.org>`__. The rest of this section will cover the various parts of the ``ens_config.json`` file and the settings they control.

Basic GEOS-Chem and ensemble settings
~~~~~~~~~~~~~

The first section of the ``ens_config.json`` file (i.e. between the first two comments) mostly controls settings analagous to those set during normal GEOS-Chem run directory creation. However, there are a few unique options in this setting particular to CHEEREIO. We'll consider these one-by-one.

* RES: The resolution of the GEOS-Chem model. Options are available on the `GEOS-Chem website <http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_horizontal_grids>`__ and include 4.0x5.0, 2.0x2.5, 0.5x0.625, 0.25x0.3125 and nested grid settings in format TwoLetterCode_MetCode (e.g. AS_MERRA2, EU_GEOSFP). Custom nested domains are not currently supported by the automated scaling factor creation utility but can be manually added by the user.
* met_name: Meteorology (chosen from MERRA2, GEOSFP, or ModelE2.1).
* LEVS: Number of levels (47 or 72).
* NEST: Is this a nested grid simulation? "T" or "F".
* REGION: Two letter region code for nested grid, or empty string ("") if not.
* ASSIM_PATH: **Full path** to the directory where the CHEEREIO repository is installed (e.g. ``/n/home12/drewpendergrass/CHEEREIO``). Directories in the ``ens_config.json`` file **should not have trailing forward slashes.**
* RUN_NAME: The name of the CHEEREIO ensemble run (will be the name of the folder containing the ensemble, template run directory, temporary files, and so on .
* MY_PATH: Path to the directory where ensembles will be created. A folder with name ``RUN_NAME`` will be created inside.
* DATA_PATH: Path to where external GEOS-Chem data is located. This can be an empty string if GEOS-Chem has already been configured on your machine (it is automatically overwritten).
* RESTART_FILE: Full path to the restart file for the simulation.
* BC_FILES: Full path to the boundary condition files for the simulation if a nested grid (empty string otherwise).
* sim_name: Simulation type. Valid options are "fullchem", "aerosol", "CH4", "CO2", "Hg", "POPs", "tagCH4", "tagCO", "tagO3", and "TransportTracers".
* chemgrid: Options are "trop+strat" and "trop_only".
* sim_extra_option: Options are "none", "benchmark", "complexSOA", "complexSOA_SVPOA", "marinePOA", "aciduptake", "TOMAS15", "TOMAS40", "APM", "RRTMG", "BaP", "PHE", and "PYR". Depending on the simulation type only some will be available. Consult the GEOS-Chem documation for more information.
* DO_SPINUP: Would you like CHEEREIO to set up a spinup directory for you? "true" or "false". The ensemble will automatically start from the end restart file produced by this run.
* SPINUP_START: Start date for spinup (YYYYMMDD). Empty string if no spinup.
* SPINUP_END: End date for spinup (YYYYMMDD).
* START_DATE: Start date for ensemble run (YYYYMMDD).
* ASSIM_START_DATE: Date where assimilation begins (YYYYMMDD). It's usually good to give a few days of spinup to create variations between ensemble members, since emissions drive the difference.
* END_DATE: End date for ensemble run (YYYYMMDD).
* nEnsemble: Number of ensemble members. 32 is usually a good number. This number of run directories will be created in the ``ensemble_runs`` folder and will be run simultaneously.
* pPERT: Range of initial emissions perturbation. For example, if "0.5" selected then scaling factors will initially range between 0.5 and 1.5 sampled from a uniform distribution.
* SIMULATE_NATURE: End users should almost always set this to "false", as this is used for testing. "true" or "false", should CHEEREIO generate an additional run directory for a run to be treated as nature (observation operators applied to create simulated observations).  

Cluster settings
~~~~~~~~~~~~~

The next section of the ``ens_config.json`` file controls settings that will be used when submitting jobs to the scheduler. These settings overwrite the template batch submission scripts included with CHEEREIO.

* NumCores: Number of cores used in each of the ensemble runs. CHEEREIO also will use these cores to parallelize assimilation computation columnwise.
* Partition: Partition of your cluster you are submitting to. At Harvard, ``huce_intel`` is a good choice.
* Memory: Memory in megabytes used by each ensemble member. CHEEREIO is quite memory intensive because it loads in restarts and history files for many ensemble members, so expect to use more than in standard GEOS-Chem runs.
* WallTime: Time allowed for the overall assimilation process (runs and assimilation) to occur in format D-HH\\:MM. Assimilation adds substantial overhead so expect it to be slow.
* SpinupWallTime: Wall time for the spinup simulation, if you're using one. Empty string otherwise.
* CondaEnv: The name of the Conda environment with all of the CHEEREIO packages installed. It is strongly recommended that you install an environment using the YAML file that ships with CHEEREIO.
* MaxPar: Maximum number of columns to assimilate in parallel using CHEEREIO, maxing out at NumCores. Setting this number smaller than NumCores saves on memory but adds to the assimilation time. 

Species in state/control/observation vectors
~~~~~~~~~~~~~

* STATE_VECTOR_CONC: Species from the restart files to be included in the state vector. It is generally recommended to include a fairly wide range of species that might affect the species you are mainly interested in, but not so large a range that you end up analyzing noise. Given as an array.
* CONTROL_VECTOR_CONC: A subset of the state vector concentration species that will be updated by assimilation. Although an update for all members of the state vector will be calculated, only these species will have that update saved. This allows a wide range of species to be considered in the update calculation process but only a smaller, more tightly coupled subset of species to actually be changed and passed to GEOS-Chem. The goal is to tamp down on noise. 
* CONTROL_VECTOR_EMIS: A dictionary linking a label for emissions scalings to the species emitted. For example, you could write ``"NO_AGR" : "NO"`` to reference agricultural NO emissions. CHEEREIO automatically will update ``HEMCO_Config.rc`` accordingly, but cannot distinguish between different emissions of the same species on its own; the user has to manually edit ``HEMCO_Config.rc`` to correct this if distinguishing between different sources of the same species. More on this in the `Template Run <template>`__ section.
* OBSERVED_SPECIES: A dictionary linking a label for observations with the species observed. For example, you could write ``"NO2_SATELLITE" : "NO2"`` to reference satellite observations of NO2. Unlike elsewhere, here the order matters. Later in the configuration file, arrays of observation operators and errors will be associated with these species according to the order they are stored. More in the next section. 


Miscellaneous LETKF settings
~~~~~~~~~~~~~

* LOCALIZATION_RADIUS_km: 1000",
* ALLOW_ADAPTIVE_LOCALIZATION": "false",
* NUM_OBS_TO_NUM_ENS_ADAPTIVE_MULTIPLIER: 5.0",
* ENABLE_4D_LETKF: false",
* OBS_ERROR_MATRICES": [
* OBS_OPERATORS": [
* NATURE_H_FUNCTIONS" : [
* INFLATION_FACTOR": "0.02",
* ASSIM_TIME": "24",
* TESTBIAS": "0.5"

The Setup Ensemble script
-------------

TKTKTK

.. _Template:

The Template Run Directory
-------------

TKTKTK

The Spinup Run Directory
-------------

TKTKTK

