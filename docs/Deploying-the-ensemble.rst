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

* RES: The resolution of the GEOS-Chem model. Options are available on the `GEOS-Chem website <http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_horizontal_grids>`__ and include 4.0x5.0, 2.0x2.5, 0.5x0.625, 0.25x0.3125 and nested grid settings in format TwoLetterCode\_MetCode (e.g. AS\_MERRA2, EU\_GEOSFP).
* met_name: GEOSFP",
* LEVS: 47",
* NEST: F",
* REGION: ",
* BUFFER: ",
* ASSIM_PATH: /n/home12/drewpendergrass/CHEEREIO",
* RUN_NAME: SHORT_FULL_TEST",
* MY_PATH: /n/holyscratch01/jacob_lab/dpendergrass/GC-LETKF",
* DATA_PATH: /n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData",
* RESTART_FILE: /n/holyscratch01/jacob_lab/dpendergrass/GC-LETKF/input_data/GEOSChem.Restart.20190101_0000z.nc4",
* BC_FILES: ",
* sim_name: fullchem",
* chemgrid: trop+strat",
* sim_extra_option: none",
* DO_SPINUP: false",
* SPINUP_START: ",
* SPINUP_END: ",
* START_DATE: 20190101",
* ASSIM_START_DATE: 20190108",
* END_DATE: 20190115",
* nEnsemble: 16",
* pPERT: 0.5",
* SIMULATE_NATURE: true",

Cluster settings
~~~~~~~~~~~~~

* NumCores: 8",
* Partition: huce_intel",
* Memory: 40000",
* WallTime: 1-08\\:00",
* SpinupWallTime: ",
* CondaEnv: cheerio",
* MaxPar: 3",

Species in state/control/observation vectors
~~~~~~~~~~~~~

* STATE_VECTOR_CONC" : [
* CONTROL_VECTOR_CONC" : [
* CONTROL_VECTOR_EMIS" : {
* OBSERVED_SPECIES" : {


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

