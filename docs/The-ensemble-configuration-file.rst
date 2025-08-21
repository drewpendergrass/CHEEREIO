.. _Configuration:

Configuring your simulation
==========

All ensemble configuration is set by the ``ens_config.json`` file in the CHEEREIO subdirectory of your ensemble folder (copied at installation time from the directory you cloned from GitHub). This file contains all the information that the user would normally input in the process of making a GEOS-Chem run directory, in addition to other settings like cluster configuration, global assimilation variables (like the localization radius), the composition of the state and control vectors, links to observations, and details about which emissions can be assimilated. This file is in JSON format, which is a format that is human readable but enforces a strict syntax. 

Because CHEEREIO requires many of the ensemble settings to be globally available to many different types of scripts and programs, written in different languages and stored in in different directories, it expects a very particular kind of installation and run process and for user settings to be made available in a strict JSON format (so that they can be read by scripts in multiple languages). This page explains how to customize the ensemble to meet your needs. This isn't merely a technical process: the user makes important scientific assumptions in this step!

**Important:** CHEEREIO constantly references the ``ens_config.json`` file throughout runtime. Changes to this file during runtime will change the settings of the ongoing assimilation calculation. However, there are some settings that are applied at the time of template run directory creation (e.g. cluster settings like memory usage) and cannot be adjusted at runtime via the ``ens_config.json`` file. **To be safe, you should redeploy the ensemble when you make changes to the configuration file and avoid making any changes during runtime.**

Let's walk through the configuration file section-by-section, after first covering some basic formatting requirements.

JSON formatting
-------------

JSON stores data in a form very similar to a Python dictionary. To have a basic idea of what CHEEREIO expects, let's take a look at the first ten lines from an example ``ens_config.json`` file:

::

	{
		"comment000" : "***************************************************************************",
		"comment001" : "****************BEGIN BASIC GEOS-CHEM AND ENSEMBLE SETTINGS****************",
		"comment002" : "***************************************************************************",
		"RES" : "4.0x5.0",
		"met_name" : "MERRA2",
		"LEVS" : "47",
		"NEST" : "F",
		"REGION" : "",
		"ASSIM_PATH" : "/n/home12/drewpendergrass/CHEEREIO_no2",
		"RUN_NAME" : "NOx_GLOBAL_v2",
		"MY_PATH" : "/n/holyscratch01/jacob_lab/dpendergrass/GC-LETKF",
		"DATA_PATH" : "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/gcdata/ExtData",
		"CH4_HEMCO_ROOT" : "/n/seasasfs02/CH4_inversion/InputData/HEMCO",
		"USE_CHEEREIO_TEMPLATE_CH4_HEMCO_Config" : "False",
		"RESTART_FILE" : "/n/holyscratch01/jacob_lab/dpendergrass/GC-LETKF/input_data/GEOSChem.Restart.20190101_0000z.nc4",

The first line of the file is a left curly bracket, while the final line is a right curly bracket; inside are a set of keys associated (left side of the colon) associated with values (right side of the colon) by the colon operator.  CHEEREIO only obtains values by key reference. For example, it will always obtain the resolution of the run (in this case 4x5) by searching for the value associated with the "RES" key. Line number does not matter, but each line must end in a comma.

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
		"NOx":["NO","NO2"]
	},

However, CHEEREIO does expect all values in ``ens_config.json`` as strings. If a setting is given by a number, wrap that number in quotation marks.  

**Important**: There is one subtlety in this particular configuration file: colons in the ``WallTime`` and ``SpinupWallTime`` entries (covered later) must be escaped with two backslashes (``\\``). The first backslash escapes the second backslash in JSON; the second backslash escapes the colon in the SED unix utility which is used in the CHEEREIO installation process. For example, to allow a one day eight hour simulation, you would write ``"WallTime" : "1-08\\:00",``.

**Important**: Follow the capitalization conventions in the template ``ens_config.json`` file! CHEEREIO is case sensitive.

More details on the JSON format are available on the JSON `website <https://www.json.org>`__. When in doubt, follow the conventions in the template ``ens_config.json`` files!

A line-by-line guide to ensemble configuration
-------------

The rest of this section will cover the various parts of the ``ens_config.json`` file and the settings they control. For a first simulation, it's usually not a bad idea to follow the settings in the template ``ens_config.json`` files.


Basic GEOS-Chem and ensemble settings
~~~~~~~~~~~~~

The first section of the ``ens_config.json`` file (i.e. between the first two comments) mostly controls settings analagous to those set during normal GEOS-Chem run directory creation. However, there are a few unique options in this setting particular to CHEEREIO. We'll consider these one-by-one.

.. option:: GC_VERSION
	
	GEOS-Chem model version (e.g. "14.1.1"). Different behaviors are required for different model versions, so this field is essential.

.. option:: RES
	
	The resolution of the GEOS-Chem model. Options are available on the `GEOS-Chem website <http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_horizontal_grids>`__ and include 4.0x5.0, 2.0x2.5, 0.5x0.625, 0.25x0.3125 and nested grid settings in format TwoLetterCode_MetCode (e.g. AS_MERRA2, EU_GEOSFP). Custom nested domains are not currently supported by the automated scaling factor creation utility but can be manually added by the user. If there is enough interest I will add more automated support in a later CHEEREIO update. 

	.. attention::

		**Once the CHEEREIO ensemble is installed, the resolution cannot be changed without re-installing the ensemble.** CHEEREIO sets up a number of LETKF-specific routines assuming a specific resolution (including the parallelization design), and will fail if the resolution is adjusted after this setup is complete.

.. option:: met_name
	
	Meteorology (chosen from MERRA2, GEOSFP, or ModelE2.1).

.. option:: LEVS
	
	Number of levels (47 or 72).

.. option:: NEST
	
	Is this a nested grid simulation? "T" or "F".

.. option:: REGION
	
	Two letter region code for nested grid, or empty string ("") if not.

.. option:: ASSIM_PATH
	
	**Full path** to the directory where the CHEEREIO repository was originally installed (e.g. ``/n/home12/drewpendergrass/CHEEREIO``). Directories in the ``ens_config.json`` file **should not have trailing forward slashes.** Again, when in doubt follow the provided templates.

	.. attention::

		CHEEREIO at runtime will not reference this directory, but rather the version which was copied into your ensemble folder (same level as ``ensemble_runs/`` or ``template_run/``) including that copy of the ens_config.json file. This entry is used only at the beginning of the installation process.

.. option:: RUN_NAME
	
	The name of the CHEEREIO ensemble run (will be the name of the folder containing the ensemble, template run directory, temporary files, and so on).

.. option:: MY_PATH
	
	Path to the directory where ensembles will be created. A folder with name ``RUN_NAME`` will be created inside.

.. option:: DATA_PATH
	
	Path to where external GEOS-Chem data is located. This can be an empty string if GEOS-Chem has already been configured on your machine (it is automatically overwritten).

.. option:: CH4_HEMCO_ROOT
	
	If the subsequent option, "USE_CHEEREIO_TEMPLATE_CH4_HEMCO_Config", is set to "True", then this is the root folder where emissions and other input files for the methane specialty simulation are located. In this case, a special CHEEREIO ``HEMCO_Config.rc`` template from the ``templates/`` folder in the code directory is used. 

	.. attention::

		This option is functional but currently causes GEOS-Chem crashes with an unknown cause.

.. option:: USE_CHEEREIO_TEMPLATE_CH4_HEMCO_Config

	See above entry.

.. option:: RESTART_FILE
	
	Full path to the restart file for the simulation. If in the initialization process you selected ``SetupSpinupRun=true``, then this restart file will be used for the classic spin up routine (getting realistic atmospheric conditions for the entire ensemble). Otherwise, this will be the restart file used to initialize all ensemble members.

.. option:: BC_FILES
	
	Full path to the boundary condition files for the simulation if you are using a nested grid (empty string otherwise).

.. option:: sim_name
	
	Simulation type. Valid options are "fullchem", "aerosol", "CH4", "CO2", "Hg", "POPs", "tagCH4", "tagCO", "tagO3", and "TransportTracers".

.. option:: chemgrid
	
	Options are "trop+strat" and "trop_only".

.. option:: sim_extra_option
	
	Options are "none", "benchmark", "complexSOA", "complexSOA_SVPOA", "marinePOA", "aciduptake", "TOMAS15", "TOMAS40", "APM", "RRTMG", "BaP", "PHE", and "PYR". Depending on the simulation type only some will be available. Consult the GEOS-Chem documation for more information.

.. option:: DO_SPINUP
	
	Would you like CHEEREIO to set up a spinup directory for you? "true" or "false". The ensemble will automatically start from the end restart file produced by this run. Note this option is for the standard GEOS-Chem spinup (run once for the whole ensemble). Note that if this is activated, you have to run the ``setup_ensemble.sh`` utility with the ``SetupSpinupRun`` switch set to ``true``.

.. option:: SPINUP_START
	
	Start date for spinup (YYYYMMDD). Empty string if no spinup.

.. option:: SPINUP_END
	
	End date for spinup (YYYYMMDD).

.. option:: DO_CONTROL_RUN
	
	The control run is a normal GEOS-Chem simulation without any assimilation. The output of this simulation can be compared with the LETKF results in the postprocessing workflow. Set to "true" if using a control run (most users). There are two ways of doing control runs in CHEEREIO, which are detailed in the next entry on this page and on :ref:`control simulation` page.

.. option:: DO_CONTROL_WITHIN_ENSEMBLE_RUNS
	
	CHEEREIO has two ways of running a control simulation. The preferred method, which is activated by setting this option to true, is to run the control simulation as an additional ensemble member with label 0 (ensemble members used for assimilation are numbered starting at 1). This allows the control simulation to match non-assimilation adjustments performed on the ensemble, such as scaling concentrations to be non-biased relative to observations. The control directory in this case is created automatically when the ``setup_ensemble.sh`` utility is used to create the ensemble. If this option is set to false, and DO_CONTROL_RUN is set to true, then the control simulation is created as an additional run directory at the top directory level (analagous to :ref:`spinup simulation`). This keeps the control simulation fully separate from the ensemble and any non-assimilation adjustments that are performed. In this case, you have to run the ``setup_ensemble.sh`` utility with the ``SetupControlRun`` switch set to ``true`` to create the control run directory. More information is available on :ref:`control simulation` page.

.. option:: CONTROL_START
	
	Start date for the control run (YYYYMMDD). (Unnecessary if ``DO_CONTROL_WITHIN_ENSEMBLE_RUNS`` is set to ``true``)

.. option:: CONTROL_END
	
	End date for the control run (YYYYMMDD).

.. option:: DO_ENS_SPINUP
	
	Do you want to use a separate job array to spin up your GEOS-Chem ensemble with randomized scaling factors applied to each ensemble member? "true" or "false". If set to "true", shell scripts entitled ``run_ensemble_spinup_simulations.sh`` and ``run_ensspin.sh`` are installed in the ``ensemble_runs/`` folder. The user should then execute ``run_ensspin.sh`` to spin up the ensemble and create variability between ensemble members before executing ``run_ens.sh`` in the normal run procedure. For more information on the ensemble spinup process, see :ref:`Run Ensemble Spinup Simulations`.

.. option:: ENS_SPINUP_FROM_BC_RESTART
	
	It is possible to start the ensemble spinup procedure using a boundary condition file, rather than a traditional restart file. Set to "true" if using a BC file, and "false" if using a normal restart file to start the ensemble spinup.

.. option:: ENS_SPINUP_START
	
	Start date for ensemble spinup run (YYYYMMDD).

.. option:: ENS_SPINUP_END
	
	End date for ensemble spinup run (YYYYMMDD).

.. option:: START_DATE
	
	Start date for main ensemble data assimilation run (YYYYMMDD).

.. option:: ASSIM_START_DATE
	
	Date where assimilation begins (YYYYMMDD). This option allows you to run the first assimilation period for an extra long time (although the assimilation window remains the same), effectively providing an ensemble-wide spinup. For more information on this ensemble spinup option, see :ref:`Run Ensemble Spinup Simulations`. If you have set ``DO_ENS_SPINUP`` to ``true``, then you should set this date to be one assimilation window later than ``START_DATE``.

.. option:: SIMPLE_SCALE_FOR_FIRST_ASSIM_PERIOD
	
	At the end of the first assimilation period, rather than doing the full LETKF calculation, CHEEREIO can scale the ensemble mean so that it matches the observational mean. This is done because if the model is biased relative to observations the LETKF will perform suboptimal updates. Set to "true" to do this scaling (recommended) or "false" to do the usual LETKF calculation. For more information, see :ref:`Simple scale`.

.. option:: END_DATE
	
	End date for ensemble run (YYYYMMDD).

.. option:: AMPLIFY_ENSEMBLE_SPREAD_FOR_FIRST_ASSIM_PERIOD
	
	At the end of the ensemble spinup period, the spread in ensemble members may still not be great enough for species in the state vector. If this option is set to "true", then CHEEREIO will multiply the standard deviation of the ensemble after ensemble spinup is complete by the factor given in ``SPREAD_AMPLIFICATION_FACTOR``. This spread amplification is done after the first assimilation period, so it will work with either spinup method. For more information, see :ref:`Spread amplification`.

.. option:: SPREAD_AMPLIFICATION_FACTOR
	
	If ``AMPLIFY_ENSEMBLE_SPREAD_FOR_FIRST_ASSIM_PERIOD`` is set to "true", then this is the factor with which CHEEREIO will multiply the ensemble standard deviation at the end of the ensemble spinup period. For more information on the burn-in period, see :ref:`Burn in period`.

.. option:: species_to_amplify_not_in_statevec

	If ``AMPLIFY_ENSEMBLE_SPREAD_FOR_FIRST_ASSIM_PERIOD`` is set to "true", then amplify the species listed here even if they aren't included in the statevector. This can be useful if you are assimilating emissions only but still would like to amplify the spread of observed species. 

.. option:: DO_BURN_IN

	Should CHEEREIO do a burn-in period? "true" or "false." A burn-in period is a time period where full LETKF assimilation is being applied, but the results will be discarded from final analysis. The idea of a burn in period is to allow CHEEREIO's emissions to "catch up" with the system, as it takes time for the updated emissions in CHEEREIO to become consistent with observations. For more information on the burn-in period, see :ref:`Burn in period`.

.. option:: SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD
	
	If this option and ``DO_BURN_IN`` are both set to "true", then at the end of the burn-in period (a date given by ``BURN_IN_END``) CHEEREIO will scale the ensemble mean to match the observational mean, as in the ``SIMPLE_SCALE_FOR_FIRST_ASSIM_PERIOD`` option. This ensures that any biases introduced in the period where CHEEREIO emissions are "catching up" with observations are corrected. For more information on the burn-in period, see :ref:`Burn in period`.

.. option:: BURN_IN_END
	
	If ``DO_BURN_IN`` is set to ``true``, then this is the date (YYYYMMDD) when the burn-in period ends.

.. option:: POSTPROCESS_START_DATE
	
	The date when the postprocessing script should start (YYYYMMDD). This should always be at least one assimilation window away from ``START_DATE``. If you are using a burn-in period, you can set this for after the burn-in period ends to ensure that all your analysis discards this period.

.. option:: POSTPROCESS_END_DATE
	
	The date when the postprocessing script should end (YYYYMMDD). Usually the same as ``END_DATE``, though the user can change the postprocess start and end dates to fit whatever application they are interested in.

.. option:: nEnsemble
	
	Number of ensemble members. 24, 32, or 48 are usually good numbers. This number of run directories will be created in the ``ensemble_runs`` folder and will be run simultaneously.

.. option:: verbose
	
	Amount of information to print out as the ensemble runs (in the letkf*out files in the ensemble_runs/logs/ folder). 1 is the default. 0 supresses most output, 2 is useful for basic debugging, and 3 for intense debugging. 


Cluster settings
~~~~~~~~~~~~~

The next section of the ``ens_config.json`` file controls settings that will be used when submitting jobs to the scheduler. These settings overwrite the template batch submission scripts included with CHEEREIO.

.. option:: NumCores
	
	Number of cores used in each of the ensemble runs. CHEEREIO also will use these cores to parallelize assimilation computation columnwise.

.. option:: NumCtrlCores
	
	Number of cores to use in the control run simulation, if using.

.. option:: Partition
	
	Partition of your cluster you are submitting to. At Harvard, ``huce_intel`` is a good choice.

.. option:: Memory
	
	Memory in megabytes used by each ensemble member. CHEEREIO can be quite memory intensive because it loads in restarts and history files for many ensemble members in addition to observations, and sometimes produces large matrices, so expect to use more than in standard GEOS-Chem runs.

.. option:: EnsCtrlSpinupMemory
	
	Memory in megabytes for ensemble spinup, control, and regular spinup simulations (i.e. those simulations without LETKF assimilation). Set as you would a normal GEOS-Chem simulation.

.. option:: WallTime
	
	Time allowed for the overall assimilation process (runs and assimilation) to occur in format D-HH\\\\:MM. Assimilation adds substantial overhead so expect it to be slow.

.. option:: EnsSpinupWallTime
	
	Time allowed for the ensemble spinup process (no assimilation, just running all ensemble members from ``ENS_SPINUP_START`` through ``ENS_SPINUP_END`` with scaling factors applied) in format D-HH\\\\:MM. If not using, you can just leave as an empty string.

.. option:: ControlWallTime
	
	Wall time for the control run simulation, if you're using one separate from the ensemble itself. Empty string otherwise.

.. option:: SpinupWallTime
	
	Wall time for the spinup simulation, if you're using one. Empty string otherwise.

.. option:: CondaEnv
	
	The name of the Conda environment with all of the CHEEREIO packages installed. It is strongly recommended that you install an environment using the ``cheereio.yaml`` file that ships with CHEEREIO in the ``environments/`` folder.

.. option:: AnimationEnv
	
	The name of the Conda environment that has the tools necessary to make animated postprocessing plots. It is strongly recommended that you install an environment using the ``animation.yml`` file that ships with CHEEREIO in the ``environments/`` folder.

.. option:: MaxPar
	
	Maximum number of cores to use while assimilating columns in parallel using CHEEREIO, maxing out at ``NumCores``. Setting this number smaller than NumCores saves on memory but adds to the assimilation time. 


Species in state/control vectors
~~~~~~~~~~~~~

One useful feature of CHEEREIO is its distinction between "control" and "state" vectors. The state vector should consist of all concentrations relevant to the problem at hand as well as the emissions of interest (e.g. large chemical families). The control vector is a subset of the state vector, and represents concentrations and the same emissions of interest that the user believes can reasonably be updated on the basis of observations. In many cases the control vector and state vector are identical. However, in some cases removing species from the control vector can help CHEEREIO handle ensemble behavior. Although the entire state vector is used to calculate the concentration and emissions update, **only the control vector is actually updated.** In practice, this distinction helps tamp down on noise and create well-behaved assimilations.

.. option:: STATE_VECTOR_CONC
	
	Species from the restart files to be included in the state vector. It is generally recommended to include a range of species that might affect the species you are mainly interested in, but not so large a range that you end up analyzing noise. Given as an array. Note that this can be left empty to do a pure emissions optimization without adjusting concentrations. This is an example for NO\ :sub:`x` data assimilation: 
	::

		"STATE_VECTOR_CONC" : [
			"NO",
			"NO2",
			"HNO3",
			"HNO4",
			"PAN",
			"MPAN",
			"N2O5"
		],


.. option:: CONTROL_VECTOR_CONC
	
	A subset of the state vector concentration species that will be updated by assimilation. Although an update for all members of the state vector will be calculated, only the species listed in this array will have that update saved. This allows a wide range of species to be considered in the update calculation process but only a smaller, more tightly coupled subset of species to actually be changed and passed to GEOS-Chem. The goal is to tamp down on noise. Again, in many simulations the state vector and control vector entries will be identical. Note that this can be left empty to do a pure emissions optimization without adjusting concentrations.

.. option:: STATE_VECTOR_CONC_REPRESENTATION
	
	How are concentrations represented within the state vector? There are several options. "3D" puts all 3D concentrations of the species in ``STATE_VECTOR_CONC`` from the restart file into the state vector. ``column_sum`` computes the partial columns of the species of interest and then sums to obtain units molec/cm2. Finally, ``trop_sum`` is identical to ``column_sum`` except that it only represents the tropospheric column. 

.. option:: CONTROL_VECTOR_EMIS
	
	A dictionary linking a label for emissions scalings to the species emitted. For example, you could write ``"CH4_WET" : "CH4"`` to reference wetland methane emissions. CHEEREIO automatically will update ``HEMCO_Config.rc`` accordingly, but cannot distinguish between different emissions of the same species on its own; the user has to manually edit ``HEMCO_Config.rc`` to correct this if distinguishing between different sources of the same species. For a more thorough explanation, see the entry in :ref:`Template`. You can also use one label to link to emissions of multiple species, meaning that all these emissions will be controlled by one scaling factor file, such as ``"NOx":["NO","NO2"]`` to indicate that NO and NO\ :sub:`2`\ are controlled by one scaling factor file. Here are a couple of examples:
	::

		"CONTROL_VECTOR_EMIS" : {
			"NOx":["NO","NO2"]
		},

		"CONTROL_VECTOR_EMIS" : {
			"CH4_WET":"CH4",
			"CH4_OTHER":"CH4"
		},


HISTORY.rc settings
~~~~~~~~~~~~~

.. option:: HISTORY_collections_to_customize
	
	A list of collections under HISTORY.rc that CHEEREIO will customize with user-specified frequency and duration settings. Here is a typical example:
	::

		"HISTORY_collections_to_customize" : [
			"SpeciesConc",
			"LevelEdgeDiags",
			"StateMet"
		],

.. option:: HISTORY_freq
	
	Frequency of data saved within history output files listed within collections in ``HISTORY_collections_to_customize``. For more information on history frequencies, see the GEOS-Chem manual.

.. option:: HISTORY_dur
	
	As in ``HISTORY_freq``, but for duration.

.. option:: SPINUP_HISTORY_freq
	
	Frequency of history output files saved during ensemble spinup (i.e. when executing ``run_ensspin.sh``). Often set to be a longer period to save memory.

.. option:: SPINUP_HISTORY_dur
	
	As in ``SPINUP_HISTORY_freq``, but for duration.

.. option:: SaveLevelEdgeDiags
	
	Should the LevelEdgeDiags collection be turned on? "True" or "False". This is mandatory for assimilating most forms of satellite data.

.. option:: SaveStateMet
	
	Should the StateMet collection be turned on? "True" or "False". This is mandatory for assimilating some forms of satellite data, like OMI NO\ :sub:`2`\ . 
 
.. option:: SaveArea
	
	Should grid cell areas be used in the assimilation process? "True" or "False".

.. option:: SaveSatDiagn
	
	Should the SatDiagn collection be turned on? "True" or "False". This setting is useful for some developmental satellite operators, like for CrIS. 

.. option:: HistorySpeciesConcToSave
	
	A list of species to save in the SpeciesConc collection. At minimum, this should encompass the concentration portion of the state vector and any concentrations needed for observation operators. Below is an example: 
	::

		"HistorySpeciesConcToSave" : [
			"NO",
			"NO2",
			"HNO3",
			"HNO4",
			"PAN",
			"MPAN",
			"N2O5"
		],

.. option:: HistoryLevelEdgeDiagsToSave
	
	A list of data to save in the LevelEdgeDiags collection. Just ``Met_PEDGE`` is sufficient for many forms of assimilation.

.. option:: HistoryStateMetToSave
	
	A list of data to save in the StateMet collection. Below is an example of necessary fields for assimilating OMI NO\ :sub:`2`\ .
	::

		"HistoryStateMetToSave" : [
			"Met_TropLev",
			"Met_BXHEIGHT",
			"Met_T"
		],

.. option:: HistoryXXXXXToSave
	
	For every collection in ``HISTORY_collections_to_customize`` that has not already been discussed (SpeciesConc, LevelEdgeDiags, StateMet, SatDiagn), the user can create a new entry in the configuration file following the pattern of ``HistoryLevelEdgeDiagsToSave`` or ``HistoryStateMetToSave`` but with the names and entries of their collection of interest. For example, if the user included ``Restart`` in ``HISTORY_collections_to_customize``, they might write the following:
	::

		"HistoryRestartToSave" : [
			"SpeciesRst_?ALL?",
			"Met_PS1DRY",
			"Met_TMPU1",
			"Met_BXHEIGHT",
			"Met_TropLev"
		],

.. _Observation settings:

Observation settings
~~~~~~~~~~~~~

.. option:: OBSERVED_SPECIES
	
	A dictionary linking a label for observations with the species observed. For example, you could write ``"NO2_SATELLITE" : "NO2"`` to reference satellite observations of NO2. 

.. option:: OBS_TYPE
	
	A dictionary linking a label for observations with the observer type, so that CHEEREIO knows how to interpret observation files. One entry is required for each entry in ``OBSERVED_SPECIES``, with every key from ``OBSERVED_SPECIES`` represented here. Valid values include "OMI" and "TROPOMI", or any other observation operator type added to CHEEREIO by you or by other users listed in ``operators.json``. Instructions on how to add an observation operator to CHEEREIO such that it can be switched on from ``OBS_TYPE`` in the configuration file are given in the :ref:`New observation:` page.

.. option:: ASSIMILATE_OBS

	A dictionary linking a label for observations (key) with a switch ("True" or "False") indicating whether or not those observations will be used in the LETKF calculation (value). If observations are not used in the LETKF (i.e. set to "False"), then CHEEREIO will use the observations only for making postprocessing plots and calculations; the user can treat these observations as an external data source for validation.

	.. attention::

		*Note: Even if a value of False is given for a given observer, entries for that observer still need to be given in any dictionary setting referencing OBSERVED_SPECIES. In many cases (e.g. error specification) the values will be ignored.*

.. option:: TROPOMI_dirs
	
	Dictionary linking observed TROPOMI species to the directory containing the observations. If you aren't using TROPOMI, this can be left blank. Here is an example, along with the corresponding ``OBSERVED_SPECIES`` settings:
	::

		"OBSERVED_SPECIES" : {
			"CH4_TROPOMI": "CH4"
		},
		"TROPOMI_dirs" : {
			"CH4" : "/n/holylfs05/LABS/jacob_lab/dpendergrass/tropomi/NO2/2019"
		},

.. option:: OMI_dirs
	
	As in TROPOMI_dirs, but for OMI. 

.. option:: XXXXXXXX_dirs

	With the exception of ObsPack observations (which are handled separately), users can link directories for any observer listed in ``operators.json`` using the above pattern from TROPOMI_dirs (e.g. OMI_dirs). Note that other observation operators can be added as separate key entries in this configuration file by following the instructions on the :ref:`New observation:` page.

.. option:: filter_obs_poleward_of_n_degrees

	A dictionary linking a label for observations with the maximum poleward extent of observations allowed. For example, a value of 60 would exclude all observations polewards of 60 degrees latitude. One entry is required for each entry in ``OBSERVED_SPECIES``, with every key from ``OBSERVED_SPECIES`` represented here. Set a given value to ``nan`` to ignore. 

.. option:: SaveDOFS
	
	Should CHEEREIO calculate and save the Degrees of Freedom for Signal (DOFS), or the trace of the observing system averaging kernel matrix? Note that since the prior error covariance matrix is not invertible because of our ensemble approach the pseudoinverse is used instead. See section 11.5.3 of Brasseur and Jacob for more information. The idea here is that if there is not enough information in a localized assimilation calculation we should set the posterior equal to the prior. 

	.. attention::

		*Note: this option is functional but DOFS values are not easily interpretable. In CHEEREIO 1.5, we will include an information metric from Voshtani et al., 2025 instead (see the citations page).*

.. option:: DOFS_filter
	
	What is the minimum DOFS for a localized region for the assimilation to be saved? If DOFS is below this threshold the posterior is set equal to the prior.

	.. attention::

		*Per previous note, leave set to "nan".*

.. option:: ACTIVATE_OBSPACK
	
	``true`` or ``false``, should we activate ObsPack? Set to ``true`` if and only if ObsPack is being used as an observed species in the assimilation (i.e. if and only if "OBSPACK" is listed as a value somewhere in the ``OBS_TYPE`` dictionary). Subsequent ObsPack entries are ignored if set to ``false``. More information on the ObsPack operator is given in the ``ObsPack tools`` entry.

.. option:: preprocess_raw_obspack_files

	``true`` or ``false``, would you like CHEEREIO to automatically process raw obspack observation files (as downloaded from NOAA) so that they are compatible to input into both CHEEREIO and GEOS-Chem? 


.. option:: raw_obspack_path
	
	If ``preprocess_raw_obspack_files`` is set to ``true``, then supply the path to the location of raw obspack files, without any preprocessing applied, as downloaded from NOAA. CHEEREIO will handle these directly.

	.. attention::

		*Note: Lee Murray reports that NOAA doesn't use consistent fields across ObsPack versions. Manual adjustments may be necessary, as discussed in the below entry.*

.. option:: gc_obspack_path
	
	Path where CHEEREIO and GEOS-Chem should find preprocessed ObsPack observation files which are compatible with GEOS-Chem. If ``preprocess_raw_obspack_files`` is ``true``, this directory can be empty and will be populated automatically by CHEEREIO during a preprocessing step at ensemble installation time.

	.. attention::

		*Note: Lee Murray reports that NOAA doesn't use consistent fields across ObsPack versions. If you get an error the preprocessing step (performed during ensemble run directory creation in the installation workflow), or you already have ObsPack files processed, you should set the preprocess_raw_obspack_files variable to false and supply an already populated directory of manually preprocessed files. Details for how to do this are provided in the ObsPack documentation for GEOS-Chem.*

.. option:: obspack_gc_input_file
	
	The file naming convention for obspack input data. This is used (1) by CHEEREIO in the pre-processing step when generating GEOS-Chem compatible ObsPack input files, and (2) by GEOS-Chem in the ObsPack diagnostic generation. For example, ``obspack_ch4.YYYYMMDD.nc`` might be appropriate for a methane simulation; the important part is that YYYYMMDD be present somewhere in the string.

.. option:: HistoryObsPackToSave
	
	A list of data to save in the ObsPack collection.Below is an example used for a methane experiment:
	::

		"HistoryObsPackToSave" : [
			"pressure",
			"CH4"
		],


.. _Scaling factor settings:

Scaling factor settings
~~~~~~~~~~~~~

.. option:: sf_initialization

	These settings control how scaling factors are intialized and how their errors are structured. It is implemented as a dictionary, and every key in ``CONTROL_VECTOR_EMIS`` must be represented here as a key. The value corresponding to these keys are another dictionary with the following required entries. 

	.. option:: sf_dim

		Must have value "1D" or "2D". Most users will pick "2D," which initalizes scaling factor errors as having some spatial dimension for both latitude and longitude. The "1D" option  initializes with error along the latitude dimension only, which may be more appropriate for cases where the user wants to optimize prescribed loss fields (like OH for methane). 

	.. option:: init_std_dev
		
		Setting for initial emissions scaling factor creation, where the number provided :math:`\sigma` is used to generate random scaling factors from a distribution specified in the ``lognormalErrors`` and ``correlatedInitialScalings`` entries. This will likely be a float less than 0.5 and greater than 0.0.

	.. option:: correlatedInitialScalings
		
		There are two possible interpretations of ``init_std_dev`` depending on the corresponding value in ``correlatedInitialScalings``. 

		.. option:: correlatedInitialScalings equals False
		
			If the corresponding value in ``correlatedInitialScalings`` is ``False``, then the  ``init_std_dev`` entry :math:`\sigma` is used to generate random scaling factors from the distribution :math:`n{\sim}N(1,\sigma)`, meaning that n is a normal random variable with mean 0 and standard deviation :math:`\sigma`. For example, if ``init_std_dev`` is "0.25" then scalings will have mean 1 and standard devation 0.25.  In this case, there is no correlation between neighboring scale factors (i.e. they are i.i.d.). If the user set to sf_dim to be 1D, then all scaling factors along a given latitude band will be equal but will be uncorrelated with scaling factors in other latitude bands.

		.. option:: correlatedInitialScalings equals True

			The other "std" option, where ``correlatedInitialScalings`` is ``True``, means that we sample scaling factors from a multivariate normal distribution. If ``speedyCorrelationApprox`` is turned off, then this is a true sampling from a multivariate normal distribution. The mean is a vector of ones and the covariance matrix is generated with exponentially decaying correlation as a function of distance. More specifically, the covariance between points :math:`a` and :math:`b` with distance :math:`d` kilometers between them is given by :math:`\exp(-d^2/(2*c))` where :math:`c` is a correlation distance constant in kilometers given by the ``corrDistances`` entry. However, in the 2D case, sampling a multivariate normal distribution like this can take an extraordinarily large amount of time and memory. CHEEREIO includes a very fast approximation to a multivariate normal distribution which is generated by applying a Gaussian blur to uncorrelated noise. To use this option, within the ``2D_settings`` dictionary, set ``speedyCorrelationApprox`` to ``True`` (strongly recommended for any resolution higher than 4x5). If the user has set ``sf_dim`` to "1D", then all sampling will be of the true normal distribution as it is quite cheap in the 1D case.

	.. option:: corrDistances
		
		See above entry for ``correlatedInitialScalings``.

	.. option:: lognormalErrors
		
		CHEEREIO supports lognormal errors for scaling factors, which can be activated using this setting. Emissions often follow a lognormal distribution, making this a reasonable choice of representation; it also naturally enforces a zero floor for emissions without squashing the ensemble spread. If activated, CHEEREIO will (1) sample the normal distributions mentioned in the entry for ``correlatedInitialScalings`` with a location parameter (mean) of 0, and then (2) take the exponential of the sample (which will then center on 1). GEOS-Chem will use these lognormally distributed samples in the physical model run, but the scaling factors will be log-transformed back into a normal distribution with mean 0 for LETKF calculations. They will again be transformed via an exponential into lognormal space for the next GEOS-Chem run. A suitable entry for an approximated multivariate lognormal sample is given below:
		::

			"sf_dim" : "2D",
			"init_std_dev" : 0.1",
			"correlatedInitialScalings" : "True",
			"corrDistances":"500",
			"lognormalErrors" : "True",
			...
			"2D_settings":{
				"speedyCorrelationApprox" : "True",
				...
			}
			

	.. option:: MaskOceanScaleFactor
		
		Should scaling factors be allowed to vary over the oceans? If "True", scaling factors for that species over the ocean are always set to 1 across all ensemble members (i.e. no assimilation calculated).

	.. option:: MaskCoastsGT25pctOcean
		
		Should we use a looser definition of ocean, including grid cells with at least 25% ocean within the definition? "True" or "False". If oceans are masked, setting this to "True" eliminates many coastal cells which can have problematic satellite retrievals for some products.

	.. option:: Mask60NScaleFactor
		
		Should scaling factors above 60 N always be set to 1? 

	.. option:: Mask60SScaleFactor
		
		Should scaling factors below 60 S always be set to 1? 

	.. option:: 2D_settings

		Settings specific to the case where the user set ``sf_dim`` to be "2D".

		.. option:: speedyCorrelationApprox
			
			See above entry for ``correlatedInitialScalings``. 

		.. option:: additional_init_perturbation_from_emis

			Sometimes users would like to include additional initial emissions perturbations to regions of interest, such as areas with higher emissions; such perturbations allow more granularity in the posterior solution in these locations. The setting is supplied as a dictionary, containing the following values:

			* ``do_add_pert``: ``True`` or ``False``, should additional initial perturbations be applied for this emission?
			* ``file``: a dictionary listing two important pieces of information:
				* ``file``: path to an input NetCDF file used to calculate additional perturbations
				* ``variable``: variable within the NetCDF input file used to calculate perturbations
			* ``max_pert``: the maximum perturbation applied to the initiaon scaling factors (applied additively). For example, a value of 0.5 means that at most a given grid cell will deviate by 0.5 from the default randomization.
			* ``saturation``: all values within the NetCDF variable supplied above greater than this value will map to ``max_pert``. For example, if ``saturation`` is ``1e-9`` then any value greater than ``1e-9`` maps to ``max_pert``, while ``0.5e-9`` would map to half of ``max_pert`` and so on.

			The idea here is that an emissions file can be supplied, and areas with greater emissions will be perturbed by a uniform random variable with a maximum range of :math:`[-\text{max_pert}, \text{max_pert}]`; in areas with lower emissions, the random variable will have a range given by the emission value over the saturation times max_pert.

			An example entry for methane, where initial perturbations are modified sich that areas with greater emissions have greater perturbations, is given below:
			:: 

				"additional_init_perturbation_from_emis" : {
					"do_add_pert":"True",
					"file":{
						"file" : "/n/holyscratch01/jacob_lab/dpendergrass/GC-LETKF/input_data/agg_control_CH4_HEMCO_diagnostics_gmd_paper.nc",
						"variable" : "EmisCH4_Total"
					},
					"max_pert":"0.5",
					"saturation":"0.5e-9"
				},

	.. option:: 1D_settings

		Settings specific to the case where the user set ``sf_dim`` to be "1D". None implemented for now, so leave as ``{}``

.. option:: MinimumScalingFactorAllowed
	
	What is the minimum scaling factor allowed? A dictionary with keys matching ``CONTROL_VECTOR_EMIS`` and with float values for each entry. Set to "nan" if no minimum scaling factor is enforced.

.. option:: MaximumScalingFactorAllowed
	
	As above, but for the maximum scaling factors allowed.

.. option:: InflateScalingsToXOfPreviousStandardDeviation
	
	CHEEREIO includes support for inflating posterior scaling factor standard deviations to a certain percentage of the initial standard deviation. A dictionary with keys matching ``CONTROL_VECTOR_EMIS`` and with float values for each entry, where "0.3" corresponds with inflating to 30% of the initial standard deviation (the recommended value). Set to "nan" to ignore. 

.. option:: MaximumScaleFactorRelativeChangePerAssimilationPeriod
	
	The maximum relative change per assimilation period allowed for scaling factors. For example, a value "0.5" means that no more than a 50% change is allowed for a given scaling factor in a given assimilation period. A dictionary with keys matching ``CONTROL_VECTOR_EMIS`` and with float values for each entry, where "nan" ignores this setting.

.. _LETKF settings:

LETKF settings
~~~~~~~~~~~~~

.. option:: REGULARIZING_FACTOR_GAMMA
	
	A dictionary of regularization factors, with a key corresponding with each key in ``OBSERVED_SPECIES``, which inflates observed error covariance by a factor of :math:`1/\gamma`.

.. option:: USE_DIFFERENT_GAMMA_FOR_BURN_IN

	``True`` or ``False``, should we use a different gamma value during the burn-in period. Supplied as a dictionary with a key corresponding with each key in ``OBSERVED_SPECIES`` and values of ``True`` or ``False``. Users might want to activate this setting to accelerate assimilation during burn-in. For more information on the burn-in period, see :ref:`Burn in period`.

.. option:: GAMMA_FOR_BURN_IN

	A dictionary of regularization factors applied during burn-in only, with a key corresponding with each key in ``OBSERVED_SPECIES``, which inflates observed error covariance by a factor of :math:`1/\gamma`. These can differ from the gamma used during assimilation, but are only applied if the corresponding entry in ``USE_DIFFERENT_GAMMA_FOR_BURN_IN`` is set to ``True``.

.. option:: OBS_ERROR
	
	An dictionary of error information, with a key corresponding with each key in ``OBSERVED_SPECIES`` and a float value. The value is interpreted as one of three categories: "relative", "absolute", or "product". This information represents uncertainty in observations. If error is relative, it is given as a decimal (0.1 means 10% relative error). If error is absolute, it is given as the same units as the observations are in (CHEEREIO will square these values for the covariance matrix). If error is "product," then CHEEREIO uses the error from the observation product. In the product case, the number recorded under OBS_ERROR will not be used. For clarity, only diagonal observational covariance matrices are supported at this time.

.. option:: OBS_ERROR_TYPE
	
	A dictionary of error types, with values given as strings reading "relative",  "absolute", or "product", and with keys corresponding to each key in ``OBSERVED_SPECIES``. This tells CHEEREIO how to interpret the error data types, as described above.

.. option:: OBS_ERROR_SELF_CORRELATION
	
	A dictionary of correlations between errors in data samples, with a key corresponding with each key in ``OBSERVED_SPECIES`` and a float value. This value is used to reduce error if the user would like to aggregate multiple observations together onto the GEOS-Chem grid ("super-observations"). More on this below in the ``AV_TO_GC_GRID`` entry. 

.. option:: MIN_OBS_ERROR
	
	 A dictionary of minimum possible errors, with a key corresponding with each key in ``OBSERVED_SPECIES`` and a float value. If the user would like to aggregate multiple observations together onto the GEOS-Chem grid ("super-observations"), this value gives the minimum possible error allowable upon error reduction. More on this below in the ``AV_TO_GC_GRID`` entry.

.. option:: OTHER_OBS_ERROR_PARAMETERS
	
	A dictionary of dictionaries, with a key corresponding with each key in ``OBSERVED_SPECIES`` and a value that itself is a dictionary with additional settings and their values. At this time, the only setting that is applied using this entry is called ``transport_error``, which is used to account for perfectly correlated model transport errors when the user aggregates multiple observations together onto the GEOS-Chem grid ("super-observations"). More information on this in the the ``AV_TO_GC_GRID`` entry. Below is valid syntax for this setting:
	::

		"OTHER_OBS_ERROR_PARAMETERS":{
			"CH4_TROPOMI":{
				"transport_error":"6.1"
			}
		},

.. option:: AV_TO_GC_GRID
	
	"True" or "False", should observations be averaged to the GEOS-Chem grid? Supplied as a dictionary with a key corresponding with each key in ``OBSERVED_SPECIES`` and values of ``True`` or ``False``. If "False", the above three entries and the below entry are all ignored for this observation. See the :ref:`New superobservation` entry for more information.

.. option:: SUPER_OBSERVATION_FUNCTION
	
	A dictionary with a key corresponding with each key in ``OBSERVED_SPECIES``, and a value corresponding to one of the super observation error reduction functions listed in the :ref:`New superobservation` entry. These include ``default``, ``sqrt``, or ``constant``. Users can add new superobservation functions within the ``produceSuperObservationFunction`` closure in the ``observation_operators.py`` file and activate them from this entry; see the :ref:`New superobservation` entry for more information. 

.. option:: SUPER_OBS_BEFORE_Hx
	
	For some observation operators (e.g. TROPOMI), users have the option to calculate superobservations *before* applying the observation operator. This leads to substantial computational savings. 


.. option:: INFLATION_FACTOR
	
	:math:`\rho-1` from Hunt et. al. (2007). A small number (start with something between 0 and 0.1 and slowly increase according to testing) that inflates the ensemble range. In ensemble Kalman filters, uncertainty usually decreases too quickly and must manually be reinflated.

.. option:: ASSIM_TIME
	
	Length in hours of assimilation window. The assimilation window refers to the period in which GEOS-Chem is run and observations are accumulated; the data assimilation update is calculated in one go within this window. The data assimilation literature contains extensive discussion of this concept.

.. option:: MAXNUMOBS
	
	Maximum number of observations used in a column assimilation calculation. If the number of observations available is greater than this value, then CHEEREIO will randomly throw out observations until only ``MAXNUMOBS`` remain.

.. option:: MINNUMOBS
	
	Minimum number of observations for a column assimilation calculation to be performed. If the number of observations is below this number, no assimilation is calculated and the posterior is set to the prior.

.. option:: LOCALIZATION_RADIUS_km
	
	When updating a column, CHEEREIO only considers data and observations within this radius (in kilometers).

.. option:: smooth_localization_with_gaspari_cohn

	``True`` or ``False``, weight observations such that more distant ones count less and smoothly transition to a weight of zero at or exceeding the ``LOCALIZATION_RADIUS_km``. This leads to smoother behavior; weighting is supplied by the Gaspari-Cohn function, which is a piecewise function resembling a Gaussian with compact support.

.. option:: AveragePriorAndPosterior
	
	"True" or "False", should the posterior be set to a weighted average of the prior and the posterior calculated in the LETKF algorithm? If set to true, the prior weight in the average is given by ``PriorWeightinPriorPosteriorAverage`` in the next setting.

.. option:: PriorWeightinPriorPosteriorAverage
	
	The prior weight if averaging with the posterior from the LETKF. A value between 0 and 1.

.. option:: AverageScaleFactorPosteriorWithPrior
	
	"True" or "False", should the posterior scaling factors be set to a weighted average of the prior (i.e. 1) and the posterior calculated in the LETKF algorithm? If set to true, the prior weight in the average is given by ``PriorWeightinSFAverage`` in the next setting.

	.. attention::

		This setting can lead your ensemble to ignore observations, and often comes up in debugging emails to you. If your simulation is not updating, try setting to "False."


.. option:: PriorWeightinSFAverage
	
	The prior weight if averaging scaling factors with the posterior from the LETKF. A value between 0 and 1.

.. option:: Activate_Relaxation_To_Prior_Spread

	"True" or "False", should we perform Relaxation to Prior Spread (RTPS) inflation in the LETKF assimilation? RTPS is a common form of ensemble inflation (can be done in lieu of ``INFLATION_FACTOR``) where the ensemble is reinflated after assimilation such that the standard deviation of the inflated analysis perturbation matrix :math:`\mathbf{X}^a_\text{infl}` equals the weighted average of the standard deviation (:math:`\sigma^b`) of the background perturbation matrix :math:`\mathbf{X}^b` and the standard deviation (:math:`\sigma^a`) of the analysis posterior ensemble standard deviation :math:`\mathbf{X}^a`, as in the following equation: :math:`\mathbf{X}^a_\text{infl}=\frac{\alpha\sigma^b+(1-\alpha)\sigma^a}{\sigma^a}\mathbf{X}^a`. Here, :math:`\alpha` is a weighting parameter ranging from 0 to 1. 

.. option:: RTPS_parameter

	If ``Activate_Relaxation_To_Prior_Spread`` is True, the parameter :math:`\alpha` from the above.

.. option:: species_not_in_statevec_to_RTPS

	If ``Activate_Relaxation_To_Prior_Spread`` is True, also inflate species not in the state vector? Normally RTPS is only applied to species in the statevector but users can optionally apply RTPS to other species simulated by GEOS-Chem and saved into the restart.

.. _Run in place settings:

Run-in-place settings
~~~~~~~~~~~~~

.. option:: DO_RUN_IN_PLACE
	
	``True`` or ``False``, should we do a run-in-place simulation. See :ref:`Run in place` for more information.

.. option:: rip_update_time
	
	In hours, the run-in-place assimilation window. In run-in-place simulations, the ``ASSIM_TIME`` variable is now interpreted as the run-in-place observation window. See :ref:`Run in place` for more information.

.. option:: DIFFERENT_RUN_IN_PLACE_FOR_BURN_IN
	
	``True`` or ``False``, should we use a different assimilation window during the burn in period for a run-in-place simulation. See :ref:`Run in place` for more information.

.. option:: rip_burnin_update_time
	
	In hours, the run-in-place assimilation window during the burn in period only.

.. option:: DO_VARON_RERUN
	
	``True`` or ``False``, should we do a rerun simulation. See :ref:`Rerun` for more information.

.. option:: APPROXIMATE_VARON_RERUN
	
	``True`` or ``False``, if using a rerun simulation with ``number_of_windows_to_rerun`` set to a value larger than 1, should we use linear regression to extrapolate GEOS-Chem results? This is an advanced setting. See :ref:`Rerun` for more information.

.. option:: species_to_approximate_for_rerun
	
	A list of species concentrations. If ``APPROXIMATE_VARON_RERUN`` is ``True``, which species concentrations should we approximate with linear regression? This is an advanced setting. See :ref:`Rerun` for more information.

.. option:: number_of_windows_to_rerun
	
	If ``DO_VARON_RERUN`` is set to ``True``, how many assimilation windows should we rerun? For most users, this value should be 1. See :ref:`Rerun` for more information.

.. _postprocessing settings:

Postprocessing settings
~~~~~~~~~~~~~

.. option:: animation_fps_scalingfactor
	
	Frames per second for movies of scaling factors and emissions made by the postprocessing workflow.

.. option:: animation_fps_concentrations
	
	Frames per second for movies of concentrations made by the postprocessing workflow.

.. option:: omit_diff_cells_with_fewer_than_n_observations
	
	Some postprocessing plots show the difference between assimilation output and observations; we can omit grid cells with very few observations to prevent noisy results. 

.. option:: hemco_diags_to_process
	
	An array of entries from HEMCO Diagnostics that you would like processed into movies. This is usually emissions totals. See below for an example entry:
	:: 

		"hemco_diags_to_process" : [
			"EmisCH4_Total"
		],

.. option:: useLogScaleForEmissionsMaps
	
	"True" or "False", use log color scale to plot emissions from HEMCO.

.. option:: min_emis_value_to_plot
	
	Minimum emission value to plot; useful if using a log scale.

.. option:: min_emis_std_value_to_plot
	
	Minimum ensemble standard deviation emission value to plot; useful if using a log scale.

.. option:: OBSERVATION_UNITS
	
	A dictionary with keys from ``OBSERVED_SPECIES`` and values representing the units will be plotted. This is governed by how the observation operator is defined. 

.. option:: EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y
	
	A dictionary with keys from ``OBSERVED_SPECIES`` and values listing entries within the ObsData class (created by observation operators) which should be saved to postprocessing. Available obsdata fields will differ by observation operators and are listed on the :ref:`Existing obs tools` page in the section corresponding with the observation operator you are using. For example, TROPOMI CH4 supports saving albedo information for diagnostic plots, while ObsPack CH4 has additional identification metrics. A user wishing to save three kinds of albedo from TROPOMI and obspack ids would supply the following:
	::

		"EXTRA_OBSDATA_FIELDS_TO_SAVE_TO_BIG_Y":{
			"CH4_TROPOMI" : [
				"albedo_swir",
				"albedo_nir",
				"blended_albedo"
			],
			"CH4_OBSPACK" : [
				"obspack_id",
				"site_code"
			]
		},


.. option:: EXTRA_OBSDATA_FIELDS_TO_REGRID_AND_PLOT
	
	A ditionary with keys from ``OBSERVED_SPECIES`` and values listing entries within the ObsData class (created by observation operators) that should be regridded and plotted in the postprocessing step. This usually is identical with the above entry unless users would like to save but not plot a given field (this is common for ObsPack data, where obspack_id may be saved but does not make sense to plot). The user in the above case might supply the following:
	::

		"EXTRA_OBSDATA_FIELDS_TO_REGRID_AND_PLOT":{
			"CH4_TROPOMI" : [
				"albedo_swir",
				"albedo_nir",
				"blended_albedo"
			],
			"CH4_OBSPACK" : []
		},

	This is because the ObsPack ID information is not useful to plot.

	If you are writing a new observation operator and would like to ensure that postprocessing routines will plot an important observational quantity (e.g. albedo), see the :ref:`New field in postprocessing` page.

.. option:: extra_plot_field_units
	
	If users are plotting additional ObsData fields supplied in the above entry, then this entry supplies the units that are the plots. Supplied as a dictionary where keys are all the ObsData fields listed above and values are units. An example for albedo is supplied below:
	::

		"extra_plot_field_units":{
			"albedo_swir":"Albedo",
			"albedo_nir":"Albedo",
			"blended_albedo":"Albedo"
		},


.. option:: OBSERVERS_TO_PLOT_AS_POINTS
	
	A dictionary listing observers which the user would like plotted as points, rather than a grid. Keys can include any keys from ``OBSERVED_SPECIES`` and values are the the field listing individual observer locations. If a value is not listed, then it is assumed to be plotted as a field (default behavior). For example, in the example simulations above, we have  both tropomi and obspack observations (keyed as ``CH4_TROPOMI`` and ``CH4_OBSPACK`` respectively). This user might like to plot individual obspack sites but leave tropomi data to be plotted as a 2D field. Obspack sites are separated by ``site_code``, so the user would supply this entry as follows:
	::

		"OBSERVERS_TO_PLOT_AS_POINTS":{
			"CH4_OBSPACK":"site_code"
		},


.. option:: scalefactor_plot_freq
	
	The CHEEREIO postprocessing routine will save out maps of scale at this temporal resolution: either "all" for save out an image for every assimilation window, or "monthly" to save out monthly means.

.. _Extensions:

Extensions
~~~~~~~~~~~~~

Additional CHEEREIO settings, usually for specific observation types, can be loaded in through extensions. Extensions in CHEEREIO are extra JSON files that store additional settings in order to prevent clutter in the ``ens_config.json`` file. Extensions can easily be added by saving a file with name ``NAME_extension.json`` within the ``extensions/`` folder. To load in the settings within the ``NAME_extension.json`` file, add the key NAME to the "Extensions" dictionary in ``ens_config.json`` with value "True". Below is an example where we load in the settings in the ``TROPOMI_CH4_extension.json`` and ``CH4_extension.json`` files. 
::

	"Extensions": {
		"TROPOMI_CH4":"True",
		"CH4":"True"
	}

Below we list the settings that you can set with extensions.

.. option:: TROPOMI_CH4 extension
	
	Specialized TROPOMI CH4 settings.

	.. option:: WHICH_TROPOMI_PRODUCT

		Which TROPOMI product should we use? "DEFAULT" for the TROPOMI operational product, 'ACMG' for the ACMG TROPOMI product, and 'BLENDED' for Belasus et al., 2023. 

	.. option:: TROPOMI_CH4_FILTERS
	
		Apply specialized filters for TROPOMI methane? Set to "True" if doing a TROPOMI methane inversion, otherwise set to "False".
	
	.. option:: TROPOMI_CH4_filter_blended_albedo
	
		Filter out TROPOMI methane observations with a blended albedo above this value. Set to "nan" to ignore.

	.. option:: TROPOMI_CH4_filter_swir_albedo_low
	
		Filter out TROPOMI methane observations with a SWIR albedo below this value. Set to "nan" to ignore.

	.. option:: TROPOMI_CH4_filter_swir_albedo_high
	
		Filter out TROPOMI methane observations with a SWIR albedo above this value. Set to "nan" to ignore.

	.. option:: TROPOMI_CH4_filter_winter_lat
	
		Filter out TROPOMI methane observations beyond this latitude in the winter hemisphere. Set to "nan" to ignore.

	.. option:: TROPOMI_CH4_filter_roughness
	
		Filter out TROPOMI methane observations with a surface roughness above this value. Set to "nan" to ignore.

	.. option:: TROPOMI_CH4_filter_swir_aot
	
		Filter out TROPOMI methane observations with a SWIR AOT above this value. Set to "nan" to ignore.

.. option:: CH4 extension
	
	Specialized settings for methane simulations.

	.. option:: USE_CUSTOM_CH4_OH_ENTRY
	
		Should we overwrite the OH field setting in ``HEMCO_Config.rc`` for the specialty CH4 simulation with a user-specified entry (given below)? "True" or "False".

	.. option:: CUSTOM_CH4_OH_ENTRY
	
		If USE_CUSTOM_CH4_OH_ENTRY is True, then we overwrite the OH field line in ``HEMCO_Config.rc`` with this entry. Note that backslashes need to be escaped with a front slash. Here is an example entry: ``* GLOBAL_OH  $ROOT\/OH\/v2014-09\/v5-07-08\/OH_3Dglobal.geos5.47L.4x5.nc OH           1985\/1-12\/1\/0 C xyz kg\/m3 * - 1 1`` 

.. option:: TROPOMI_CO extension
	
	Specialized TROPOMI CO settings.

	.. option:: WHICH_TROPOMI_PRODUCT

		Which TROPOMI product should we use? Currently, only "DEFAULT" is supported. 

	.. option:: TROPOMI_CO_FILTERS
	
		Apply specialized filters for TROPOMI CO? Set to "True" if doing a TROPOMI CO inversion, otherwise set to "False".
	
	.. option:: TROPOMI_CO_filter_blended_albedo
	
		Filter out TROPOMI CO observations with a blended albedo above this value. Set to "nan" to ignore.

	.. option:: TROPOMI_CO_filter_swir_albedo_low
	
		Filter out TROPOMI CO observations with a SWIR albedo below this value. Set to "nan" to ignore.

	.. option:: TROPOMI_CO_filter_swir_albedo_high
	
		Filter out TROPOMI CO observations with a SWIR albedo above this value. Set to "nan" to ignore.

	.. option:: TROPOMI_CO_filter_winter_lat
	
		Filter out TROPOMI CO observations beyond this latitude in the winter hemisphere. Set to "nan" to ignore.

	.. option:: TROPOMI_CO_filter_roughness
	
		Filter out TROPOMI CO observations with a surface roughness above this value. Set to "nan" to ignore.

	.. option:: TROPOMI_CO_filter_swir_aot
	
		Filter out TROPOMI CO observations with a SWIR AOT above this value. Set to "nan" to ignore.

.. option:: TCCON_CO extension
	
	Specialized TCCON CO settings. None currently implemented, but reserved for future use.

	.. option:: TCCON_CO_FILTERS
	
		Apply specialized filters for TROPOMI CO? Set to "True" if doing a TROPOMI CO inversion, otherwise set to "False". Leave as "False" as none currently implemented.

.. option:: N2O extension
	
	Specialized N2O settings for a custom GC simulation, as yet unreleased. Other users can ignore.

	.. option:: INTERPRET_GC_CH4_AS_N2O
	
		In the custom GC simulation for N2O I am currently using, we save out N2O as CH4 (it is a modification of the CH4 simulation). If this setting is set to True, it overwrites this with "N2O" in internal CHEEREIO logic.

.. option:: TCCON_N2O extension
	
	Specialized TCCON N2O settings.

	.. option:: PT700_N2O_CORRECTION
	
		Apply Josh Laughner's temperature-dependent bias correction to N2O for GGG2020.0? Set to "True" if so, otherwise set to "False". Will no longer be necessary with GGG2020.1.

	.. option:: TCCON_N2O_FILTERS
	
		Apply specialized filters for TCCON N2O? Set to "True" if doing a TCCON N2O inversion, otherwise set to "False". Leave as "False" as none currently implemented.

.. option:: Obspack_N2O extension 

	Specialized ObsPack N2O settings

	.. option:: Max_Obspack_N2O_Error
	
		Maximum allowable difference (in ppb) between observations and GEOS-Chem. If it exceeds this threshold, discard. Set to nan to ignore.

.. option:: OMI_NO2 extension
	
	Specialized settings for the OMI NO2 operator.

	.. option:: OMI_NO2_FILTERS
	
		Apply specialized filters for OMI NO2? Set to "True" if doing an OMI NO2 inversion, otherwise set to "False".

	.. option:: OMI_NO2_filter_sza
	
		Filter out OMI NO2 observations with a solar zenith angle above this value. Set to "nan" to ignore.

	.. option:: OMI_NO2_filter_cloud_radiance_frac
	
		Filter out OMI NO2 observations with a cloud radiance fraction above this value. Set to "nan" to ignore.
	
	.. option:: OMI_NO2_filter_surface_albedo
	
		Filter out OMI NO2 observations with a surface albedo above this value. Set to "nan" to ignore.

