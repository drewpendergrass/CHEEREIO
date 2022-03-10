.. _Configuration:

The ensemble configuration file
==========

Because CHEEREIO requires many of the ensemble settings to be globally available to many different types of scripts and programs, written in different languages and stored in in different directories, it expects a very particular kind of installation process and for user settings to be made available in a strict format (so that they can be read by scripts in multiple languages). The overall install process is not too different from the way GEOS-Chem run directories are created as of version 13.0.0, but it does require some care on the part of the user. This page explains how to customize the ensemble to meet your needs. This isn't merely a technical process: the user makes important scientific assumptions in this step!

All ensemble configuration is set by the ``ens_config.json`` file in the main CHEEREIO code directory. This file contains all the information that the user would normally input in the process of making a GEOS-Chem run directory, in addition to other settings like cluster configuration, global assimilation variables (like the localization radius), the composition of the state and control vectors, links to observations, and details about which emissions can be assimilated. This file is in JSON format, which is a format that is human readable but enforces a strict syntax. 

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
		"NOx":["NO","NO2"]
	},

However, CHEEREIO does expect all values in ``ens_config.json`` as strings. If a setting is given by a number, wrap that number in quotation marks.  

**Important**: There is one subtlety in this particular configuration file: colons in the ``WallTime`` and ``SpinupWallTime`` entries (covered later) must be escaped with two backslashes (``\\``). The first backslash escapes the second backslash in JSON; the second backslash escapes the colon in the SED unix utility which is used in the CHEEREIO installation process. For example, to allow a one day eight hour simulation, you would write ``"WallTime" : "1-08\\:00",``.

More details on the JSON format are available on the JSON `website <https://www.json.org>`__. When in doubt, follow the template ``ens_config.json`` file!

A line-by-line guide to ensemble configuration
-------------

The rest of this section will cover the various parts of the ``ens_config.json`` file and the settings they control.


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
* CH4_HEMCO_ROOT: If the subsequent option, "USE_CHEEREIO_TEMPLATE_CH4_HEMCO_Config", is set to "True", then this is the root folder where emissions and other input files for the methane specialty simulation are located. In this case, a special CHEEREIO ``HEMCO_Config.rc`` template from the ``templates/`` folder in the code directory is used. *Note: this option is functional but currently causes GEOS-Chem crashes with an unknown cause (DP, 2022/03/09).*
* RESTART_FILE: Full path to the restart file for the simulation.
* BC_FILES: Full path to the boundary condition files for the simulation if a nested grid (empty string otherwise).
* sim_name: Simulation type. Valid options are "fullchem", "aerosol", "CH4", "CO2", "Hg", "POPs", "tagCH4", "tagCO", "tagO3", and "TransportTracers".
* chemgrid: Options are "trop+strat" and "trop_only".
* sim_extra_option: Options are "none", "benchmark", "complexSOA", "complexSOA_SVPOA", "marinePOA", "aciduptake", "TOMAS15", "TOMAS40", "APM", "RRTMG", "BaP", "PHE", and "PYR". Depending on the simulation type only some will be available. Consult the GEOS-Chem documation for more information.
* DO_SPINUP: Would you like CHEEREIO to set up a spinup directory for you? "true" or "false". The ensemble will automatically start from the end restart file produced by this run. Note this option is for the standard GEOS-Chem spinup (run once for the whole ensemble).
* SPINUP_START: Start date for spinup (YYYYMMDD). Empty string if no spinup.
* SPINUP_END: End date for spinup (YYYYMMDD).
* DO_ENS_SPINUP: Do you want to use a separate job array to spin up your GEOS-Chem ensemble with randomized scaling factors applied to each ensemble member? "true" or "false". If set to "true", shell scripts entitled ``run_ensemble_spinup_simulations.sh`` and ``run_ensspin.sh`` are installed in the ``ensemble_runs/`` folder. The user should then execute ``run_ensspin.sh`` to spin up the ensemble and create variability between ensemble members before executing ``run_ens.sh`` in the normal run procedure. 
* ENS_SPINUP_FROM_BC_RESTART: It is possible to start the ensemble spinup procedure using a boundary condition file, rather than a traditional restart file. Set to "true" if using a BC file, and "false" if using a normal restart file to start the ensemble spinup.
* ENS_SPINUP_START: Start date for ensemble spinup run (YYYYMMDD).
* ENS_SPINUP_END: End date for ensemble spinup run (YYYYMMDD).
* START_DATE: Start date for main ensemble data assimilation run (YYYYMMDD).
* ASSIM_START_DATE: Date where assimilation begins (YYYYMMDD). After GEOS-Chem version 13.4 this option can be used in lieu of ``DO_ENS_SPINUP``; just set this date to be sufficiently far away from ``START_DATE``. Prior to version 13.4, it is buggy to run GEOS-Chem for a non-standard length of time (e.g. 4 months and a week) which is usually desired for the ensemble spinup. For these versions, the separate ensemble spinup script installed by ``DO_ENS_SPINUP`` is a good work-around.  
* END_DATE: End date for ensemble run (YYYYMMDD).
* nEnsemble: Number of ensemble members. 32 is usually a good number. This number of run directories will be created in the ``ensemble_runs`` folder and will be run simultaneously.
* pPERT: Setting for initial emissions scaling factor creation, where the number provided :math:`p` is used to generate random scaling factors from the distribution :math:`p^u,\\ u\\sim\\mathcal{U}(-1,1)`, meaning that u is a uniform random variable ranging from -1 to 1. For example, if ``pPERT`` is "4" then scalings will range from 0.25 to 4, centered on 1.
* SIMULATE_NATURE: *Deprecated: will be removed before official release (DP, 2022/03/09)*. End users should leave this set to "false", as this was used for testing in early CHEEREIO development.  

Cluster settings
~~~~~~~~~~~~~

The next section of the ``ens_config.json`` file controls settings that will be used when submitting jobs to the scheduler. These settings overwrite the template batch submission scripts included with CHEEREIO.

* NumCores: Number of cores used in each of the ensemble runs. CHEEREIO also will use these cores to parallelize assimilation computation columnwise.
* Partition: Partition of your cluster you are submitting to. At Harvard, ``huce_intel`` is a good choice.
* Memory: Memory in megabytes used by each ensemble member. CHEEREIO is quite memory intensive because it loads in restarts and history files for many ensemble members in addition to observations, so expect to use more than in standard GEOS-Chem runs.
* WallTime: Time allowed for the overall assimilation process (runs and assimilation) to occur in format D-HH\\:MM. Assimilation adds substantial overhead so expect it to be slow.
* EnsSpinupWallTime: Time allowed for the ensemble spinup process (no assimilation, just running all ensemble members from ``ENS_SPINUP_START`` through ``ENS_SPINUP_END`` with scaling factors applied) in format D-HH\\:MM. If not using, you can just leave as an empty string.
* SpinupWallTime: Wall time for the spinup simulation, if you're using one. Empty string otherwise.
* CondaEnv: The name of the Conda environment with all of the CHEEREIO packages installed. It is strongly recommended that you install an environment using the YAML file that ships with CHEEREIO in the ``environments/`` folder.
* MaxPar: Maximum number of cores to use while assimilating columns in parallel using CHEEREIO, maxing out at ``NumCores``. Setting this number smaller than NumCores saves on memory but adds to the assimilation time. 

Species in state/control/observation vectors
~~~~~~~~~~~~~

* STATE_VECTOR_CONC: Species from the restart files to be included in the state vector. It is generally recommended to include a fairly wide range of species that might affect the species you are mainly interested in, but not so large a range that you end up analyzing noise. Given as an array. This is an example for NO\ :sub:`x`\ data assimilation: 
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

* CONTROL_VECTOR_CONC: A subset of the state vector concentration species that will be updated by assimilation. Although an update for all members of the state vector will be calculated, only the species listed in this array will have that update saved. This allows a wide range of species to be considered in the update calculation process but only a smaller, more tightly coupled subset of species to actually be changed and passed to GEOS-Chem. The goal is to tamp down on noise. 
* CONTROL_VECTOR_EMIS: A dictionary linking a label for emissions scalings to the species emitted. For example, you could write ``"CH4_WET" : "CH4"`` to reference wetland methane emissions. CHEEREIO automatically will update ``HEMCO_Config.rc`` accordingly, but cannot distinguish between different emissions of the same species on its own; the user has to manually edit ``HEMCO_Config.rc`` to correct this if distinguishing between different sources of the same species. More on this in :ref:`Template`. You can also use one label to link to emissions of multiple species, meaning that all these emissions will be controlled by one scaling factor file, such as ``"NOx":["NO","NO2"]`` to indicate that NO and NO\ :sub:`2`\ are controlled by one scaling factor file. Here are a couple of examples:
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

* HISTORY_collections_to_customize: A list of collections under HISTORY.rc that CHEEREIO will customize with user-specified frequency and duration settings. Here is a typical example:
::

	"HISTORY_collections_to_customize" : [
		"SpeciesConc",
		"LevelEdgeDiags",
		"StateMet"
	],

* HISTORY_freq: Frequency of data saved within history output files listed within collections in ``HISTORY_collections_to_customize``. For more information on history frequencies, see the GEOS-Chem manual.
* HISTORY_dur: As in ``HISTORY_freq``, but for duration.
* SPINUP_HISTORY_freq: Frequency of history output files saved during ensemble spinup (i.e. when executing ``run_ensspin.sh``).
* SPINUP_HISTORY_dur: As in ``SPINUP_HISTORY_freq``, but for duration.
* SaveLevelEdgeDiags: Should the LevelEdgeDiags collection be turned on? "True" or "False". This is mandatory for assimilating TROPOMI data.
* SaveStateMet: Should the StateMet collection be turned on? "True" or "False". This is mandatory for assimilating TROPOMI NO\ :sub:`2`\ .
* SaveArea: Should grid cell areas be used in the assimilation process? "True" or "False". This is mandatory for assimilating TROPOMI NO\ :sub:`2`\ .
* HistorySpeciesConcToSave: A list of species to save in the SpeciesConc collection. At minimum, this should encompass the concentration portion of the state vector. Below is an example: 
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

* HistoryLevelEdgeDiagsToSave: A list of data to save in the LevelEdgeDiags collection. Just ``Met_PEDGE`` is sufficient for assimilation.
* HistoryStateMetToSave: A list of data to save in the StateMet collection. Just ``Met_AD`` is sufficient for assimilation.

Observation settings
~~~~~~~~~~~~~

* OBSERVED_SPECIES: A dictionary linking a label for observations with the species observed. For example, you could write ``"NO2_SATELLITE" : "NO2"`` to reference satellite observations of NO2. Unlike elsewhere, here the order matters. Later in the configuration file, arrays of errors and regularization factors will be associated with these species according to the order they are stored. More in the next section. 
* OBS_4D: *Deprecated: until removed, make an array of "True" of length OBSERVED_SPECIES*.
* OBS_TYPE_TROPOMI: Array of length ``OBSERVED_SPECIES``, with a value of "True" or "False" if the observation is from TROPOMI.
* TROPOMI_dirs: Dictionary linking observed TROPOMI species to the directory containing the observations. Here is an example, along with the corresponding ``OBSERVED_SPECIES`` settings:
::

	"OBSERVED_SPECIES" : {
		"NO2_TROPOMI": "NO2"
	},
	"TROPOMI_dirs" : {
		"NO2" : "/n/holylfs05/LABS/jacob_lab/dpendergrass/tropomi/NO2/2019"
	},

* OMI_dirs: As in TROPOMI_dirs, but for OMI.
* LOW_MEMORY_TROPOMI_AVERAGING_KERNEL_CALC: For TROPOMI observations, should CHEEREIO use the "low memory" algorithm to apply the averaging kernel ("True") or the default fast algorithm ("False"). For highly dense observations like TROPOMI NO\ :sub:`2`\ set to "True", otherwise leave as "False"
* LOW_MEMORY_TROPOMI_AVERAGING_KERNEL_BATCH_SIZE: Batch size for the low memory TROPOMI averaging kernel algorithm. Users should probably leave as default (100000), but the user can increase this value to speed up the algorithm at the cost of memory, or vice versa.
* TROPOMI_CH4_FILTERS: Apply specialized filters for TROPOMI methane? Set to "True" if doing a TROPOMI methane inversion, otherwise set to "False".
* TROPOMI_CH4_filter_blended_albedo: Filter out TROPOMI methane observations with a blended albedo above this value. Set to "nan" to ignore.
* TROPOMI_CH4_filter_swir_albedo_low: Filter out TROPOMI methane observations with a SWIR albedo below this value. Set to "nan" to ignore.
* TROPOMI_CH4_filter_swir_albedo_high: Filter out TROPOMI methane observations with a SWIR albedo above this value. Set to "nan" to ignore.
* TROPOMI_CH4_filter_winter_lat: Filter out TROPOMI methane observations beyond this latitude in the winter hemisphere. Set to "nan" to ignore.
* TROPOMI_CH4_filter_roughness: Filter out TROPOMI methane observations with a surface roughness above this value. Set to "nan" to ignore.
* TROPOMI_CH4_filter_swir_aot: Filter out TROPOMI methane observations with a SWIR AOT above this value. Set to "nan" to ignore.

Scaling factor settings
~~~~~~~~~~~~~

* MaskOceanScaleFactor: Should scaling factors be allowed to vary over the oceans? An array of "True" or "False" of length ``OBSERVED_SPECIES`` where each boolean corresponds with the species in the same index in ``OBSERVED_SPECIES``. If "True", scaling factors for that species over the ocean are always set to 1 across all ensemble members.
* MaskCoastsGT25pctOcean: Should we use a looser definition of ocean, including grid cells with at least 25% ocean within the definition? "True" or "False". If oceans are masked, setting this to "True" eliminates many coastal cells which can have problematic satellite retrievals for some products.
* Mask60NScaleFactor: Should scaling factors above 60 N always be set to 1? An array of "True" or "False" of length ``OBSERVED_SPECIES`` where each boolean corresponds with the species in the same index in ``OBSERVED_SPECIES``.
* Mask60SScaleFactor: Should scaling factors below 60 S always be set to 1? An array of "True" or "False" of length ``OBSERVED_SPECIES`` where each boolean corresponds with the species in the same index in ``OBSERVED_SPECIES``.
* MinimumScalingFactorAllowed: What is the minimum scaling factor allowed? An array of numbers of length ``OBSERVED_SPECIES`` where each number corresponds with the species in the same index in ``OBSERVED_SPECIES``. Set to "nan" if no minimum scaling factor is enforced.
* MaximumScalingFactorAllowed: As above, but for the maximum scaling factors allowed.
* InflateScalingsToXOfPreviousStandardDeviation: Following the work of Kazuyuki Miyazaki, CHEEREIO includes support for inflating posterior scaling factor standard deviations to a certain percentage of the prior standard deviation. An array of numbers of length ``OBSERVED_SPECIES``, where "0.3" corresponds with inflating to 30% of the prior standard deviation (the recommended value). Set to "nan" to ignore. 
* MaximumScaleFactorRelativeChangePerAssimilationPeriod: The maximum relative change per assimilation period allowed for scaling factors. For example, a value "0.5" means that no more than a 50% change is allowed for a given scaling factor in a given assimilation period. An array of numbers of length ``OBSERVED_SPECIES``, where "nan" ignores this setting.

LETKF settings
~~~~~~~~~~~~~

* REGULARIZING_FACTOR_GAMMA: An array of regularization factors, corresponding with ``OBSERVED_SPECIES``, which inflates observed error covariance by a factor of :math:`\gamma`. 
* OBS_COVARIANCE: An array of errors, either relative or absolute, representing uncertainty in observations. If error is relative, it is given as a decimal (0.1 means 10% relative error). If error is absolute, it is given in the square of whatever units the observations are in. Order is the same as ``OBSERVED_SPECIES``. For clarity, only diagonal observational covariance matrices are supported at this time.
* OBS_COVARIANCE_TYPE: An array of error types, given as strings reading either "relative" or "absolute", corresponding to each error value in ``OBS_COVARIANCE``. This tells CHEEREIO how to interpret the covariance data type. 
* INFLATION_FACTOR: :math:`\rho-1` from Hunt et. al. (2007). A small number (start with something between 0 and 0.1 and slowly increase according to testing) that inflates the ensemble range. In ensemble Kalman filters, uncertainty usually decreases too quickly and must manually be reinflated.
* ASSIM_TIME: Length in hours of assimilation window. The assimilation window refers to the period in which GEOS-Chem is run and observations are accumulated; the data assimilation update is calculated in one go within this window. The data assimilation literature contains extensive discussion of this concept.
* MAXNUMOBS: Maximum number of observations used in a column assimilation calculation. If the number of observations available is greater than this value, then CHEEREIO will randomly throw out observations until only ``MAXNUMOBS`` remain.
* MINNUMOBS: Minimum number of observations for a column assimilation calculation to be performed. If the number of observations is below this number, no assimilation is calculated and the posterior is set to the prior.
* AV_TO_GC_GRID: "True" or "False", should observations be averaged to the GEOS-Chem grid? *Note that errors are not currently updated after averaging; this is a top development priority*.
* LOCALIZATION_RADIUS_km: When updating a column, CHEEREIO only considers data and observations within this radius (in kilometers).
* AveragePriorAndPosterior: "True" or "False", should the posterior be set to a weighted average of the prior and the posterior calculated in the LETKF algorithm? If set to true, the prior weight in the average is given by ``PriorWeightinPriorPosteriorAverage`` in the next setting.
* PriorWeightinPriorPosteriorAverage: The prior weight if averaging with the posterior from the LETKF. A value between 0 and 1.

Miscellaneous settings
~~~~~~~~~~~~~

* animation_fps_scalingfactor: Frames per second for movies of scaling factors and emissions made by the postprocessing workflow.
* animation_fps_concentrations: Frames per second for movies of concentrations made by the postprocessing workflow.
* postprocess_save_albedo: Should the postprocessing workflow save out albedo? "True" or "False".
