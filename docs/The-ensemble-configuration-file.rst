.. _Configuration:

Configuring your simulation
==========

All ensemble configuration is set by the ``ens_config.json`` file in the main CHEEREIO code directory. This file contains all the information that the user would normally input in the process of making a GEOS-Chem run directory, in addition to other settings like cluster configuration, global assimilation variables (like the localization radius), the composition of the state and control vectors, links to observations, and details about which emissions can be assimilated. This file is in JSON format, which is a format that is human readable but enforces a strict syntax. 

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

More details on the JSON format are available on the JSON `website <https://www.json.org>`__. When in doubt, follow the conventions in the template ``ens_config.json`` file!

A line-by-line guide to ensemble configuration
-------------

The rest of this section will cover the various parts of the ``ens_config.json`` file and the settings they control.


Basic GEOS-Chem and ensemble settings
~~~~~~~~~~~~~

The first section of the ``ens_config.json`` file (i.e. between the first two comments) mostly controls settings analagous to those set during normal GEOS-Chem run directory creation. However, there are a few unique options in this setting particular to CHEEREIO. We'll consider these one-by-one.

* RES: The resolution of the GEOS-Chem model. Options are available on the `GEOS-Chem website <http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_horizontal_grids>`__ and include 4.0x5.0, 2.0x2.5, 0.5x0.625, 0.25x0.3125 and nested grid settings in format TwoLetterCode_MetCode (e.g. AS_MERRA2, EU_GEOSFP). Custom nested domains are not currently supported by the automated scaling factor creation utility but can be manually added by the user. If there is enough interest I will add more automated support in a later CHEEREIO update.
* met_name: Meteorology (chosen from MERRA2, GEOSFP, or ModelE2.1).
* LEVS: Number of levels (47 or 72).
* NEST: Is this a nested grid simulation? "T" or "F".
* REGION: Two letter region code for nested grid, or empty string ("") if not.
* ASSIM_PATH: **Full path** to the directory where the CHEEREIO repository is installed (e.g. ``/n/home12/drewpendergrass/CHEEREIO``). Directories in the ``ens_config.json`` file **should not have trailing forward slashes.** Again, when in doubt follow the provided templates.
* RUN_NAME: The name of the CHEEREIO ensemble run (will be the name of the folder containing the ensemble, template run directory, temporary files, and so on).
* MY_PATH: Path to the directory where ensembles will be created. A folder with name ``RUN_NAME`` will be created inside.
* DATA_PATH: Path to where external GEOS-Chem data is located. This can be an empty string if GEOS-Chem has already been configured on your machine (it is automatically overwritten).
* CH4_HEMCO_ROOT: If the subsequent option, "USE_CHEEREIO_TEMPLATE_CH4_HEMCO_Config", is set to "True", then this is the root folder where emissions and other input files for the methane specialty simulation are located. In this case, a special CHEEREIO ``HEMCO_Config.rc`` template from the ``templates/`` folder in the code directory is used. *Note: this option is functional but currently causes GEOS-Chem crashes with an unknown cause (DP, 2022/03/09).*
* RESTART_FILE: Full path to the restart file for the simulation. If in the initialization process you selected ``SetupSpinupRun=true``, then this restart file will be used for the classic spin up routine (getting realistic atmospheric conditions for the entire ensemble). Otherwise, this will be the restart file used to initialize all ensemble members.
* BC_FILES: Full path to the boundary condition files for the simulation if you are using a nested grid (empty string otherwise).
* sim_name: Simulation type. Valid options are "fullchem", "aerosol", "CH4", "CO2", "Hg", "POPs", "tagCH4", "tagCO", "tagO3", and "TransportTracers".
* chemgrid: Options are "trop+strat" and "trop_only".
* sim_extra_option: Options are "none", "benchmark", "complexSOA", "complexSOA_SVPOA", "marinePOA", "aciduptake", "TOMAS15", "TOMAS40", "APM", "RRTMG", "BaP", "PHE", and "PYR". Depending on the simulation type only some will be available. Consult the GEOS-Chem documation for more information.
* DO_SPINUP: Would you like CHEEREIO to set up a spinup directory for you? "true" or "false". The ensemble will automatically start from the end restart file produced by this run. Note this option is for the standard GEOS-Chem spinup (run once for the whole ensemble). Note that if this is activated, you have to run the ``setup_ensemble.sh`` utility with the ``SetupSpinupRun`` switch set to ``true``.
* SPINUP_START: Start date for spinup (YYYYMMDD). Empty string if no spinup.
* SPINUP_END: End date for spinup (YYYYMMDD).
* DO_CONTROL_RUN: The control run is a normal GEOS-Chem simulation without any assimilation. The output of this simulation can be compared with the LETKF results in the postprocessing workflow. Set to "true" if using a control run. Note that if this is activated, you have to run the ``setup_ensemble.sh`` utility with the ``SetupControlRun`` switch set to ``true``.
* CONTROL_START: Start date for the control run (YYYYMMDD).
* CONTROL_END: End date for the control run (YYYYMMDD).
* DO_ENS_SPINUP: Do you want to use a separate job array to spin up your GEOS-Chem ensemble with randomized scaling factors applied to each ensemble member? "true" or "false". If set to "true", shell scripts entitled ``run_ensemble_spinup_simulations.sh`` and ``run_ensspin.sh`` are installed in the ``ensemble_runs/`` folder. The user should then execute ``run_ensspin.sh`` to spin up the ensemble and create variability between ensemble members before executing ``run_ens.sh`` in the normal run procedure. For more information on the ensemble spinup process, see :ref:`Run Ensemble Spinup Simulations`.
* ENS_SPINUP_FROM_BC_RESTART: It is possible to start the ensemble spinup procedure using a boundary condition file, rather than a traditional restart file. Set to "true" if using a BC file, and "false" if using a normal restart file to start the ensemble spinup.
* ENS_SPINUP_START: Start date for ensemble spinup run (YYYYMMDD).
* ENS_SPINUP_END: End date for ensemble spinup run (YYYYMMDD).
* START_DATE: Start date for main ensemble data assimilation run (YYYYMMDD).
* ASSIM_START_DATE: Date where assimilation begins (YYYYMMDD). This option allows you to run the first assimilation period for an extra long time (although the assimilation window remains the same), effectively providing an ensemble-wide spinup. For more information on this ensemble spinup option, see :ref:`Run Ensemble Spinup Simulations`. If you have set ``DO_ENS_SPINUP`` to ``true``, then you should set this date to be one assimilation window later than ``START_DATE``.
* SIMPLE_SCALE_FOR_FIRST_ASSIM_PERIOD" : At the end of the first assimilation period, rather than doing the full LETKF calculation, CHEEREIO can scale the ensemble mean so that it matches the observational mean. This is done because if the model is biased relative to observations the LETKF will perform suboptimal updates. Set to "true" to do this scaling (recommended) or "false" to do the usual LETKF calculation. For more information, see :ref:`Simple scale`.
* END_DATE: End date for ensemble run (YYYYMMDD).
* AMPLIFY_ENSEMBLE_SPREAD_FOR_FIRST_ASSIM_PERIOD: At the end of the ensemble spinup period, the spread in ensemble members may still not be great enough. If this option is set to "true", then CHEEREIO will multiply the standard deviation of the ensemble after ensemble spinup is complete by the factor given in ``SPREAD_AMPLIFICATION_FACTOR``. This spread amplification is done after the first assimilation period, so it will work with either spinup method. For more information, see :ref:`Spread amplification`.
* SPREAD_AMPLIFICATION_FACTOR: If ``AMPLIFY_ENSEMBLE_SPREAD_FOR_FIRST_ASSIM_PERIOD`` is set to "true", then this is the factor with which CHEEREIO will multiply the ensemble standard deviation at the end of the ensemble spinup period. For more information on the burn-in period, see :ref:`Burn in period`.
* SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD" : Should CHEEREIO do a burn-in period? "true" or "false." A burn-in period is a time period where full LETKF assimilation is being applied, but the results will be discarded from final analysis. The idea of a burn in period is to allow CHEEREIO's emissions to "catch up" with the system, as it takes time for the updated emissions in CHEEREIO to become consistent with observations. If this option is set to "true", then at the end of the burn-in period (given by ``BURN_IN_END``) CHEEREIO will scale the ensemble mean to match the observational mean, as in the ``SIMPLE_SCALE_FOR_FIRST_ASSIM_PERIOD`` option. This ensures that any biases introduced in the period where CHEEREIO emissions are "catching up" with observations are corrected. For more information on the burn-in period, see :ref:`Burn in period`.
* BURN_IN_END: If ``SIMPLE_SCALE_AT_END_OF_BURN_IN_PERIOD`` is set to ``true``, then this is the date (YYYYMMDD) when the burn-in period ends
* POSTPROCESS_START_DATE: The date when the postprocessing script should start (YYYYMMDD). This should always be at least one assimilation window away from ``START_DATE``. If you are using a burn-in period, you can set this for after the burn-in period ends to ensure that all your analysis discards this period.
* POSTPROCESS_END_DATE: The date when the postprocessing script should end (YYYYMMDD). Usually the same as ``END_DATE``, though the user can change the postprocess start and end dates to fit whatever application they are interested in.
* nEnsemble: Number of ensemble members. 32 or 48 are usually good numbers. This number of run directories will be created in the ``ensemble_runs`` folder and will be run simultaneously.
* SIMULATE_NATURE: *Deprecated: will be removed before official release (DP, 2022/03/09)*. End users should leave this set to "false", as this was used for testing in early CHEEREIO development. 
* verbose: Amount of information to print out as the ensemble runs. 1 is the default. 0 supresses most output, 2 is useful for basic debugging, and 3 for intense debugging. 

Cluster settings
~~~~~~~~~~~~~

The next section of the ``ens_config.json`` file controls settings that will be used when submitting jobs to the scheduler. These settings overwrite the template batch submission scripts included with CHEEREIO.

* NumCores: Number of cores used in each of the ensemble runs. CHEEREIO also will use these cores to parallelize assimilation computation columnwise.
* NumCtrlCores: Number of cores to use in the control run simulation, if using.
* Partition: Partition of your cluster you are submitting to. At Harvard, ``huce_intel`` is a good choice.
* Memory: Memory in megabytes used by each ensemble member. CHEEREIO can be quite memory intensive because it loads in restarts and history files for many ensemble members in addition to observations, and sometimes produces large matrices, so expect to use more than in standard GEOS-Chem runs.
* EnsCtrlSpinupMemory: Memory in megabytes for ensemble spinup, control, and regular spinup simulations (i.e. those simulations without LETKF assimilation). Set as you would a normal GEOS-Chem simulation.
* WallTime: Time allowed for the overall assimilation process (runs and assimilation) to occur in format D-HH\\\\:MM. Assimilation adds substantial overhead so expect it to be slow.
* EnsSpinupWallTime: Time allowed for the ensemble spinup process (no assimilation, just running all ensemble members from ``ENS_SPINUP_START`` through ``ENS_SPINUP_END`` with scaling factors applied) in format D-HH\\\\:MM. If not using, you can just leave as an empty string.
* ControlWallTime: Wall time for the control run simulation, if you're using one. Empty string otherwise.
* SpinupWallTime: Wall time for the spinup simulation, if you're using one. Empty string otherwise.
* CondaEnv: The name of the Conda environment with all of the CHEEREIO packages installed. It is strongly recommended that you install an environment using the YAML file that ships with CHEEREIO in the ``environments/`` folder.
* AnimationEnv: The name of the Conda environment that has the tools necessary to make animated postprocessing plots. A YAML file will be added to the ``environments/`` folder before release giving this Conda environment.
* MaxPar: Maximum number of cores to use while assimilating columns in parallel using CHEEREIO, maxing out at ``NumCores``. Setting this number smaller than NumCores saves on memory but adds to the assimilation time. 

Species in state/control vectors
~~~~~~~~~~~~~

One useful feature of CHEEREIO is its distinction between "control" and "state" vectors. The state vector should consist of all concentrations relevant to the problem at hand as well as the emissions of interest (e.g. large chemical families). The control vector is a subset of the state vector, and represents concentrations and the same emissions of interest that the user believes can reasonably be updated on the basis of observations. In many cases the control vector and state vector are identical. However, in some cases removing species from the control vector can help CHEEREIO handle ensemble behavior. Although the entire state vector is used to calculate the concentration and emissions update, **only the control vector is actually updated.** In practice, this distinction helps tamp down on noise and create well-behaved assimilations.

* STATE_VECTOR_CONC: Species from the restart files to be included in the state vector. It is generally recommended to include a range of species that might affect the species you are mainly interested in, but not so large a range that you end up analyzing noise. Given as an array. This is an example for NO\ :sub:`x` data assimilation: 
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

* CONTROL_VECTOR_CONC: A subset of the state vector concentration species that will be updated by assimilation. Although an update for all members of the state vector will be calculated, only the species listed in this array will have that update saved. This allows a wide range of species to be considered in the update calculation process but only a smaller, more tightly coupled subset of species to actually be changed and passed to GEOS-Chem. The goal is to tamp down on noise. Again, in many simulations the state vector and control vector entries will be identical.
* STATE_VECTOR_CONC_REPRESENTATION: How are concentrations represented within the state vector? There are several options. "3D" puts all 3D concentrations of the species in ``STATE_VECTOR_CONC`` from the restart file into the state vector. ``column_sum`` computes the partial columns of the species of interest and then sums to obtain units molec/cm2.``trop_sum`` is identical to ``column_sum`` except that it only represents the tropospheric column. Finally, ``column_sum_noheight`` is an option if the ``Met_BXHEIGHT`` is not saved in the GEOS-Chem restart. A simple sum of partial pressures is taken in this case, without weighting by box height. 
* CONTROL_VECTOR_EMIS: A dictionary linking a label for emissions scalings to the species emitted. For example, you could write ``"CH4_WET" : "CH4"`` to reference wetland methane emissions. CHEEREIO automatically will update ``HEMCO_Config.rc`` accordingly, but cannot distinguish between different emissions of the same species on its own; the user has to manually edit ``HEMCO_Config.rc`` to correct this if distinguishing between different sources of the same species. For a more thorough explanation, see the entry in :ref:`Template`. You can also use one label to link to emissions of multiple species, meaning that all these emissions will be controlled by one scaling factor file, such as ``"NOx":["NO","NO2"]`` to indicate that NO and NO\ :sub:`2`\ are controlled by one scaling factor file. Here are a couple of examples:
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
* SPINUP_HISTORY_freq: Frequency of history output files saved during ensemble spinup (i.e. when executing ``run_ensspin.sh``). Often set to be a longer period to save memory.
* SPINUP_HISTORY_dur: As in ``SPINUP_HISTORY_freq``, but for duration.
* SaveLevelEdgeDiags: Should the LevelEdgeDiags collection be turned on? "True" or "False". This is mandatory for assimilating most forms of satellite data.
* SaveStateMet: Should the StateMet collection be turned on? "True" or "False". This is mandatory for assimilating some forms of satellite data, like OMI NO\ :sub:`2`\ .
* SaveArea: Should grid cell areas be used in the assimilation process? "True" or "False".
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

* HistoryLevelEdgeDiagsToSave: A list of data to save in the LevelEdgeDiags collection. Just ``Met_PEDGE`` is sufficient for many forms of assimilation.
* HistoryStateMetToSave: A list of data to save in the StateMet collection. Below is an example of necessary fields for assimilating OMI NO\ :sub:`2`\ .
::

	"HistoryStateMetToSave" : [
		"Met_TropLev",
		"Met_BXHEIGHT",
		"Met_T"
	],

Observation settings
~~~~~~~~~~~~~

* OBSERVED_SPECIES: A dictionary linking a label for observations with the species observed. For example, you could write ``"NO2_SATELLITE" : "NO2"`` to reference satellite observations of NO2. 
* OBS_TYPE: A dictionary linking a label for observations with the observer type, so that CHEEREIO knows how to interpret observation files. One entry is required for each entry in ``OBSERVED_SPECIES``, with every key from ``OBSERVED_SPECIES`` represented here. Valid values include "OMI" and "TROPOMI", or any other observation operator type added to CHEEREIO by you or by other users. Instructions on how to add an observation operator to CHEEREIO such that it can be switched on from ``OBS_TYPE`` in the configuration file are given in the :ref:`New observation:` page.
* TROPOMI_dirs: Dictionary linking observed TROPOMI species to the directory containing the observations. If you aren't using TROPOMI, this can be left blank. Here is an example, along with the corresponding ``OBSERVED_SPECIES`` settings:
::

	"OBSERVED_SPECIES" : {
		"CH4_TROPOMI": "CH4"
	},
	"TROPOMI_dirs" : {
		"CH4" : "/n/holylfs05/LABS/jacob_lab/dpendergrass/tropomi/NO2/2019"
	},

* OMI_dirs: As in TROPOMI_dirs, but for OMI. Note that other observation operators can be added as separate key entries in this configuration file by following the instructions on the :ref:`New observation:` page. 
* SaveDOFS: Should CHEEREIO calculate and save the Degrees of Freedom for Signal (DOFS), or the trace of the observing system averaging kernel matrix? Note that since the prior error covariance matrix is not invertible because of our ensemble approach the pseudoinverse is used instead. See section 11.5.3 of Brasseur and Jacob for more information. The idea here is that if there is not enough information in a localized assimilation calculation we should set the posterior equal to the prior. *Note: this option is functional but DOFS values are not easily interpretable; hold off use for now while we think of alternative definitions in our rank-deficient space (DP, 2022/11/07).*
* DOFS_filter: What is the minimum DOFS for a localized region for the assimilation to be saved? If DOFS is below this threshold the posterior is set equal to the prior.

Scaling factor settings
~~~~~~~~~~~~~

* pPERT: Setting for initial emissions scaling factor creation, where the number provided :math:`p` is used to generate random scaling factors from a distribution specified in the ``pertType`` and ``correlatedInitialScalings`` entries. This value is specified by a dictionary where the key is one of the keys in ``CONTROL_VECTOR_EMIS`` and the value is a float. Every key in ``CONTROL_VECTOR_EMIS`` must be represented here.
* pertType: How ``pPERT`` is interpreted, specified as a dictionary where the key one of the keys in ``CONTROL_VECTOR_EMIS``. Valid entries are "exp", "percent", and "std". 

   * For "exp", the ``pPERT`` entry :math:`p` is used to generate random scaling factors from the distribution :math:`p^u,\\ u{\sim}U(-1,1)`, meaning that u is a uniform random variable ranging from -1 to 1. For example, if ``pPERT`` is "4" then scalings will range from 0.25 to 4, centered on 1. Note that for this option, there is no correlation between neighboring scale factors (i.e. they are i.i.d.). 
   * For "percent",  the ``pPERT`` entry :math:`p` is used to generate random scaling factors from the distribution :math:`u{\sim}U(1-p,1+p)`, meaning that u is a uniform random variable ranging from :math:`1-p` to :math:`1+p`. For example, if ``pPERT`` is "0.5" then scalings will range from 0.5 to 1, centered on 1.  Like "exp", for this option there is no correlation between neighboring scale factors (i.e. they are i.i.d.). 
   * For "std",  there are two possible interpretations depending on the corresponding value in ``correlatedInitialScalings``. If the corresponding value in ``correlatedInitialScalings`` is ``False``, then the  ``pPERT`` entry :math:`p` is used to generate random scaling factors from the distribution :math:`n{\sim}N(1,p)`, meaning that n is a normal random variable with mean 0 and standard deviation :math:`p`. For example, if ``pPERT`` is "0.25" then scalings will have mean 1 and standard devation 0.25.  Like "exp", for this option there is no correlation between neighboring scale factors (i.e. they are i.i.d.). 
   * The other "std" option, where ``correlatedInitialScalings`` is ``True``, means that we sample scaling factors from a multivariate normal distribution. If ``speedyCorrelationApprox`` is turned off, then this is a true sampling from a multivariate normal distribution. The mean is a vector of ones and the covariance matrix is generated with exponentially decaying correlation as a function of distance. More specifically, the covariance between points :math:`a` and :math:`b` with distance :math:`d` kilometers between them is given by :math:`\exp(-d^2/(2*c))` where :math:`c` is a correlation distance constant in kilometers given by the ``corrDistances`` entry (specified as a dictionary with keys from ``CONTROL_VECTOR_EMIS``). However, sampling a multivariate normal distribution like this can take an extraordinarily large amount of time and memory. CHEEREIO includes a very fast approximation to a multivariate normal distribution which is generated by applying a Gaussian blur to uncorrelated noise. To use this option, set ``speedyCorrelationApprox`` to ``True`` (strongly recommended for any resolution higher than 4x5).

* correlatedInitialScalings: See above entry for ``pertType``.
* corrDistances: See above entry for ``pertType``.
* speedyCorrelationApprox: See above entry for ``pertType``. A suitable entry for an approximated multivariate normal sample is given below:
::

	"pPERT" : {
		"CH4":"0.1"
	},
	"pertType" : {
		"CH4":"std"
	},
	"correlatedInitialScalings" : {
		"CH4":"True"
	},
	"corrDistances" : {
		"CH4":"500"
	},
	"speedyCorrelationApprox" : "True",

* MaskOceanScaleFactor: Should scaling factors be allowed to vary over the oceans? A dictionary with keys matching ``CONTROL_VECTOR_EMIS`` and with values of "True" or "False" for each entry. If "True", scaling factors for that species over the ocean are always set to 1 across all ensemble members (i.e. no assimilation calculated).
* MaskCoastsGT25pctOcean: Should we use a looser definition of ocean, including grid cells with at least 25% ocean within the definition? "True" or "False". If oceans are masked, setting this to "True" eliminates many coastal cells which can have problematic satellite retrievals for some products.
* Mask60NScaleFactor: Should scaling factors above 60 N always be set to 1? A dictionary with keys matching ``CONTROL_VECTOR_EMIS`` and with values of "True" or "False" for each entry.
* Mask60SScaleFactor: Should scaling factors below 60 S always be set to 1? A dictionary with keys matching ``CONTROL_VECTOR_EMIS`` and with values of "True" or "False" for each entry.
* MinimumScalingFactorAllowed: What is the minimum scaling factor allowed? A dictionary with keys matching ``CONTROL_VECTOR_EMIS`` and with float values for each entry. Set to "nan" if no minimum scaling factor is enforced.
* MaximumScalingFactorAllowed: As above, but for the maximum scaling factors allowed.
* InflateScalingsToXOfPreviousStandardDeviation: CHEEREIO includes support for inflating posterior scaling factor standard deviations to a certain percentage of the initial standard deviation. A dictionary with keys matching ``CONTROL_VECTOR_EMIS`` and with float values for each entry, where "0.3" corresponds with inflating to 30% of the initial standard deviation (the recommended value). Set to "nan" to ignore. 
* MaximumScaleFactorRelativeChangePerAssimilationPeriod: The maximum relative change per assimilation period allowed for scaling factors. For example, a value "0.5" means that no more than a 50% change is allowed for a given scaling factor in a given assimilation period. A dictionary with keys matching ``CONTROL_VECTOR_EMIS`` and with float values for each entry, where "nan" ignores this setting.

LETKF settings
~~~~~~~~~~~~~

* REGULARIZING_FACTOR_GAMMA: A dictionary of regularization factors, with a key corresponding with each key in ``OBSERVED_SPECIES``, which inflates observed error covariance by a factor of :math:`1/\gamma`.
* OBS_ERROR: An dictionary of error information, with a key corresponding with each key in ``OBSERVED_SPECIES`` and a float value. The value is interpreted as one of three categories: "relative", "absolute", or "product". This information represents uncertainty in observations. If error is relative, it is given as a decimal (0.1 means 10% relative error). If error is absolute, it is given as the same units as the observations are in (CHEEREIO will square these values for the covariance matrix). If error is "product," then CHEEREIO uses the error from the observation product. In the product case, the number recorded under OBS_ERROR will not be used. For clarity, only diagonal observational covariance matrices are supported at this time.
* OBS_ERROR_TYPE: A dictionary of error types, with values given as strings reading "relative",  "absolute", or "product", and with keys corresponding to each key in ``OBSERVED_SPECIES``. This tells CHEEREIO how to interpret the error data types, as described above.
* OBS_ERROR_SELF_CORRELATION: A dictionary of correlations between errors in data samples, with a key corresponding with each key in ``OBSERVED_SPECIES`` and a float value. This value is used to reduce error if the user would like to aggregate multiple observations together onto the GEOS-Chem grid ("super-observations"). More on this below in the ``AV_TO_GC_GRID`` entry. 
* MIN_OBS_ERROR: A dictionary of minimum possible errors, with a key corresponding with each key in ``OBSERVED_SPECIES`` and a float value. If the user would like to aggregate multiple observations together onto the GEOS-Chem grid ("super-observations"), this value gives the minimum possible error allowable upon error reduction. More on this below in the ``AV_TO_GC_GRID`` entry. 
* OTHER_OBS_ERROR_PARAMETERS: A dictionary of dictionaries, with a key corresponding with each key in ``OBSERVED_SPECIES`` and a value that itself is a dictionary with additional settings and their values. At this time, the only setting that is applied using this entry is called ``transport_error``, which is used to account for perfectly correlated model transport errors when the user aggregates multiple observations together onto the GEOS-Chem grid ("super-observations"). More information on this in the the ``AV_TO_GC_GRID`` entry. Below is valid syntax for this setting:
::

	"OTHER_OBS_ERROR_PARAMETERS":{
		"CH4_TROPOMI":{
			"transport_error":"6.1"
		}
	},

* AV_TO_GC_GRID: "True" or "False", should observations be averaged to the GEOS-Chem grid? If "false", the above three entries and the below entry are all ignored. The use of "super observations" is a useful technique to balance prior and observational errors while also reducing the computational complexity of the optimization (by reducing the size of the observational vectors and matrices in the LETKF calculation). The main subtlety that needs to be handled for this super observation aggregation is the adjustment of observational error. Users can specify one of several error reduction functions listed below, specified in the ``SUPER_OBSERVATION_FUNCTION`` entry.

   * "sqrt": A modified version of the familiar square root law, where if we aggregate :math:`n` observations (indexed by :math:`i`) with errors :math:`\sigma_i` together, the new error is :math:`\bar{\sigma}/\sqrt{n}` where :math:`\bar{\sigma}` is the mean of the :math:`\sigma_i`. The modification accounts for correlations :math:`c` between errors (e.g. due to correlated retrieval errors from shared surface type or similar albedo), and for a user-specified minimum error :math:`\sigma_{\min}`. Thus the equation that is actually applied is given by :math:`\max\left[\left(\bar{\sigma}\cdot\sqrt{\frac{1-c}{n}+c}\right),\sigma_{\min}\right]`. The correlation :math:`c` is taken from ``OBS_ERROR_SELF_CORRELATION`` with default value 0, and the minimum error :math:`\sigma_{\min}` is taken from ``MIN_OBS_ERROR`` with default value 0 (i.e. the normal square root law).
   * "default": As with "sqrt", but with an additional term accounting for the fact that GEOS-Chem transport errors are perfectly correlated. Because perfectly correlated errors are irriducible no matter how many realizations are averaged, the resulting equation is given by :math:`\max\left[\sqrt{\bar{\sigma}^2\cdot\left(\frac{1-c}{n}+c\right)+\sigma_t^2},\sigma_{\min}\right]` where :math:`\sigma_t` is transport error supplied by the "transport_error" entry from ``OTHER_OBS_ERROR_PARAMETERS`` etnry. 
   * "constant": No error reduction applied. In other words, no matter how many observations are averaged, this function just returns :math:`\bar{\sigma}`. 

* SUPER_OBSERVATION_FUNCTION: A dictionary with a key corresponding with each key in ``OBSERVED_SPECIES``, and a value corresponding to one of the super observation error reduction functions listed in the ``AV_TO_GC_GRID`` entry. Users can add new superobservation functions within the ``produceSuperObservationFunction`` closure in the ``observation_operators.py`` file and activate them from this entry; see the :ref:`New superobservation` entry for more information. 
* INFLATION_FACTOR: :math:`\rho-1` from Hunt et. al. (2007). A small number (start with something between 0 and 0.1 and slowly increase according to testing) that inflates the ensemble range. In ensemble Kalman filters, uncertainty usually decreases too quickly and must manually be reinflated.
* ASSIM_TIME: Length in hours of assimilation window. The assimilation window refers to the period in which GEOS-Chem is run and observations are accumulated; the data assimilation update is calculated in one go within this window. The data assimilation literature contains extensive discussion of this concept.
* MAXNUMOBS: Maximum number of observations used in a column assimilation calculation. If the number of observations available is greater than this value, then CHEEREIO will randomly throw out observations until only ``MAXNUMOBS`` remain.
* MINNUMOBS: Minimum number of observations for a column assimilation calculation to be performed. If the number of observations is below this number, no assimilation is calculated and the posterior is set to the prior.
* LOCALIZATION_RADIUS_km: When updating a column, CHEEREIO only considers data and observations within this radius (in kilometers).
* AveragePriorAndPosterior: "True" or "False", should the posterior be set to a weighted average of the prior and the posterior calculated in the LETKF algorithm? If set to true, the prior weight in the average is given by ``PriorWeightinPriorPosteriorAverage`` in the next setting.
* PriorWeightinPriorPosteriorAverage: The prior weight if averaging with the posterior from the LETKF. A value between 0 and 1.

Postprocessing settings
~~~~~~~~~~~~~

* animation_fps_scalingfactor: Frames per second for movies of scaling factors and emissions made by the postprocessing workflow.
* animation_fps_concentrations: Frames per second for movies of concentrations made by the postprocessing workflow.
* hemco_diags_to_process: An array of entries from HEMCO Diagnostics that you would like processed into movies. This is usually emissions totals. See below for an example entry:
:: 

	"hemco_diags_to_process" : [
		"EmisCH4_Total"
	],

* OBSERVATION_UNITS: a dictionary with keys from ``OBSERVED_SPECIES`` and values representing the units will be plotted. This is governed by how the observation operator is defined. 
* scalefactor_plot_freq: The CHEEREIO postprocessing routine will save out maps of scale at this temporal resolution: either "all" for save out an image for every assimilation window, or "monthly" to save out one per month.

Extensions
~~~~~~~~~~~~~

Additional CHEEREIO settings, usually for specific observation types, can be loaded in through extensions. Extensions in CHEEREIO are extra JSON files that store additional settings in order to prevent clutter in the ``ens_config.json`` file. Extensions can easily be added by saving a file with name ``NAME_extension.json`` within the ``extensions/`` folder. To load in the settings within the ``NAME_extension.json`` file, add the key NAME to the "Extensions" dictionary in ``ens_config.json`` with value "True". Below is an example where we load in the settings in the ``TROPOMI_ALL_extension.json``, ``TROPOMI_CH4_extension.json``, and ``CH4_extension.json`` files. 
::

	"Extensions": {
		"TROPOMI_ALL":"True",
		"TROPOMI_CH4":"True",
		"CH4":"True"
	}

Below we list the settings that you can set with extensions.

* TROPOMI_CH4 extension.

   * TROPOMI_CH4_FILTERS: Apply specialized filters for TROPOMI methane? Set to "True" if doing a TROPOMI methane inversion, otherwise set to "False".
   * TROPOMI_CH4_filter_blended_albedo: Filter out TROPOMI methane observations with a blended albedo above this value. Set to "nan" to ignore.
   * TROPOMI_CH4_filter_swir_albedo_low: Filter out TROPOMI methane observations with a SWIR albedo below this value. Set to "nan" to ignore.
   * TROPOMI_CH4_filter_swir_albedo_high: Filter out TROPOMI methane observations with a SWIR albedo above this value. Set to "nan" to ignore.
   * TROPOMI_CH4_filter_winter_lat: Filter out TROPOMI methane observations beyond this latitude in the winter hemisphere. Set to "nan" to ignore.
   * TROPOMI_CH4_filter_roughness: Filter out TROPOMI methane observations with a surface roughness above this value. Set to "nan" to ignore.
   * TROPOMI_CH4_filter_swir_aot: Filter out TROPOMI methane observations with a SWIR AOT above this value. Set to "nan" to ignore.

* TROPOMI_ALL extension.

   * postprocess_save_albedo: Should the postprocessing workflow save out albedo? "True" or "False".

* CH4 extension.

   * USE_CUSTOM_CH4_OH_ENTRY: Should we overwrite the OH field setting in ``HEMCO_Config.rc`` for the specialty CH4 simulation with a user-specified entry (given below)? "True" or "False".
   * CUSTOM_CH4_OH_ENTRY: If USE_CUSTOM_CH4_OH_ENTRY is True, then we overwrite the OH field line in ``HEMCO_Config.rc`` with this entry. Note that backslashes need to be escaped with a front slash. Here is an example entry: ``* GLOBAL_OH  $ROOT\/OH\/v2014-09\/v5-07-08\/OH_3Dglobal.geos5.47L.4x5.nc OH           1985\/1-12\/1\/0 C xyz kg\/m3 * - 1 1`` 

* OMI_NO2 extension.

   * OMI_NO2_FILTERS: Apply specialized filters for OMI NO2? Set to "True" if doing an OMI NO2 inversion, otherwise set to "False".
   * OMI_NO2_filter_sza: Filter out OMI NO2 observations with a solar zenith angle above this value. Set to "nan" to ignore.
   * OMI_NO2_filter_cloud_radiance_frac: Filter out OMI NO2 observations with a cloud radiance fraction above this value. Set to "nan" to ignore.
   * OMI_NO2_filter_surface_albedo: Filter out OMI NO2 observations with a surface albedo above this value. Set to "nan" to ignore.

