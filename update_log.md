# Update log

## Version 1.4.0 (rc)

**Action items for users updating to version 1.4.0:**

* Update ens_config.json by comparing with templates, as new fields are now expected!

**Updates under version 1.4.0:**

* Added TCCON and TROPOMI CO operators (thanks to Sina Voshtani)
* Added IASI NH3 operator.
* Bug fixes for SatDiagn assimilation, and for scale factor initialization (thanks to Yunxiao Tang).
* Significant performance and memory improvements (halving my methane assimilation time)
* Add option, in super observation case, to average data before applying observation operator for performance improvements.
* Streamlined running nested grid simulations.
* Removed DOFS pseudoinverse approach (long deprecated, never used) and replaced with an information quantification metric from Voshtani et al., 2025 (in review).
* Improved Python scripting so that code works across a wider variety of shells. 
* Additional bug fixes and improvements.

## Version 1.3.1

**Action items for users updating to version 1.2.1:**

* None. This is a small update to fix bugs.

**Updates under version 1.3.1:**

* Bug leading some assimilations to fail is fixed.

## Version 1.3.0

**Action items for users updating to version 1.3.0:**

* Update ens_config.json by comparing with templates, as new fields are now expected!

**Updates under version 1.3.0:**

* Bug fixes for regional inversions (thanks to Yunxiao Tang)
* Added experimental support for inversions utilizing the SatDiagn collection in GEOS-Chem
* Added relaxation to prior spread (RTPS) inflation technique
* Added ability to amplify spreads and execute inflation methods like RTPS for species not in the state vector.
* Added convenience restore_backup.batch script (stored in scratch) to overwrite a failed run with the spun-up backup.
* Brought tests folder up to date and added more comprehensive testing of assimilation workflow.
* Improved CHEEREIO error catching and run control.
* Fixed CHEEREIO conda environment to replace removed packages (thanks to Lee Murray)
* CHEEREIO now installs a copy of LETKF model code in the run folder, so multiple runs can be tested at once and modified independently from one CHEEREIO code folder (analagous to GC run directory system).
* Added convenience copy_backup_into_new_ensemble.batch script (stored in scratch) to duplicate a spun-up backup into a new ensemble, for easy sensitivity simulations

## Version 1.2.1

**Action items for users updating to version 1.2.1:**

* None. This is a small update to fix bugs.

**Updates under version 1.2.1:**

* Bug leading some assimilations to fail is fixed.


## Version 1.2.0

**Action items for users updating to version 1.2:**

* Update ens_config.json by comparing with templates, as new fields are now expected!

**Updates under version 1.2.0:**

* TROPOMI CO operator added for the tagCO simulation (thanks to Sina Voshtani)
* Revamped the postprocess workflow
  * Add option to plot sites (e.g. obspack) rather than gridding everything.
  * Remove separate TROPOMI albedo postprocessing routines and replace with generalized method for adding additional data fields to be plotted. 
  * Postprocess workflow debugged for experiments with multiple observers.
* Extended rerun option to support a reprocessing an arbitrary number of assimilation windows.
* Additional options for burn-in periods
* Option to add additional uncorrelated initial perturbations to scaling factor fields according to prior emissions.
* Option to smooth included observations by distance weighting according to the gaspari cohn function.
* Simplified CHEEREIO use of HEMCO_Config.rc (only one version generated in template_run, rather than two as in previous versions).
* Improved and expanded documentation.
  * Compatibility updates with ReadTheDocs (thanks to Bob Yantosca)
* Bug fixes, including to the ObsPack operator (thanks to Lee Murray)


## Version 1.1.0

**Action items for users updating to version 1.1.0:**

* Update ens_config.json by comparing with templates, as new fields are now expected!
* Update CHEEREIO conda environment using the cheereio.yml file in the environments folder, as new packages are required!

**Updates under version 1.1.0:**

* Add validation script for ens_config.json, to be filled out over time. Checks at install time if ensemble settings are valid, including some syntax!
* Add rerun support from *Varon et. al.*, 2023, which is a lower cost approximation of run-in-place type methods.
* Initial YAML support for eventual transition from JSON.
* Support for GC 14.1 and newer (SpeciesConc labels are changed, causing errors in CHEEREIO 1.0).
* Made it possible to use observations for validating the prior/posterior results on the fly without assimilation (e.g. for validation).
* Improvements to the postprocessing workflow. 
* Add ObsPack support for GC >=14.1.0
  * New ObsPack observation operator released and ready for either use in assimilation or validation.

## Version 1.0.0

* Initial release, as described in the *Pendergrass et. al.*, 2023 model description paper.
