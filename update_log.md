# Update log

## Version 1.1

**Action items for users updating to version 1.1:**

* Update ens_config.json by comparing with templates, as new fields are now expected!
* Update CHEEREIO conda environment using the cheereio.yml file in the environments folder, as new packages are required!

**Updates under version 1.1:**

* Add validation script for ens_config.json, to be filled out over time. Checks at install time if ensemble settings are valid, including some syntax!
* Add rerun support from *Varon et. al.*, 2023, which is a lower cost approximation of run-in-place type methods.
* Initial YAML support for eventual transition from JSON.
* Support for GC 14.1 and newer (SpeciesConc labels are changed, causing errors in CHEEREIO 1.0).
* Made it possible to use observations for validating the prior/posterior results on the fly without assimilation (e.g. for validation).
* Improvements to the postprocessing workflow. 
* Add ObsPack support for GC >=14.1.0
  * New ObsPack observation operator released and ready for either use in assimilation or validation.

## Version 1.0

* Initial release, as described in the *Pendergrass et. al.*, 2023 model description paper.
