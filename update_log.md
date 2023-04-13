# Update log

## Version 1.1 goals (under development)

* Add validation script for ens_config.json, to be filled out over time. Checks at install time if ensemble settings are valid, including some syntax!
* Add rerun support from *Varon et. al.*, 2023, which is a lower cost approximation of run-in-place type methods.
* Initial YAML support for eventual transition.
* Made it possible to compare some observations to prior/posterior on the fly without assimilation (e.g. for validation).
* Add ObsPack support for GC >=14.0.0
  * New ObsPack observation operator released and ready for either use in assimilation or validation.
  * New ObsPack postprocess workflow for validation.

## Version 1.0

* Initial release, as described in the *Pendergrass et. al.*, 2023 model description paper.