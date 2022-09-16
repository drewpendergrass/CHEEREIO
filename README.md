[![Documentation Status](https://readthedocs.org/projects/cheereio/badge/?version=latest)](https://cheereio.readthedocs.io/en/latest/?badge=latest)

# CHEEREIO

The CHEmistry and Emissions REanalysis Interface with Observations (CHEEREIO) is a set of Python and shell scripts that support data assimilation and emissions inversions for arbitrary runs of the GEOS-Chem chemical transport model via an ensemble approach (i.e. without the model adjoint). CHEEREIO follows four design principles:

1. **Easy to customize**: Assimilate anything, in any GEOS-Chem configuration or simulation.
2. **Easy to maintain**: Science automatically aligned with latest model version.
3. **Easy to deploy**: One configuration file controls installation and settings
4. **Easy to link observations**: Object-oriented observation operator implementation allows the user to rapidly add new kinds of data with minimal programming required.

The beta release for CHEEREIO and an accompanying model description is currently in prep. Until then, this codebase is under constant revision and re-engineering. At any one time the code here may be broken by an update only to be fixed later that day. Don't treat it as stable yet!

## Documentation
Detailed documentation is available on [ReadTheDocs](https://cheereio.readthedocs.io), including installation instructions. If you encounter any problems not covered by the documentation, please open an [issue](https://github.com/drewpendergrass/CHEEREIO/issues) on GitHub.

## About CHEEREIO
The CHEmistry and Emissions REanalysis Interface with Observations (CHEEREIO) is a package that wraps the [GEOS-Chem](https://github.com/geoschem) chemical transport model source code. After a simple modification of a single configuration file (`ens_config.json`), CHEEREIO automatically produces and compiles a template GEOS-Chem run directory, which it then copies into an ensemble. Each ensemble member comes with a randomized set of gridded emissions scaling factors for species specified by the user. As the ensemble of runs progresses, CHEEREIO will periodically pause the ensemble, compare with a set of observations (i.e. satellite, surface, and/or aircraft), and update relevant emissions scaling factors and chemical concentrations to best match reality given the uncertainties of measurements and model. CHEEREIO calculates this update via the 4D Asynchronous Localized Ensemble Transform Kalman Filter (4D-LETKF) as described in [Hunt. et. al., \[2007\]](https://doi.org/10.1016/j.physd.2006.11.008). Because this approach is model agnostic (specifically, it does not rely on the adjoint), CHEEREIO supports emissions updates and chemical concentration corrections for arbitrary configurations of the GEOS-Chem model. However, the current CHEEREIO codebase assumes that GEOS-Chem code is version 13.0.0 or later.

