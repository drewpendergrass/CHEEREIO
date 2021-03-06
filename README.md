[![Documentation Status](https://readthedocs.org/projects/cheereio/badge/?version=latest)](https://cheereio.readthedocs.io/en/latest/?badge=latest)

# CHEEREIO

The CHEmistry and Emissions REanalysis Interface with Observations (CHEEREIO) is a set of Python and shell scripts that support data assimilation and emissions inversions for arbitrary runs of the GEOS-Chem chemical transport model via an ensemble approach (i.e. without the model adjoint). CHEEREIO follows five design principles:

1. **Easy to customize**: Assimilate anything, in any GEOS-Chem configuration or simulation.
2. **Easy to maintain**: Science automatically aligned with latest model version.
3. **Easy to deploy**: One configuration file controls installation and settings
4. **Easy to link observations**: Object-oriented observation operator implementation allows the user to rapidly add new kinds of data with minimal programming required.
5. **Quick to run**: Wall runtime should be no more than 2x longer than vanilla GEOS-Chem (4D-Var limit) assuming resources are available.

CHEEREIO is currently fully tested and operational for 3D LETKF data assimilation. Development of 4D assimilation and additional tools for observation operator creation are currently being tested for TROPOMI CH4 and in development for TROPOMI/OMI NO2. I am adding new satellite data as rapidly as possible. After some validation runs to confirm the workflow's scientific value, I will attempt to reduce memory and I/O time complexity by switching from GNU parallel parallelization to MPI via the ray package (allowing the sharing of large matrices in memory and removing the need to write columns to file).

## Documentation
Detailed documentation is available on [ReadTheDocs](https://cheereio.readthedocs.io), including installation instructions. If you encounter any problems not covered by the documentation, please open an [issue](https://github.com/drewpendergrass/CHEEREIO/issues) on GitHub.

## About CHEEREIO
The CHEmistry and Emissions REanalysis Interface with Observations (CHEEREIO) is a package that wraps the [GEOS-Chem](https://github.com/geoschem) chemical transport model source code. After a simple modification of a single configuration file (`ens_config.json`), CHEEREIO automatically produces and compiles a template GEOS-Chem run directory, which it then copies into an ensemble. Each ensemble member comes with a randomized set of gridded emissions scaling factors for species specified by the user. As the ensemble of runs progresses, CHEEREIO will periodically pause the ensemble, compare with a set of observations (i.e. satellite, surface, and/or aircraft), and update relevant emissions scaling factors and chemical concentrations to best match reality given the uncertainties of measurements and model. CHEEREIO calculates this update via the 4D Asynchronous Localized Ensemble Transform Kalman Filter (4D-LETKF) as described in [Hunt. et. al., \[2007\]](https://doi.org/10.1016/j.physd.2006.11.008). Because this approach is model agnostic (specifically, it does not rely on the adjoint), CHEEREIO supports emissions updates and chemical concentration corrections for arbitrary configurations of the GEOS-Chem model. However, the current CHEEREIO codebase assumes that GEOS-Chem code is version 13.0.0 or later.

