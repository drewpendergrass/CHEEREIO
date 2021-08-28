About CHEEREIO
==========

What CHEEREIO is
-------------

The CHEmistry and Emissions REanalysis Interface with Observations (CHEEREIO) is a package that wraps the `GEOS-Chem <https://github.com/geoschem>`__ chemical transport model source code. After a simple modification of a single configuration file (``ens_config.json``), CHEEREIO automatically produces and compiles a template GEOS-Chem run directory, which it then copies into an ensemble. Each ensemble member comes with a randomized set of gridded emissions scaling factors for species specified by the user. As the ensemble of runs progresses, CHEEREIO will periodically pause the ensemble, compare with a set of observations (i.e. satellite, surface, and/or aircraft), and update relevant emissions scaling factors and chemical concentrations to best match reality given the uncertainties of measurements and model. CHEEREIO calculates this update via the 4D Asynchronous Localized Ensemble Transform Kalman Filter (4D-LETKF) as described in `Hunt. et. al., [2007] <https://doi.org/10.1016/j.physd.2006.11.008>`__. Because this approach is model agnostic (specifically, it does not rely on the adjoint), CHEEREIO supports emissions updates and chemical concentration corrections for arbitrary configurations of the GEOS-Chem model. However, the current CHEEREIO codebase assumes that GEOS-Chem code is version 13.0.0 or later.

How CHEEREIO works
-----------------------------

TKTKTKTKT

