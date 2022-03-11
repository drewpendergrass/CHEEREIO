Field guide to all core CHEEREIO scripts
==========

The ``core/`` folder of the main CHEEREIO code directory is where the, well, core of CHEEREIO resides. All the main code utilities, assimilation tools, and observation integration scripts live here. This page gives an overview of all the files in the ``core`` folder; the most important files and files which users might interact with are explored on other pages in the documentation, and are linked to when appropriate.

Installation and spinup scripts
-------------

change_histcollections_durfreq.sh
~~~~~~~~~~~~~

A shell script, called in the course of running ``run_ensemble_spinup_simulations.sh`` in the ensemble spinup process, that updates the duration and frequency by which GEOS-Chem output is saved from spinup mode (monthly mean output) to assimilation mode (much more frequent output).

change_histrst_durfreq.sh
~~~~~~~~~~~~~

A shell script, called in the course of running ``run_ensemble_spinup_simulations.sh`` in the ensemble spinup process, that updates the duration and frequency by which GEOS-Chem restart files are saved from spinup mode to assimilation mode.

check_for_all_ensemble_spinup_restarts.sh
~~~~~~~~~~~~~

A shell script, called in the course of running ``run_ensemble_spinup_simulations.sh`` in the ensemble spinup process, that checks if all restart files expected to be present at the end of ensemble spinup have in fact been created. 

hemco_config_updater.py
~~~~~~~~~~~~~

A very short Python script that invokes ``hemco_utils.py`` and updates ``HEMCO_Config.rc`` to match user settings specified in ``ens_config.json``. This is called in the template run directory creation stage of ``setup_ensemble.sh``.

hemco_utils.py
~~~~~~~~~~~~~

A set of Python classes and functions designed to help parse, modify, and save ``HEMCO_Config.rc`` so that it matches user settings; in particular, these utilities are designed for linking CHEEREIO-generated scaling factors to key emissions.

initialize_scaling_factors.py
~~~~~~~~~~~~~

A Python script that creates a randomized initial set of scaling factors for each emissions grouping the user would like to assimilate, incorporating relevant user settings from ``ens_config.json``. This reflects the prior emissions scaling distribution and is called in the ensemble run directory creation stage of ``setup_ensemble.sh``.

prep_par.py
~~~~~~~~~~~~~

prepare_template_hemco_config.sh
~~~~~~~~~~~~~

update_history.py
~~~~~~~~~~~~~

Run management scripts
-------------

advance_timestep.py
~~~~~~~~~~~~~

check_and_complete_assimilation.sh
~~~~~~~~~~~~~

check_for_all_columns.py
~~~~~~~~~~~~~

check_for_all_restarts.sh
~~~~~~~~~~~~~

cleanup.sh
~~~~~~~~~~~~~

update_current_time.sh
~~~~~~~~~~~~~

update_input_geos.sh
~~~~~~~~~~~~~

Assimilation scripts
-------------

combine_columns_and_update.py
~~~~~~~~~~~~~

letkf_utils.py
~~~~~~~~~~~~~

par_assim.sh
~~~~~~~~~~~~~

par_letkf.py
~~~~~~~~~~~~~

toolbox.py
~~~~~~~~~~~~~

Observation scripts
-------------

omi_tools.py
~~~~~~~~~~~~~

tropomi_tools.py
~~~~~~~~~~~~~

Deprecated scripts
-------------

The following scripts have been deprecated and will be removed before the official release of CHEEREIO:

* diff_col.py
* observation_operators.py
* randomize_restarts.py
* regrid_landmask_fraction.py
* tropomi_loader.py
