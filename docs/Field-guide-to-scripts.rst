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

hemco_config_updater.py
~~~~~~~~~~~~~

hemco_utils.py
~~~~~~~~~~~~~

initialize_scaling_factors.py
~~~~~~~~~~~~~

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
