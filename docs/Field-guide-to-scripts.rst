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

A Python script that prepares LETKF parallelization in advance of any assimilation. This is done by dividing up the columns that will be assimilated by each core in each ensemble run job (LETKF is an "embarassingly parallel" algorithm and requires no coordination between columns at assimilation time). This division of columns is stored in the ``scratch/`` directory and is consulted by each core at run time to ensure each column is processed exactly once. The script is called in the ensemble run directory creation stage of ``setup_ensemble.sh``.

prepare_template_hemco_config.sh
~~~~~~~~~~~~~

A simple wrapper shell script, called by ``setup_ensemble.sh`` in the template run directory creation stage, that in turn calls ``hemco_config_updater.py`` within an appropriate conda environment.

update_history.py
~~~~~~~~~~~~~

A Python toolkit and set of scripts designed to align the ``HISTORY.rc`` output settings with CHEEREIO's needs at various stages of the installation, spinup, and assimilation processes. This script is called in several places by ``setup_ensemble.sh``, ``change_histcollections_durfreq.sh``, and ``change_histrst_durfreq.sh`` to update CHEEREIO output settings at different stages of execution.

Run management scripts
-------------

advance_timestep.py
~~~~~~~~~~~~~

This short Python script called by ``update_input_geos.sh`` at the end of assimilation, which advances the ensemble timestep stored in the ``scratch/`` directory. It also checks if the simulation is complete, and if so produces the file ``ENSEMBLE_COMPLETE`` stored in ``scratch/``, which terminates assimilation.

check_and_complete_assimilation.sh
~~~~~~~~~~~~~

A shell script that calls the Python script ``check_for_all_columns.py`` to see if all expected assimilated columns (with extension ``.npy`` are present in the ``scratch/`` folder. If they are, execute the Python script ``combine_columns_and_update.py`` to update NetCDF files.

check_for_all_columns.py
~~~~~~~~~~~~~

A brief Python script which counts the number of ``.npy`` files present in the ``scratch/`` folder, and checks if it matches the total number of columns that need to be assimilated. If all expected files are present, it writes a file called ``ALL_COLUMNS_FOUND`` into the ``scratch/`` folder, signalling to all runs that it is time to complete assimilation.

check_for_all_restarts.sh
~~~~~~~~~~~~~

A shell script which checks if all expected restarts are present with a timestamp corresponding to the end of the current GEOS-Chem run period. If all expected restarts are present, the script writes a file called ``ALL_RUNS_COMPLETE`` into the ``scratch/`` folder. This file's presence means that all ensemble members have finished running their respective GEOS-Chem simulations and the assimilation step can begin. 

cleanup.sh
~~~~~~~~~~~~~

A shell script called after assimilation has fully completed (i.e., all restart files and scaling factors are updated with the posterior results). This script (1) removes all assimilated columns and signal files from ``scratch/``, and (2) calls ``update_current_time.sh`` and ``update_input_geos.sh`` which prepare the GEOS-Chem input files for the next run. The removal of signal files like ``ALL_RUNS_COMPLETE`` indicate to the ensemble run script that it is safe to start GEOS-Chem again.

update_current_time.sh
~~~~~~~~~~~~~

A very brief shell script called at the very end of assimilation by ``cleanup.sh``, which updates the file ``CURRENT_DATE_TIME`` in ``scratch/`` so that it contains the start date for the upcoming GEOS-Chem run. 

update_input_geos.sh
~~~~~~~~~~~~~

A shell script which (1) calls ``advance_timestep.py`` to update the internal time stored in the ``scratch/`` directory, and (2) uses that updated internal time to update the ``input.geos`` file across the ensemble.

Assimilation scripts
-------------

combine_columns_and_update.py
~~~~~~~~~~~~~

If the script ``check_and_complete_assimilation.sh`` finds that all expected ``.npy`` files containing assimilated columns are present in ``scratch/``, then this Python script is called. This script gathers the assimilated columns and loads in all the ensemble restarts and scaling factors, uses the contents of the columns to update restarts and scaling factors, and then writes the updated data to disk.

letkf_utils.py
~~~~~~~~~~~~~

This long Python file is the core of CHEEREIO, and is described in detail in the REMOVED, UPDATE page. It contains complex classes and associated methods that do the IO and associated calculations required for the LETKF algorithm.

par_assim.sh
~~~~~~~~~~~~~

A wrapper shell script that calls ``par_letkf.py`` within the appropriate conda environment, passes information to the Python script ensuring that the appropriate set of columns are assimilated, and logs errors that occur in the assimilation process.

par_letkf.py
~~~~~~~~~~~~~

A short Python script, many instantiations of which are run in parallel, that creates relevant objects and calls methods from ``letkf_utils.py`` to assimilate the set of columns assigned to a particular core or set of cores.

toolbox.py
~~~~~~~~~~~~~

Basic utilities including distance calculations, JSON file I/O, and indexing support that are used across CHEEREIO Python scripts. 

Observation scripts
-------------

omi_tools.py
~~~~~~~~~~~~~

This long Python file includes tools and classes necessary for interfacing with OMI satellite products, and is described in detail in the :ref:`OMI tools` entry.

tropomi_tools.py
~~~~~~~~~~~~~

This long Python file includes tools and classes necessary for interfacing with TROPOMI satellite products, and is described in detail in the :ref:`TROPOMI tools` entry.

Deprecated scripts
-------------

The following scripts have been deprecated and will be removed before the official release of CHEEREIO:

* diff_col.py
* observation_operators.py
* randomize_restarts.py
* regrid_landmask_fraction.py
* tropomi_loader.py
