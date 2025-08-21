Guide to CHEEREIO core files
==========

The ``core/`` folder of the main CHEEREIO code directory is where the, well, core of CHEEREIO resides. All the main code utilities, assimilation tools, and observation integration scripts live here. This page gives an overview of all the files in the ``core`` folder; the most important files and files which users might interact with are explored on other pages in the documentation, and are linked to when appropriate.

Installation and spinup scripts
-------------

change_hemcodiag_freq.sh
~~~~~~~~~~~~~

A shell script, called in the course of running ``run_ensemble_spinup_simulations.sh`` in the ensemble spinup process or ``setup_ensemble.sh`` in the installation process depending on user settings, that updates the duration and frequency by which HEMCO Diagnostic output (e.g. emissions) is saved from spinup mode (less output) to assimilation mode (much more frequent output).

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

setup_obs_dates.py
~~~~~~~~~~~~~

To save time during assimilation, CHEEREIO produces a Python dictionary linking each observational file to the time period it covers. This script produces that dictionary and saves it to the ``scratch/`` folder. 

update_history.py
~~~~~~~~~~~~~

A Python toolkit and set of scripts designed to align the ``HISTORY.rc`` output settings with CHEEREIO's needs at various stages of the installation, spinup, and assimilation processes. This script is called in several places by ``setup_ensemble.sh``, ``change_histcollections_durfreq.sh``, and ``change_histrst_durfreq.sh`` to update CHEEREIO output settings at different stages of execution.


validate_ensconfig.py
~~~~~~~~~~~~~

A Python script that checks for common errors in the ens_config.json file before installation.

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

check_for_all_runs_complete.py
~~~~~~~~~~~~~

A Python script which checks if all GEOS-Chem runs are complete. The function is very similar to ``check_for_all_restarts.sh`` above, but it is robust to more exotic running-in-place settings; CHEEREIO automatically uses this script if necessary. 

cleanup.sh
~~~~~~~~~~~~~

A shell script called after assimilation has fully completed (i.e., all restart files and scaling factors are updated with the posterior results). This script (1) removes all assimilated columns and signal files from ``scratch/``, and (2) calls ``update_input_geos.sh`` which prepare the GEOS-Chem input files for the next run. The removal of signal files like ``ALL_RUNS_COMPLETE`` indicate to the ensemble run script that it is safe to start GEOS-Chem again.


update_input_geos.sh
~~~~~~~~~~~~~

A shell script which (1) calls ``advance_timestep.py`` to update the internal time stored in the ``scratch/`` directory, and (2) uses that updated internal time to update the ``input.geos`` file across the ensemble.

Assimilation support scripts
-------------

combine_columns_and_update.py
~~~~~~~~~~~~~

If the script ``check_and_complete_assimilation.sh`` finds that all expected ``.npy`` files containing assimilated columns are present in ``scratch/``, then this Python script is called. This script gathers the assimilated columns and loads in all the ensemble restarts and scaling factors, uses the contents of the columns to update restarts and scaling factors, and then writes the updated data to disk.

par_assim.sh
~~~~~~~~~~~~~

A wrapper shell script that calls ``par_letkf.py`` within the appropriate conda environment, passes information to the Python script ensuring that the appropriate set of columns are assimilated, and logs errors that occur in the assimilation process.

par_letkf.py
~~~~~~~~~~~~~

A short Python script, many instantiations of which are run in parallel, that creates relevant objects and calls methods from ``Assimilator.py`` to assimilate the set of columns assigned to a particular core or set of cores.

toolbox.py
~~~~~~~~~~~~~

Basic mathematical tools and utilities that are used across CHEEREIO Python scripts, including distance calculations, indexing support, and prior error covariance sampling. 

settings_interface.py
~~~~~~~~~~~~~

Basic utilities that interact with user settings and other global parameters and pass them to other Python scripts. 

LETKF classes
-------------

Assimilator.py
~~~~~~~~~~~~~

Contains the Assimilator class, which actually performs the LETKF operation. More details in the :ref:`Assimilator` entry.

GC_Translator.py
~~~~~~~~~~~~~

Contains the GC_Translator class and a few other support classes, which wraps around GEOS-Chem restarts and scaling factors and translates them into state vectors for use in CHEEREIO, and vice versa. More details in the :ref:`GC Translator` entry.

GT_Container.py
~~~~~~~~~~~~~

Contains the GT_Container class, which is used to combine assimilated columns and update GEOS-Chem after the LETKF operations complete. More details in the :ref:`GT Container` entry.

HIST_Translator.py
~~~~~~~~~~~~~

Contains the HIST_Translator class, which wraps around GEOS-Chem history output (concentrations saved out over time) and passes them on to the HIST_Ensemble object. More details in the :ref:`HIST Translator` entry.

HIST_Ens.py
~~~~~~~~~~~~~

Contains the HIST_Ensemble class, which combines GEOS-Chem history data from HIST_Translator objects and passes it to observation operators to create vectors of simulated observations for use in the LETKF algorithm. More details in the :ref:`HIST Ensemble` entry.


Observation operators
-------------

observation_operators.py
~~~~~~~~~~~~~

This Python file contains tools used to create observation operators. It also includes the parent class for all observation operators. Use of this file is described in detail in the :ref:`New observation` entry.

omi_tools.py
~~~~~~~~~~~~~

This Python file includes tools and classes necessary for interfacing with OMI satellite products, and is described in detail in the :ref:`OMI tools` entry.

iasi_tools.py
~~~~~~~~~~~~~

This Python file includes tools and classes necessary for interfacing with IASI satellite products, and is described in detail in the :ref:`iasi` entry.

tccon_tools.py
~~~~~~~~~~~~~

This Python file includes tools and classes necessary for interfacing with TCCON products, and is described in detail in the :ref:`tccon` entry.

prep_tccon_aggregated.py
~~~~~~~~~~~~~

This Python file produces TCCON aggregated files ready for CHEEREIO, and is described in detail in the :ref:`tccon` entry.

obspack_tools.py
~~~~~~~~~~~~~

This Python file includes tools and classes necessary for interfacing with NOAA ObsPack products, and is described in detail in the :ref:`ObsPack tools` entry.

preprocess_obspack.py
~~~~~~~~~~~~~

This Python file produces ObsPack aggregated files ready for CHEEREIO/GEOS-Chem, and is described in detail in the :ref:`ObsPack tools` entry.

tropomi_tools.py
~~~~~~~~~~~~~

This Python file includes tools and classes necessary for interfacing with TROPOMI satellite products, and is described in detail in the :ref:`TROPOMI tools` entry.

Utilities for the user
-------------

testing_tools.py
~~~~~~~~~~~~~

A suite of utilities used by CHEEREIO's ``pytest`` suite, as well as some utility functions for generating some assimilation objects for debugging CHEEREIO. See commentary in the script for details.

testing_workflow.py
~~~~~~~~~~~~~

A convenient script for diagnosing errors in assimilation. It takes command line arguments and uses them to test building an Assimilator object and walking through one assimilation step, with lots of print out for debugging. Essentially a convenience wrapper for part of testing_tools.py.

cleanup_after_kill_ens.sh
~~~~~~~~~~~~~

If the ensemble fails at runtime for a relatively simple reason, like a cluster hiccup or a minor bug, then you can use this script to clean up the ensemble and prepare it for resubmission. See the :ref:`Fix Kill Ens` entry for more information.

Deprecated scripts
-------------

The following scripts have been deprecated and will be removed:

* regrid_landmask_fraction.py

