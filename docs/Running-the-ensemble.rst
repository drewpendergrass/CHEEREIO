Running the ensemble
==========

Starting the run
-------------

Once the ensemble is fully installed, initiating the run is simple. Navigate to the Ensemble Runs Directory (see the :ref:`Guide to the Ensemble Directory` entry for details on ensemble directory structure). Then, execute the ``run_ens.sh`` script with the command ``nohup bash run_ens.sh &``. That's it!

Monitoring the run progress
-------------

Information about the ensemble state is continuously recorded during run time, and is stored in several locations.

* Overall run status for the entire job array (one job per ensemble member) is available from the SLURM scheduler. A good command is ``sacct``, which will display run status (it is especially useful as jobs are pending while resources become available). Each run of GEOS-Chem will be recorded as a separate entry under the ``time`` sub-job label, as each GEOS-Chem run (initialized by the outcome of the previous assimilation step) is submitted with a separate (timed) ``srun`` command. 
* GEOS-Chem run status for individual ensemble members are available in the ``GC.log`` file in each ensemble member run directory.
* Additional log files, including shell-level ``.err`` and ``.out`` files and log files containing data about assimilation, are all available in the ``log`` folder under :ref:`Ensemble Runs`.

A taxonomy of scripts called during runtime
-------------

TKTKTK

.. _Run Ensemble Simulations:

About the Run Ensemble Simulations script
-------------

TKTKTK

