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

.. _Run Ensemble Simulations:

About the Run Ensemble Simulations script
-------------

Although the user will not ever execute the ``run_ensemble_simulations.sh`` script manually, it is very useful in terms of debugging to understand how the script works. We'll walk through it step-by-step.

SBATCH header
~~~~~~~~~~~~~

The first part of the ``run_ensemble_simulations.sh`` script should be familiar to anyone who has worked with the SLURM job scheduler. Entries in curly braces are replaced with user settings from ``ens_config.json`` at Template Run Directory creation time. In short, this part of the script tells the scheduler about the resources required by a single ensemble member, what partition the job should use, and names two files to store shell-level output.  

.. code-block:: bash

	#!/bin/bash

	#SBATCH -J {RunName}
	#SBATCH -c {NumCores}
	#SBATCH -N 1
	#SBATCH -p {Partition}
	#SBATCH --mem {Memory}
	#SBATCH -t {WallTime}
	#SBATCH -o logs/ensemble_slurm_%j.out    # File to which STDOUT will be written, %j inserts jobid       
	#SBATCH -e logs/ensemble_slurm_%j.err    # File to which STDERR will be written, %j inserts jobid

One-time initializations
~~~~~~~~~~~~~

Before GEOS-Chem is run or assimilation can be calculated, a few global settings have to be handled for the ensemble member. First, the software environment must be configured correctly because CHEEREIO requires many modules that can conflict with one another. This is accomplished with the following lines, which (1) purge and load appropriate modules, and (2) configures Anaconda for the subshell that is running the job.

.. code-block:: bash

	#Source clean environment with compatible netcdf and compiler environments and packages like GNU parallel:
	source {ASSIM}/environments/cheereio.env #This is specific to the Harvard cluster; rewrite for yours
	eval "$(conda shell.bash hook)"

Next, a few global variables are set. The end user will likely not ever make use of CHEEREIO's testing suite, so the variable ``TESTING`` is set to false automatically. A few key directories are stored in the variables ``$ENSDIR`` (the Ensemble Runs directory), ``$MY_PATH`` (path to the directory containing all CHEEREIO ensembles), and ``$RUN_NAME`` (the name of this ensemble). The latter two are grabbed from the ``ens_config.json`` file using the ``jq`` command, which allows shell scripts to access data stored in JSON format on the disk. The variable ``$x`` includes the ensemble member ID (ranging from 1 to the total number of ensemble members).  

.. code-block:: bash

	### Run directory
	TESTING={TESTBOOL}
	ENSDIR=$(pwd -P)

	if [ "${TESTING}" = true ]; then
	  MY_PATH="$(jq -r ".MY_PATH" {ASSIM}/testing/test_config.json)"
	  RUN_NAME="$(jq -r ".RUN_NAME" {ASSIM}/testing/test_config.json)"
	else
	  MY_PATH="$(jq -r ".MY_PATH" {ASSIM}/ens_config.json)"
	  RUN_NAME="$(jq -r ".RUN_NAME" {ASSIM}/ens_config.json)"
	fi

	### Get current task ID
	x=${SLURM_ARRAY_TASK_ID}

	### Add zeros to the current task ID
	if [ $x -lt 10 ]; then
	  xstr="000${x}"
	elif [ $x -lt 100 ]; then
	  xstr="00${x}"
	elif [ $x -lt 1000 ]; then
	  xstr="0${x}"
	else
	  xstr="${x}"
	fi

After all these variables are set, then CHEEREIO navigates to the particular ensemble run directory it has been assigned and exports the proper number of OpenMP threads so that GEOS-Chem can exploit parallelization.

.. code-block:: bash

	### Run GEOS-Chem in the directory corresponding to the cluster Id
	cd  {RunName}_${xstr}

	# Set the proper # of threads for OpenMP
	# SLURM_CPUS_PER_TASK ensures this matches the number you set with NumCores in the ens_config file
	export OMP_NUM_THREADS={NumCores}

While loop part 1: Run GEOS-Chem
~~~~~~~~~~~~~

With global settings taken care of, we can now proceed to the while loop that repeatedly runs GEOS-Chem and assimilation routines until the entire ensemble is completed. The whole process starts off with the following while statement:

.. code-block:: bash

	#Run GC; hang until assimilation complete (later also will do assimilation).
	#This will loop until a file appears in scratch signalling assimilation is complete.
	while [ ! -f ${MY_PATH}/${RUN_NAME}/scratch/ENSEMBLE_COMPLETE ]; do

This means that the while loop will continue until a file named ``ENSEMBLE_COMPLETE`` appears in the Scratch folder. The first thing that happens in the while loop is that GEOS-Chem is submitted. Again, since ``$TESTING`` will always be ``false`` for end users, we can consider only this part of the code block.

.. code-block:: bash

	  # Run GEOS_Chem.  The "time" command will return CPU and wall times.
	  # Stdout and stderr will be directed to the "GC.log" log file
	  # If just testing assimilation, skip all this
	  if [ "${TESTING}" = false ]; then
	    srun -c $OMP_NUM_THREADS time -p ./gcclassic >> GC.log
	    wait

This runs GEOS-Chem for one assimilation period (often just 24 hours). The ``wait`` command means that the job will hang until the ``srun`` job managing GEOS-Chem is complete. When ``srun`` terminates, CHEEREIO looks at the last line of the ``GC.log`` file and checks if it indicates that GEOS-Chem terminated successfully.

.. code-block:: bash

	    taillog="$(tail -n 1 GC.log)"
	    #Check if GC finished.
	    if [[ ${taillog:0:1} != "*" ]]; then
	      printf "GEOS-Chem did not complete successfully\n" > ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS #This file's presence will break loop
	    fi
	      #If there is a problem, the KILL_ENS file will be produced. Break then
	    if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
	      break
	    fi

If GEOS-Chem fails, then the file ``KILL_ENS`` is created and stored in the Scratch directory. Other ensemble members will detect this file and terminate themselves as well, because the ensemble can only continue if all GEOS-Chem runs terminate successfully. However, if everything works, the ensemble member 1 takes on the role as the "job coordinator." This ensemble member checks every second if a restart with the appropriate time stamp is present in each ensemble run directory. If everything is present, then a file labeled ``ALL_RUNS_COMPLETE`` is created and stored in the Scratch directory and we can continue to the assimilation phase (as it breaks the until loop). The rest of the ensemble members just check for the ``KILL_ENS`` file repeatedly and terminate themselves if necessary.

.. code-block:: bash

	    #Ensemble member 1 handles checking. CD to core.
	    if [ $x -eq 1 ]; then
	      cd {ASSIM}/core
	    fi
	    #Hang until ALL_RUNS_COMPLETE found in scratch folder
	    until [ -f ${MY_PATH}/${RUN_NAME}/scratch/ALL_RUNS_COMPLETE ]
	    do
	      #If this is ensemble member 1, look for all restarts and flag if found. Otherwise do nothing.
	      if [ $x -eq 1 ]; then
	        bash check_for_all_restarts.sh #Check if restarts exist; if they do, create ALL_RUNS_COMPLETE.
	      fi
	      #If there is a problem, the KILL_ENS file will be produced. Break then
	      if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
	        break 2
	      fi
	      sleep 1
	    done

If testing mode is active, GEOS-Chem is not run and we proceed to assimilation immediately. Again, this is not relevant to the end-user.

.. code-block:: bash

	  else
	    #Create done signal
	    if [ $x -eq 1 ]; then
	      echo "Done" > ${MY_PATH}/${RUN_NAME}/scratch/ALL_RUNS_COMPLETE
	    fi
	  fi

While loop part 2: Execute parallelized assimilation
~~~~~~~~~~~~~

With all GEOS-Chem runs completed successfully, we can now begin the assimilation process. All ensemble members navigate to the CHEEREIO code directory in order to submit the ``par_assim.sh`` script via GNU Parallel.

.. code-block:: bash

	  #CD to core
	  cd {ASSIM}/core
	  #Use GNU parallel to submit parallel sruns, except nature
	  if [ $x -ne 0 ]; then
	    parallel -j {MaxPar} "bash par_assim.sh ${TESTING} ${x} {1}" ::: {1..{NumCores}}
	  fi

The GNU parallel line works as follows. Up to ``MaxPar`` jobs in a single ensemble member will run the command ``bash par_assim.sh ${TESTING} ${x} {1}`` simultaneously. The ``par_assim.sh`` takes three command line inputs: a boolean signal for whether or not we are in testing mode; and ensemble ID number; and a core ID number. The first two inputs are supplied by global settings, while the third is supplied by a special GNU Parallel substitution line. Each core will then compute the LETKF data assimilation for each of its assigned columns and save them in ``.npy`` format to the scratch directory.

While loop part 3: Clean-up and ensemble completion 
~~~~~~~~~~~~~

Once this parallelized assimilation is complete, a fair amount of clean up must be done before the entire while loop can repeat. Before CHEEREIO can update the NetCDFs containing restarts and scaling factors in each ensemble member run directory, we have to wait for all columns to be saved to the Scratch directory. The loop thus hangs until a file labeled ``ASSIMILATION_COMPLETE`` appears in the Scratch directory. While we hang, the script ``check_and_complete_assimilation.sh`` is run every second by ensemble member 1. If the number of ``*.npy`` files in scratch matches the expected number of columns, then ensemble member 1 will load all the columns in and update the relevant NetCDFs (and create the ``ASSIMILATION_COMPLETE`` file).

.. code-block:: bash

	  #Hang until assimilation completes or cleanup completes (in case things go too quickly)
	  until [ -f ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE ] || [ ! -f ${MY_PATH}/${RUN_NAME}/scratch/ALL_RUNS_COMPLETE ]; do
	    #If this is ensemble member 1, check if assimilation is complete; if it is, do the final overwrites.
	    if [ $x -eq 1 ]; then
	      bash check_and_complete_assimilation.sh ${TESTING}
	    fi
	    sleep 1
	  done

Next, we do our usual check for the ``KILL_ENS`` file, which will be generated if any of the CHEEREIO Python scripts fail. If everything works correctly, then ensemble member 1 runs the script ``cleanup.sh``. This clean-up script deletes flag files and column files particular to this assimilation run from the Scratch folder, updates internal state files with the new correct dates, and prepares the ``input.geos`` file in each ensemble member run directory so that we can correctly run GEOS-Chem for the next assimilation period.

.. code-block:: bash

	  #If there is a problem, the KILL_ENS file will be produced. Break then
	  if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
	    break
	  fi
	  #If this is ensemble member 1, execute cleanup. This is because we only want it to run once.
	  if [ $x -eq 1 ]; then
	    bash cleanup.sh ${TESTING} #This also will break us out of this loop when assimilation complete.
	  fi
	  #Hang until cleanup complete, as determined by temp file deletion.
	  until [ ! -f ${MY_PATH}/${RUN_NAME}/scratch/ASSIMILATION_COMPLETE ]; do
	    sleep 1
	  done

	  #CD back to run directory
	  cd ${ENSDIR}/{RunName}_${xstr}
	  #Everything cleaned up; we can head back to the beginning.
	done

If, after clean-up completes, our current date now exceeds the ensemble end time stored in the ``ens_config.json`` file, then the file ``ENSEMBLE_COMPLETE`` will appear in the Scratch directory and the while loop will terminate. However, the loop will also terminate if the ``KILL_ENS`` file is created due to a script failure. CHEEREIO picks the exit code for the job accordingly.

.. code-block:: bash

	#If there is a problem, the KILL_ENS file will be produced. Exit with code 1 in that case
	if [ -f ${MY_PATH}/${RUN_NAME}/scratch/KILL_ENS ]; then
	  exit 1
	else
	  exit 0
	fi



