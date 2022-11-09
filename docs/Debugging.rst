Debugging CHEEREIO
==========

The testing suite
-------------

CHEEREIO comes with a testing suite that can be run using the pytest package. This testing workflow verifies that the most important aspects of the assimilation module are working properly. For example, it checks that state vectors are formed and subsetted properly, and that assimilation calculations are correct. It also verifies that satellite operators are working as expected and that data is being passed across modules correctly.

However, not every scenario can be tested by this workflow, so please be careful! This workflow also does not test the full CHEEREIO run, which involves running GEOS-Chem. It only works on assimilation. Run this workflow every time you make modifications to the assimilation workflow, just to make sure you didn't break it!

To ensure that testing results are reproducible in different computing environments, CHEEREIO includes data for testing in the ``data_for_tests/`` folder. CHEEREIO also will temporarily override the existing ``ens_config.json`` file with standardized settings for test stability. Within the ``data_for_tests/`` folder, there is are several files ending with ``_settings_to_override.json`` that list the settings that will be used to override ``ens_config.json`` for various test types. CHEEREIO comes with a cleanup routine to ensure that the code directory is returned to its previous state after testing is complete.

To execute the testing suite, navigate to the ``tests`` folder in the CHEEREIO code directory and activate your CHEEREIO conda environment. Then, run the command 'pytest' at the command line within this directory. Test results will be printed to the terminal.

To add new tests to the testing suite, add functions with no input that with names that follow the format ``test_*()``. You can make use of the ``testing_tools.py`` module within the ``core/`` folder to ensure stable testing environments and to generate useful data structures. Every testing function should end with an ``assert`` command that takes a boolean. See the `Pytest documentation <https://docs.pytest.org/en/7.1.x/contents.html>`__ for more information.


Help! CHEEREIO killed my ensemble
-------------

So you've found the dreaded KILL_ENS file in your ensemble's ``scratch`` folder. This file is produced by CHEEREIO if something goes wrong with either (1) GEOS-Chem, or (2) the assimilation routine. The presence of this file signals all runs to cease operations and exit with nonzero status. Here we discuss the most common causes of ensemble failure and how to fix them. First, if the ensemble is still running, cancel it. Next, examine the contents of KILL_ENS. The file will give you an idea of where to look for more detailed error information. 

.. image:: kill_ens.png
  :width: 600
  :alt: Picture of terminal with KILL_ENS file present. 


GEOS-Chem problems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you examine the contents of KILL_ENS, you might find the message "GEOS-Chem in ensemble member x did not complete successfully." This means that at least one of our GEOS-Chem runs failed to complete successfully. Navigate to the corresponding run directory in the ``ensemble_runs`` folder and look at the end of the GEOS-Chem log file ``GC.log``. This file will have a traceback that describes the GEOS-Chem error. You should also look at the corresponding shell error file, which will be stored in the ``ensemble_runs/logs`` folder in a file named ``ensemble_slurm_JOBNUMBER.err``, which might offer further information like an out-of-memory error. Consult the `GEOS-Chem debugging page <http://wiki.seas.harvard.edu/geos-chem/index.php/Debugging_GEOS-Chem>`__ for more guidance on how to debug GEOS-Chem. 

Some rules of thumb: If GEOS-Chem fails immediately after the runs begin, the problem probably has to do with either the `Unix environment <https://geos-chem.readthedocs.io/en/latest/gcc-guide/01-startup/login-env.html#>`__ or a problem with settings files, such as HEMCO_Config.rc or input.geos, or a missing or incompatible data files (e.g. no restart, no meteorology). The GC.log file will give you a clue as to what your problem is. If you are missing a library file (e.g. libnetcdf.so) this means your Unix environment is not correctly specified. We recommend that you modify the ``cheereio.env`` file provided in the ``environments`` folder to match your cluster environment. It is also very common in CHEEREIO to hit out-of-memory errors, though this is more often seen in the assimilation phase. Finally, like any piece of complicated software in a shared environment, sometimes GEOS-Chem fails randomly. 

Assimilation problems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you examine the contents of KILL_ENS, you might find the message "Python assimilation script exited without code 0 in ensemble x and core y." This means that at least one of our assimilation routines failed to complete successfully. To figure out what went wrong, look at the corresponding shell error file, which will be stored in the ``ensemble_runs/logs`` folder in a file named ``ensemble_slurm_JOBNUMBER.err``. This will include the Python traceback that caused the problem. For fully-tested CHEEREIO routines, such as TROPOMI CH4, the most common cause of problems are out-of-memory errors. Occasionally assimilation will also fail for random reasons related to cluster hiccups. For less tested CHEEREIO experiments, the error may be more substantial and emerge from either errors in user implementation, or from a bug in the software. Please open an `issue <https://github.com/drewpendergrass/CHEEREIO/issues>`__ on Github if you think that your error is due to the base CHEEREIO software, or if you have a question best addressed by the developer. 

I think I solved it!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you conclude that the KILL_ENS file was produced by a minor problem, such as an out-of-memory error, bad Unix environment, or an easily-resolved fix to source code, you can run the ``cleanup_after_kill_ens.sh`` script from within the ``core/`` folder. This script will take care of some of the main technicalities that can cause a resubmitted CHEEREIO run to fail, such as making sure that all signal files are aligned. From there, you can simply resubmit the job array. However, CHEEREIO is a complex piece of code and sometimes errors can root themselves deeply in the ensemble in ways that are difficult to remove. If your second submission fails, it is easiest to (1) use the clean spun-up backup ensemble, or (2) fully redeploy the ensemble. This is the CHEEREIO equivalent of turning it off and back on again.



