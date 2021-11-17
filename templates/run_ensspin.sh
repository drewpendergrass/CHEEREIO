#!/bin/bash

sbatch --array={START}-{END} -W run_ensemble_spinup_simulations.sh

exit 0