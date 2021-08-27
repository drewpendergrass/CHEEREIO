#!/bin/bash

sbatch --array={START}-{END} -W run_ensemble_simulations.sh

exit 0