#!/bin/bash

#Usage: bash snapshot.sh [--gossipgirl]
#Optional gossip girl flag changes print output.

source ../environments/cheereio.env #This is specific to the Harvard cluster; rewrite for yours
eval "$(conda shell.bash hook)"

conda activate $(jq -r ".CondaEnv" ../ens_config.json) 
python snapshot_pt1.py $1
conda deactivate 
conda activate $(jq -r ".AnimationEnv" ../ens_config.json) 
python snapshot_pt2.py $1
conda deactivate 


