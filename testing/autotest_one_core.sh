#!/bin/bash

cd ../core
bash par_assim.sh true 1 1
cd ../testing
bash diffcolumn.sh 0 0
bash diffcolumn.sh