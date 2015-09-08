#!/bin/bash
# 01/11/2014
# Team PREDICTS-PD
# Run and log all stages
# Usage: sh run.sh &

echo 'Running pipeline'
sh setup.sh > setup_log.txt 2>&1
Rscript calculate.R > calculate_log.txt 2>&1
echo 'Complete'
