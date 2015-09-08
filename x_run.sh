#!/bin/bash
# 01/11/2014
# Team PREDICTS-PD
# Run and log all stages
# Usage: sh x_run.sh &

echo 'Running pipeline'
echo '[Make sure 0_pglt exists]'
sh x_setup.sh > setup_log.txt 2>&1
Rscript x_calculate.R > calculate_log.txt 2>&1
echo 'Complete'
