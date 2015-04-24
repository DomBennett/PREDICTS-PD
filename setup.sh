#!/bin/bash
# 24/04/2014
# Team PREDICTS-PD
# Run all setup R scripts for R pipeline
# Usage: sh setup.sh

echo 'Running pipeline setup'
Rscript stages/setup_preresolve.R >& setup_preresolve_log.txt
Rscript stages/setup_parse.R >& setup_parse_log.txt
echo 'Complete'
