#!/bin/bash
# 24/04/2014
# Team PREDICTS-PD
# Run all setup R scripts
# Usage: sh setup.sh

echo 'Running setup'
#Rscript stages/setup_preresolve.R
Rscript stages/pub_parse.R
Rscript stages/map.R
echo 'Complete'
