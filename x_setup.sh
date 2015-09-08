#!/bin/bash
# 24/04/2014
# Team PREDICTS-PD
# Run all setup R scripts
# Usage: sh x_setup.sh

echo 'Running setup'
#Rscript stages/A1_preresolve.R
Rscript stages/A2_parse.R
Rscript stages/A3_pgltsetup.R
Rscript stages/A4_map.R
echo 'Complete'
echo '[Now generate 0_pglt with A3_pgltsetup]'
