#!/bin/bash
# 01/11/2014
# Team PREDICTS-PD
# Run and log all stages -- except stages 2 and 8
# Usage: sh run.sh > run_log.txt 2>&1 &

echo 'Running pipeline'
Rscript stages/pGltsetup.R > 1_pGltsetup_log.txt 2>&1
Rscript stages/map.R > 3_map_log.txt 2>&1
#Rscript stages/parse.R 2>& 4_parse_log.txt
#Rscript stages/compare.R 2>& 5_compare_log.txt
#Rscript stages/metrics.R 2>& 6_metrics_log.txt
#Rscript stages/commplots.R 2>& 7_commplots_log.txt
echo 'Complete'
