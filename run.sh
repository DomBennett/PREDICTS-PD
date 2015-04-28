#!/bin/bash
# 01/11/2014
# Team PREDICTS-PD
# Run and log all stages -- except stages 2 and 8
# Usage: sh run.sh

echo 'Running pipeline'
Rscript stages/pGltsetup.R >& 1_pGltsetup_log.txt
Rscript stages/map.R >& 3_map_log.txt
Rscript stages/parse.R >& 4_parse_log.txt
Rscript stages/compare.R >& 5_compare_log.txt
Rscript stages/metrics.R >& 6_metrics_log.txt
Rscript stages/commplots.R >& 7_commplots_log.txt
echo 'Complete'
