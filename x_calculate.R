# 01/11/2014
# Team PREDICTS-PD
# Run calculate
# Usage: RScript x_calculate.R

# Check pglt res
if (!file.exists ('0_pglt')) {
  stop ('No 0_pglt/ found')
}

# Parameters
ncpus <- 8
use.unconstrained.vec <- c (TRUE, FALSE)
normalise.vec <- c ('age', 'total')

# Pipeline
cat ('\nRunning pipeline ....')
# loop through compare stages with constrained and unconstrained
for (use.unconstrained in use.unconstrained.vec) {
  source ('stages/B1_parse.R', print.eval=TRUE, local=TRUE)
  source ('stages/B2_compare.R', print.eval=TRUE, local=TRUE)
}
# loop metrics for the two ways of normalising branch lengths
for (normalise in normalise.vec) {
  source ('stages/B3_metrics.R', print.eval=TRUE, local=TRUE)
}
source ('stages/B4_plots.R', print.eval=TRUE, local=TRUE)
cat ('\nComplete....')