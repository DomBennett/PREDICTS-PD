# 01/11/2014
# Team PREDICTS-PD
# Run calculate
# Usage: RScript run.R

# Parameters
ncpus <- 8
use.unconstrained.vec <- c (TRUE, FALSE)
normalise.vec <- c ('age', 'total')

# Pipeline
cat ('\nRunning pipeline ....')
source ('stages/pGltsetup.R', print.eval=TRUE, local=TRUE)
# loop through compare stages with constrained and unconstrained
for (use.unconstrained in use.unconstrained.vec) {
  source ('stages/parse.R', print.eval=TRUE, local=TRUE)
  source ('stages/compare.R', print.eval=TRUE, local=TRUE)
}
# loop metrics for the two ways of normalising branch lengths
for (normalise in normalise.vec) {
  source ('stages/metrics.R', print.eval=TRUE, local=TRUE)
}
source ('stages/commplots.R', print.eval=TRUE, local=TRUE)
cat ('\nComplete....')