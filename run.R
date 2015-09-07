# 01/11/2014
# Team PREDICTS-PD
# Run pipeline (e.g. RScript run.R > run_log.txt 2>&1)

# Parameters
ncpus <- 8
parameters <- list ('normalise'=c ('age', 'total'),
                    'use.unconstrained'=c (TRUE, FALSE))

# Pipeline
cat ('\nRunning pipeline ....')
source ('stages/pGltsetup.R', print.eval=TRUE, local=TRUE)
# loop through compare stages with constrained and unconstrained
for (use.unconstrained in parameters[['use.unconstrained']]) {
  source ('stages/parse.R', print.eval=TRUE, local=TRUE)
  source ('stages/compare.R', print.eval=TRUE, local=TRUE)
}
# loop metrics for the two ways of normalising branch lengths
for (normalise in parameters[['normalise']]) {
  source ('stages/metrics.R', print.eval=TRUE, local=TRUE)
}
source ('stages/commplots.R', print.eval=TRUE, local=TRUE)
cat ('\nComplete....')