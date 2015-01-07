# 01/11/2014
# Team PREDICTS-PD
# Read in pG-lt trees and make ultrametric

# START
cat (paste0 ('\nStage 3 started at [', Sys.time (), ']'))

# LIBS
library (ape)

# DIRS
input.dir <- '2_pGltRun'
output.dir <- '3_parse'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading trees for each study ....')
tree.counter <- 0
studies <- list.files (input.dir)
trees <- list ()
for (study in studies) {
  fpath <- file.path (input.dir, study, '4_phylogeny', 'distribution.tre')
  treedist <- suppressWarnings (try (read.tree (fpath), silent = TRUE))
  if (class (treedist) != 'try-error') {
    trees <- c (trees, list (treedist))
    names (trees)[length (trees)] <- study
    tree.counter <- tree.counter + 1
  }
}
cat ('\nDone.')

# MAKE ULTRAMETRIC
cat ('\nRate-smoothing ....')
smooth.counter <- 0
# TODO: chronos or chrnonsML -- which is better?
cat ('\nDone.')

# OUTPUT
cat ('\nWriting out ....')
for (i in 1:length (trees)) {
  treedist <- trees[[i]]
  filename <- paste0 (names (trees)[i], '.tre')
  write.tree (file = file.path (output.dir, filename), phy = treedist)
}
cat ('Done. Discovered [', length (studies),
     '] studies of which [', tree.counter,
     '] had tree distirbutions that were read in and [',
     smooth.counter, '] were rate-smoothed.',
     sep='')

# FINISH
cat (paste0 ('\nStage 3 finished at [', Sys.time (), ']'))