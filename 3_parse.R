## 01/11/2014
## Team PREDICTS-PD
## Read in pG-lt trees and make ultrametric

# START
cat (paste0 ('\nStage 1 started at [', Sys.time (), ']'))

## Libs
library (ape)

## Dirs
input.dir <- '2_pGltRun'
output.dir <- '3_parse'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

## Input
# read in dist for each study
studies <- list.files (input.dir)
trees <- list ()
for (study in studies) {
  treedist <- read.tree (file.path (input.dir, study, '4_phylogeny', 'distribution.tre'))
  trees <- c (trees, list (treedist))
  names (trees)[length (trees)] <- study
}

## Make ultrametric (TODO)
# chronos or chrnonsML -- which is better?

## Output
for (i in 1:length (trees)) {
  treedist <- trees[[i]]
  filename <- paste0 (names (trees)[i], '.tre')
  write.tree (file = file.path (output.dir, filename), phy = treedist)
}