# 08/04/2015
# Team PREDICTS-PD
# Map PREDICTS names onto published trees

# START
cat (paste0 ('\nStage 3 started at [', Sys.time (), ']\n'))

# LIBS
source (file.path ('tools', 'tree_tools.R'))

# DIRS
data.dir <- '0_data'
tree.dir <- file.path (data.dir, 'parsed_trees')
input.dir <- '1_pGltsetup'
output.dir <- '3_map'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading published trees ....')
# read in parsed and aged pub trees
load (file.path (tree.dir, 'trees_with_ages.RData'))
for (i in 1:length (pub.trees)) {
  if ('chronos' %in% class (pub.trees[[i]])) {
    class (pub.trees[[i]]) <- 'phylo'
  }
  pub.trees[[i]] <- dropUnderscore (pub.trees[[i]])
}
# read in preresolved names
load (file.path (tree.dir, 'preresolved.RData'))
cat ('\nMapping names to published trees and outputting for each study ....')
# counter and nnames
counter <- nnames <- 0
studies <- list.files (input.dir)
for (i in 1:length (studies)) {
  study <- studies[i]
  cat ('\n.... study [', study, '], [', i, '/', length (studies), ']',
       sep = '')
  # get names
  names <- read.delim (file.path (input.dir, study, 'names.txt'), header=TRUE,
                       stringsAsFactors=FALSE)[ ,1]
  # find best published tree using string matching
  best.tree <- findBestRef (tip.labels=names, ref.trees=pub.trees)
  if (!is.null (best.tree[[1]])) {
    cat ('\n........ mapping')
    # generate a distribution of trees from random name mapping
    mapped.trees <- NULL
    try ({
      mapped.trees <- mapNames (tree=best.tree[[1]], names=names,
                                iterations=100, resolve.list=resolve.list)
    }, silent=TRUE)
    if (!is.null (mapped.trees) & !is.na (mapped.trees[[1]])) {
      cat ('\n........ outputting')
      filename <- paste0 (file.path (output.dir, study), '.tre')
      write.tree (file=filename, phy=mapped.trees)
      counter <- counter + 1
      nnames <- nnames + length (names)
    }
  }
}
cat ('\nDone.')

# OUTPUT
cat ('Done. Identified [', length (studies), '] studies of which [', counter,
     '] had names mapped to published trees representing a total [', nnames ,
     '] names.', sep='')

# FINISH
cat (paste0 ('\nStage 3 finished at [', Sys.time (), ']'))