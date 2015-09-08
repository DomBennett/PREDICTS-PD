# 08/04/2015
# Team PREDICTS-PD
# Map PREDICTS names onto published trees

# START
cat (paste0 ('\nMap started at [', Sys.time (), ']\n'))

# LIBS (UNIX ONLY)
library (foreach)
library (doMC)
source (file.path ('tools', 'tree_tools.R'))

# PARAMETERS
ncpus <- 8
registerDoMC (ncpus)

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
studies <- list.files (input.dir)
counter <- foreach (i=1:length (studies)) %dopar% {
  counter <- NULL  # record m and n
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
      counter <- c (counter, 'm', rep ('n', length (names)))
    }
  }
  counter
}
cat ('\nDone.')

# OUTPUT
counter <- table (unlist (counter))
map.counter <- ifelse (is.na (counter['m']), 0, counter['m'])
names.counter <- ifelse (is.na (counter['n']), 0, counter['n'])
cat ('Done. Identified [', length (studies), '] studies of which [', map.counter,
     '] had names mapped to published trees representing a total [', names.counter ,
     '] names.', sep='')

# FINISH
cat (paste0 ('\nStage 3 finished at [', Sys.time (), ']'))