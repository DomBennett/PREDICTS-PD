# 08/04/2015
# Team PREDICTS-PD
# Map PREDICTS names onto published trees

# START
cat (paste0 ('\nStage 3 started at [', Sys.time (), ']'))

# LIBS
source (file.path ('tools', 'tree_tools.R'))

# DIRS
data.dir <- '0_data'
input.dir <- '1_pGltsetup'
output.dir <- '3_map'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading published trees ....')
pub.trees <- readInTrees (folder=file.path (data.dir, 'pub_phylos'))
# add ages to pub trees to speed up mapNames
for (i in 1:length (pub.trees)) {
  cat('\n .... [', names (pub.trees)[i], ']', sep='')
  # ensure no node labels
  pub.trees[[i]]$node.label <- NULL
  if (!is.null (pub.trees[[i]]$edge.length) && is.ultrametric (pub.trees[[i]])) {
    cat ('\n ........ adding node ages')
    pub.trees[[i]]$node.ages <- getAge (tree=pub.trees[[i]])[ ,2]
  }
}
cat('\nDone.')
# create subject environment -- holds name resolutions of subject names
# for mapNames. This will prevent searching same names multiple times
sbjctenv <- new.env (parent=emptyenv ())
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
    mapped.trees <- mapNames (tree=best.tree[[1]], names=names,
                              iterations=100)
    if (!is.null (mapped.trees)) {
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