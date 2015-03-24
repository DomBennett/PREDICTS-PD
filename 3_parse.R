# 01/11/2014
# Team PREDICTS-PD
# Parse published and pG-lt trees

# START
cat (paste0 ('\nStage 3 started at [', Sys.time (), ']'))

# LIBS
source (file.path ('tools', 'tree_tools.R'))

# DIRS
data.dir <- '0_data'
input.dir <- '2_pGltRun'
output.dir <- '3_parse'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading trees for each study ....')
pub.trees <- readInTrees (folder=file.path (data.dir, 'pub_phylos'))
# add ages to pub trees to speed up mapNames
for (i in 1:length (pub.trees)) {
  if (!is.null (pub.trees[[i]]$edge.length) && is.ultrametric (pub.trees[[i]])) {
    pub.trees[[i]]$node.ages <- getAge (tree=pub.trees[[i]])[ ,2]
    # ensure no node labels
    pub.trees[[i]]$node.label <- NULL
  }
}
# create subject environment -- holds name resolutions of subject names
# for mapNames. This will prevent searching same names multiple times
sbjctenv <- new.env (parent=emptyenv ())
# counters
pglt.counter <- pub.counter <- 0
studies <- list.files (input.dir)
trees <- list ()
for (study in studies) {
  # init container
  study.tree <- list ()
  # pglt-trees
  fpath <- file.path (input.dir, study, '4_phylogeny', 'distribution.tre')
  pglt.trees <- suppressWarnings (try (read.tree (fpath), silent = TRUE))
  if (class (pglt.trees) != 'try-error') {
    study.tree['pglt.trees'] <- list (pglt.trees)
    pglt.counter <- pglt.counter + 1
  }
  # published trees
  names <- read.delim (file.path (input.dir, study, 'names.txt'), header=TRUE,
                       stringsAsFactors=FALSE)[ ,1]
  best.tree <- findBestRef (tip.labels=names, ref.trees=pub.trees)
  if (!is.null (best.tree[[1]]) && !is.null (best.tree[[1]]$edge.length) && is.ultrametric (best.tree[[1]])) {
    # generate a distribution of trees from random name mapping
    mapped.trees <- mapNames (tree=best.tree[[1]], names=names,
                              iterations=1)
    if (!is.null (mapped.trees)) {
      study.tree['mapped.trees'] <- list (mapped.trees)
      pub.counter <- pub.counter + 1
    }
  }
  if (length (study.tree) > 0) {
    # add study container to all container
    trees[study] <- list (study.tree)
  }
}
cat ('\nDone.')

# MAKE ULTRAMETRIC
cat ('\nRate-smoothing ....')
smooth.counter <- 0
# TODO: chronos or chrnonsML -- which is better?
cat ('\nDone.')

studies

# OUTPUT
cat ('\nWriting out ....')
for (i in 1:length (trees)) {
  pglt.trees <- trees[[i]][['pglt.trees']]
  mapped.trees <- trees[[i]][['mapped.trees']]
  folder.path <- file.path (output.dir, paste0 (names (trees)[i]))
  if (!file.exists (folder.path)) {
    dir.create (folder.path)
  }
  if (!is.null (mapped.trees)) {
    write.tree (file = file.path (folder.path, 'mapped.tre'), phy = mapped.trees)
  }
  write.tree (file = file.path (folder.path, 'pglt.tre'), phy = pglt.trees)
}
cat ('Done. Discovered [', length (studies),
     '] studies of which [', pglt.counter,
     '] had pglt distirbutions and [', pub.counter ,
     '] had names mapped to published trees that were read in and [',
     smooth.counter, '] were rate-smoothed.',
     sep='')

# FINISH
cat (paste0 ('\nStage 3 finished at [', Sys.time (), ']'))