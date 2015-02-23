# 04/11/2014
# Team PREDICTS-PD
# Tree wrangling tools

# LIBS
library (MoreTreeTools)

# FUNCTIONS
readInTrees <- function (folder) {
  # point at a folder, it will return a list of trees in that folder
  filenames <- list.files (path = folder, pattern = '\\.tre')
  studies <- sub ('\\.tre', '', filenames)
  all.trees <- list ()
  for (i in 1:length (studies)) {
    treedist <- read.tree (file.path (folder, filenames[i]))
    # drop underscores from tip names
    treedist <- dropUnderscore(treedist)
    all.trees <- c (all.trees, list (treedist))
    names (all.trees)[i] <- studies[i]
  }
  all.trees
}

dropUnderscore <- function (phylos) {
  # drop _ in names of trees in a multiphylo
  if (class (phylos) == 'multiPhylo') {
    .drop <- function (i) {
      phylo <- phylos[[i]]
      phylo$tip.label <- gsub ('_', ' ', phylo$tip.label)
      res <<- c (res, list (phylo))
    }
    res <- list ()
    m_ply (.data=data.frame (i=1:length (phylos)), .fun=.drop)
    class (res) <- 'multiPhylo'
    return (res)
  } else {
    phylos$tip.label <- gsub ('_', ' ', phylos$tip.label)
    return (phylos)
  }
}

findBestRef <- function (tip.labels, ref.trees) {
  # find the best reference tree based on given tip.labels
  ptips <- rep (NA, length (ref.trees))
  for (i in 1:length (ref.trees)) {
    ptips[i] <- sum (ref.trees[[i]]$tip.label %in% tip.labels)/length (tip.labels)
  }
  besti <- which (ptips == max (ptips))[1]
  ref.trees[besti]
}

getNames <- function(phylos) {
  # get tip names from a multiphylo
  .get <- function (i) {
    res <<- c (res, phylos[[i]]$tip.label)
  }
  res <- NULL
  m_ply (.data=data.frame (i=1:length (phylos)), .fun=.get)
  unique (res)
}