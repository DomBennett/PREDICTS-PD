# 24/04/2015
# Team PREDICTS-PD
# Resolve ALL names in published trees and PREDICTS data

# START
cat (paste0 ('\nPre-resolve started at [', Sys.time (), ']\n'))

# LIBRARY
source (file.path ('tools', 'tree_tools.R'))

# DIRS
data.dir <- '0_data'
tree.dir <- file.path (data.dir, 'raw_trees')
predicts.dir <- file.path (data.dir, 'PREDICTS-DATA')
output.dir <- 'A1_preresolve'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# GET NAMES
cat ('\nReading in names ....')
# tree names
trees <- readInTrees (tree.dir)
names <- c ()
for (i in 1:length (trees)) {
  tree <- trees[[i]]
  if (class (tree) == 'multiPhylo') {
    tree <- tree[[1]]
  }
  names <- append (names, tree$tip.label)
}
# PREDICTS names
x <- readRDS (file.path (predicts.dir,
                         'diversity-2014-10-29-03-40-20.rds'))
names <- append (names, x$Parsed_name)
# clean up
names <- unique (names)
rm (trees, x)
cat('\nDone. Read in [', length (names),
    '] unique names', sep='')

# RESOLVE
cat ('\nResolving ....')
resolve.list <- mapNamesPreDownload (names)
cat ('\nDone.')

# SAVE
outfile <- file.path (output.dir, 'preresolved.RData')
save (resolve.list, file=outfile)

# FINISH
cat (paste0 ('\nFinished at [', Sys.time (), ']'))