# 24/04/2015
# Team PREDICTS-PD
# Read in raw trees, rate-smooth and output

# START
cat (paste0 ('\nSetup parse started at [', Sys.time (), ']\n'))

# LIBS
library (ape)

# DIRS
input.dir <- file.path ('0_data', 'raw_trees')
output.dir <- file.path ('0_data', 'parsed_trees')

# PROCESS
cat('\nReading in trees ....')
treefiles <- list.files (input.dir, pattern = '.tre')
trees <- list ()
for (i in 1:length (treefiles)) {
  cat ('\n........ working on [', i, '/', length (treefiles), ']', sep='')
  tree <- read.tree (file.path (input.dir, treefiles[i]))
  if (class (tree) == 'multiPhylo') {
    # if multiphylo, select a random 100
    random <- sample (1:length (tree), 100)
    tree <- tree[random]
  }
  # add to list
  trees[treefiles[i]] <- list (tree)
}
cat('\nDone.')

# RATE-SMOOTH
# TODO

# WRITE OUT
for (i in 1:length (trees)) {
  outfile <- file.path (output.dir, names (trees)[i])
  write.tree (phy=trees[[i]], file=outfile)
}

# FINISH
cat (paste0 ('\nFinished at [', Sys.time (), ']'))