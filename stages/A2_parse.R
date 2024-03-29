# 24/04/2015
# Team PREDICTS-PD
# Read in raw trees, rate-smooth and output

# START
cat (paste0 ('\nPub parse started at [', Sys.time (), ']\n'))

# LIBS
source (file.path ('tools', 'tree_tools.R'))

# DIRS
input.dir <- file.path ('0_data', 'raw_trees')
output.dir <- 'A2_parse'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

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
cat ('\nRate-smoothing where needed ....')
for (i in 1:length (trees)) {
  trees[[i]] <- runRateSmoother (trees[[i]], i)
}
cat ('\nDone.')

# WRITE OUT
for (i in 1:length (trees)) {
  if (length (trees[[i]]) > 1) {
    outfile <- file.path (output.dir, names (trees)[i])
    write.tree (phy=trees[[i]], file=outfile)
  }
}

# WRITE OUT W/ AGES
cat ('\n Adding node ages ....')
for (i in 1:length (trees)) {
  cat('\n .... [', names (trees)[i], ']', sep='')
  trees[[i]] <- addAges (phylos=trees[[i]])
}
pub.trees <- trees
save (pub.trees, file=file.path (output.dir, 'trees_with_ages.RData'))
cat('\nDone.')

# FINISH
cat (paste0 ('\nA2 finished at [', Sys.time (), ']'))
