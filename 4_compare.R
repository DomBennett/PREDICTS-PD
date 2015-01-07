# 01/11/2014
# Team PREDICTS-PD
# Compare pG-lt and published trees.

# PARAMETERS
min.tree <- 5  # minimum number of tips in a tree for reference

# START
cat (paste0 ('\nStage 4 started at [', Sys.time (), ']'))

# LIBS
library (ape)
source (file.path ('functions', 'tools.R'))

# DIRS
input.dirs <- c ('0_data', '3_parse')
output.dir <- '4_compare'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in trees ....')
# read in pub trees
pub.trees <- readInTrees (folder = file.path (input.dirs[1], 'pub_phylos'))
# read in pG-lt trees
pglt.trees <- readInTrees (folder = input.dirs[2])
cat ('\nDone.')

# PROCESS
cat ('\nCalculating shared nodes ....')
# calculate the proportion of nodes in the pglt tree that differ from the best pub tree
different.nodes <- tot.nodes <- rep (NA, length (pglt.trees))
for (i in 1:length (pglt.trees)) {
  cat ('\n.... tree [', i, '/', length (pglt.trees), ']', sep = '')
  # get dist
  treedist <- pglt.trees[[i]]
  tip.labels <- getNames (treedist)
  # find best reference tree
  ref.tree <- findBestRef (tip.labels, pub.trees)
  # for each in dist, break down to same sized tree and calc topo dist
  temp.tot.nodes <- temp.different.nodes <- rep (NA, length (treedist))
  for (j in 1:length (treedist)) {
    tree <- treedist[[j]]
    if (sum (tree$tip.label %in% ref.tree$tip.label) >= min.tree) {
      tree1 <- drop.tip (ref.tree, tip = ref.tree$tip.label[!ref.tree$tip.label %in% tree$tip.label])
      tree2 <- drop.tip (tree, tip = tree$tip.label[!tree$tip.label %in% tree1$tip.label])
      temp.different.nodes[j] <- dist.topo (tree1, tree2)
      temp.tot.nodes[j] <- (tree1$Nnode + tree2$Nnode)/2
    } else {
      temp.different.nodes[j] <- NA
      temp.tot.nodes[j] <- NA
    }
  }
  different.nodes[i] <- mean (temp.different.nodes, na.rm = TRUE)
  tot.nodes[i] <- mean (temp.tot.nodes, na.rm = TRUE)
}
p.shared.nodes <- 1 - mean (different.nodes/tot.nodes, na.rm = TRUE)
p.shared.nodes.sd <- sd (different.nodes/tot.nodes, na.rm = TRUE)
cat ('\nDone. [', p.shared.nodes, '±', p.shared.nodes.sd, '] shared nodes.', sep = '')
# plot (tot.nodes, different.nodes)

# TODO: do same with RF dist:
#  -- essentially the same as topo.dist but each node difference is also weighted by its branch -- should give a better result
#  -- requires both the reference and the subject trees to be ultramteric

# FINISH
cat (paste0 ('\nStage 4 finished at [', Sys.time (), ']'))