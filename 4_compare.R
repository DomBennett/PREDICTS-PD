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
# calculate the proportion of nodes/branch in the pglt tree that differ from the best pub tree
different.branch <- different.nodes <- tot.nodes <- ref.trees <- rep (NA, length (pglt.trees))
for (i in 1:length (pglt.trees)) {
  cat ('\n.... tree [', i, '/', length (pglt.trees), ']', sep = '')
  # get dist
  treedist <- pglt.trees[[i]]
  tip.labels <- getNames (treedist)
  # find best reference tree
  ref.tree <- findBestRef (tip.labels, pub.trees)
  ref.trees[i] <- names (ref.tree)[1]
  ref.tree <- ref.tree[[1]]
  # for each in dist, break down to same sized tree and calc topo dist
  temp.tot.nodes <- temp.different.nodes <- temp.different.branch <-
    rep (NA, length (treedist))
  for (j in 1:length (treedist)) {
    tree <- treedist[[j]]
    if (sum (tree$tip.label %in% ref.tree$tip.label) >= min.tree) {
      # break down to comparable trees
      tree1 <- drop.tip (ref.tree, tip = ref.tree$tip.label[!ref.tree$tip.label %in% tree$tip.label])
      tree2 <- drop.tip (tree, tip = tree$tip.label[!tree$tip.label %in% tree1$tip.label])
      # calculate topo.dist
      temp.different.nodes[j] <- dist.topo (tree1, tree2)
      temp.tot.nodes[j] <- (tree1$Nnode + tree2$Nnode)/2
      # calculate branch dist
      if (!is.null (tree1$edge.length)) {
        # scale both to 1
        tree1$edge.length <- tree1$edge.length/sum (tree1$edge.length)
        tree2$edge.length <- tree2$edge.length/sum (tree2$edge.length)
        temp.different.branch[j] <- dist.topo (tree1, tree2, 'score')
      }
    }
  }
  different.nodes[i] <- mean (temp.different.nodes, na.rm = TRUE)
  tot.nodes[i] <- mean (temp.tot.nodes, na.rm = TRUE)
  different.branch[i] <- mean (temp.different.branch, na.rm = TRUE)
}
p.different.branch <- mean (different.branch, na.rm = TRUE)
p.different.branch.sd <- sd (different.branch, na.rm = TRUE)
p.different.nodes <- mean (different.nodes/tot.nodes, na.rm = TRUE)
p.different.nodes.sd <- sd (different.nodes/tot.nodes, na.rm = TRUE)
cat ('\nDone. [', p.different.nodes, '±', p.different.nodes.sd, '] different nodes and [',
     p.different.branch, '±',  p.different.branch.sd, '] different branch.', sep = '')
# plot (tot.nodes, different.nodes)

# OUTPUT
pdf (file.path (output.dir, 'comparisons.pdf'))
pull <- !is.na (different.branch)
plot (different.branch[pull] ~ factor(ref.trees[pull]), xlab = 'Comparison Tree',
      ylab = 'Proportion of different branch')
pull <- !is.na (different.nodes)
plot ((different.nodes[pull]/tot.nodes[pull]) ~ factor(ref.trees[pull]), xlab = 'Comparison Tree',
      ylab = 'Proportion of different nodes')
dev.off ()
save (different.branch, different.nodes, tot.nodes, ref.trees,
      file = file.path (output.dir, 'results.RData'))

# FINISH
cat (paste0 ('\nStage 4 finished at [', Sys.time (), ']'))